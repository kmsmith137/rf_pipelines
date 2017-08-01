import zmq
import numpy as np
from scipy.ndimage import maximum_filter
from itertools import product
import L1_event

simple_dtype = [('itree', int), ('dm', int), ('sm', int),
                ('beta', int), ('time', int), ('snr', float)]


class L1Grouper(object):
    """
    A class for grouping coarse-grained triggers and forming L1 events

    Parameters
    ----------
    dedisp : instance of bonsai.Dedisperser
        Used primarily to scrape config parameters
    thr : float, optional
        SNR threshold used to identify candidate events (default: 7)
    beam : int, optional
        ID of source beam (0 -> 1023 for CHIME, default: 0)
    addr : str, optional
        Identified events can be forwarded via a zmq PUB socket.
        String should be compliant with a ``socket.connect()`` call.
        (default: None)
    """
    def __init__(self, dedisp, thr=7, beam=0, addr=None):
        self.dedisp = dedisp
        self.thr = thr
        self.beam = beam

        self.L1_HDR_DTYPE = L1_event.get_dtypes(dedisp)[0]

        self.start_time = np.datetime64(0, 'us')
        self.cur_events = []
        self.prev_events = []

        self._init_scalings()
        self._init_bowtie_footprints()
        self._init_trigger_buffers()
        self._init_networking(addr)

    @property
    def _ichunk(self):
        return self.dedisp.nchunks_processed - 1

    def _init_scalings(self):
        """
        Relative coarse time & DM scalings for each tree.
        Most likely, these scalings will both be 1, 2, 4, 8...
        """
        nt = self.dedisp.nt_coarse_per_chunk
        ndm = self.dedisp.ndm_coarse
        self.t_scale = nt[0]/np.array(nt)
        self.dm_scale = np.array(ndm)/ndm[0]*2**np.arange(self.dedisp.ntrees)

    def _init_bowtie_footprints(self):
        """
        Defines neighborhoods used for local-maxima condition.

        These are boolean arrays used as the `footprint` for
        ``scipy.ndimage.maximum_filter()`` call in ``_find_peaks()``.

        The smallest bowtie is...
            ``np.array([[1, 1, 1, 0, 0],
                        [0, 1, 1, 1, 0],
                        [0, 0, 1, 1, 1], bool)``
        """
        def bowtie(ndm):
            B = np.zeros((ndm*2+1, ndm*2+3), bool)
            for i in 0, 1, 2:
                np.fill_diagonal(B[:, i:], True)
            return B
        self.bowties = [bowtie(ndm) for ndm in 8, 4, 2, 1, 1, 1, 1]

    def _init_trigger_buffers(self):
        """Create buffers necessary to accommodate local max neighborhoods"""
        shapes = [x.shape for x in self.dedisp.get_triggers()]
        self.bufs = [np.zeros(shape[:-1]+(len(b)+1,))
                     for shape, b in zip(shapes, self.bowties)]

    def _init_networking(self, addr):
        """Setup zmq context and socket"""
        if addr:
            context = zmq.Context.instance()
            self.pub = context.socket(zmq.PUB)
            self.pub.connect(addr)
        else:
            self.pub = None

    def process_triggers(self, triggers=None):
        """Form L1 events from bonsai coarse-grained triggers

        Parameters
        ----------
        triggers: list of array_like, optional
            Coarse-grained triggers generated by bonsai,
            (default: None, i.e. triggers are obtained via
            ``self.dedisp.get_triggers())

        Returns
        -------
        array_like
            L1 Events, with dtype obtained from ``L1_event.get_dtypes()``
        """
        self.chunk = triggers or self.dedisp.get_triggers()
        self._manage_buffers()
        self._find_peaks()
        return self._dump_events()

    def _manage_buffers(self):
        """Prepends buffers to trigger arrays and then updates buffers"""
        for itree, snrs in enumerate(self.chunk):
            self.chunk[itree] = np.concatenate([self.bufs[itree], snrs], -1)
            self.bufs[itree] = snrs[..., -self.bufs[itree].shape[-1]:]

    def _find_peaks(self):
        """
        For each downsampling factor, trial scattering measure, and trial
        spectral index, searches for valid SNR peaks in the (DM, time)
        plane.  A *peak* is defined as a local maxima in a bowtie-shaped
        neighborhood that is above threshold.  Peaks are handled via
        ``_insert_event()``
        """
        def img_peaks(img, footprint):
            above_thr = img > self.thr
            local_max = img == maximum_filter(img, footprint=footprint,
                                              mode='constant', cval=1e3)
            return np.argwhere(np.logical_and(above_thr, local_max))

        for itree, snrs in enumerate(self.chunk):
            sm_beta = product(xrange(snrs.shape[1]), xrange(snrs.shape[2]))
            for sm, beta in sm_beta:
                img = snrs[:, sm, beta, :]
                peaks = img_peaks(img, self.bowties[itree])
                snrs = img[peaks[:, 0], peaks[:, 1]]
                peaks[:, 0] *= self.dm_scale[itree]
                peaks[:, 1] -= self.bufs[itree].shape[-1]
                peaks[:, 1] *= self.t_scale[itree]
                for snr, (dm, t) in zip(snrs, peaks):
                    e = np.array((itree, dm, sm, beta, t, snr), simple_dtype)
                    self._insert_event(e.view(np.recarray))

    def _insert_event(self, e1):
        """
        Adds an event, grouping with other events that are close in DM & time.
        These other events will be attributed to different itree, sm, beta.
        """
        for batch in self.cur_events:
            for e2 in batch:
                E, e = (e1, e2) if e1.itree > e2.itree else (e2, e1)
                E_ddm, E_dt = self.dm_scale[E.itree], self.t_scale[E.itree]
                if max(abs(E.time-e.time)/E_dt, abs(E.dm-e.dm)/E_ddm) <= 2:
                    batch.append(e1)
                    return
        self.cur_events.append([e1])

    def _dump_events(self):
        """Prepare, sift, and send collected events"""
        if len(self.cur_events) is 0 or not self.cur_events:
            if self.pub:
                self.pub.send_multipart([str(self.beam),
                                         str(self._ichunk),
                                         ''])
            return
        dump = np.hstack([self._prepare_event(group)
                          for group in self.cur_events if len(group) != 0])
        dump = self._sift_repeats(dump.view(np.recarray))
        if self.pub:
            self.pub.send_multipart([str(self.beam),
                                     str(self._ichunk),
                                     dump.tostring()])
        self.cur_events = []
        return dump

    def _sift_repeats(self, events, t_thr=0.256, dm_thr=7):
        """
        Removes events coincident with events from the previous chunk.
        This is a fix for a rare corner case where detections of a pulse
        in different trees straddling the chunking seam. This is caused by...
          1. buffers not covering the same length of time
          2. events detected in more downsampled trees peaking slightly later
        """
        if len(self.prev_events) is 0:
            self.prev_events = events
            return events
        prev = self.prev_events
        tmax_prev = prev.time.max()
        tmin_cur = events.time.min()
        icands = np.where((events.time-tmax_prev).astype(float)/1e3 < t_thr)[0]
        prev = prev[(tmin_cur - prev.time).astype(float)/1e3 < t_thr]
        repeats = []
        for i in icands:
            for p in prev:
                if abs(events[i]['dm']-p['dm']) > dm_thr:
                    continue
                if (events[i]['time'] - p['time']).astype(float)/1e3 > t_thr:
                    continue
                repeats.append(i)
                break
        events = events[list(set(xrange(len(events))).difference(repeats))]
        self.prev_events = events
        return events

    def _prepare_event(self, group):
        """
        Transforms a group of *simple headers* into a single L1 event.
        The header with the best snr characterizes the pulse while the
        other headers are used to form the `snr_vs_itree` curve.
        """
        group = np.array(group).view(np.recarray)
        best = group[group.snr.argmax()]
        event = np.zeros(1, self.L1_HDR_DTYPE).view(np.recarray)
        event.beam = self.beam
        event.itree = best['itree']
        event.snr = best['snr']

        # Arrival Time
        buf_len = self.bufs[best['itree']].shape[-1]
        it = best['time']/self.t_scale[best['itree']]+buf_len
        t0 = self._ichunk * self.dedisp.nt_chunk * self.dedisp.dt_sample
        dt = (self.dedisp.nt_chunk * self.dedisp.dt_sample /
              self.dedisp.nt_coarse_per_chunk[best['itree']])
        t = t0 + (it-buf_len)*dt - self.dedisp.trigger_lag_dt[best['itree']]
        event.time = self.start_time + np.timedelta64(int(1e6*t), 'us')

        # DM
        idm = best['dm']/self.dm_scale[best['itree']]
        event.dm = (idm * self.dedisp.max_dm[best['itree']] /
                    self.dedisp.ndm_coarse[best['itree']])

        # SNR curves
        dm_offs = np.arange(1-buf_len/2, buf_len/2).repeat(3)
        dm_offs = dm_offs.reshape(buf_len-1, 3)
        t_offs = np.tile(np.arange(-1, 2), buf_len-1).reshape(buf_len-1, 3)
        t_offs += dm_offs
        bowtie = self.chunk[best['itree']][idm+dm_offs, best['sm'],
                                        best['beta'], it+t_offs]
        event.snr_vs_dm[0, :len(bowtie)] = bowtie.max(1)

        if self.dedisp.nsm[best['itree']] > 1:
            event.snr_vs_sm = self.chunk[best['itree']][idm, :, event.beta[0], it]

        if self.dedisp.nbeta[best['itree']] > 1:
            event.snr_vs_beta = self.chunk[best['itree']][idm, event.sm[0], :, it]

        if self.dedisp.ntrees > 1:
            for itree in xrange(self.dedisp.ntrees):
                try:
                    snr = group[group['itree'] == itree].snr.max()
                    event.snr_vs_itree[0, itree] = snr
                except:
                    pass
        return event[0]
