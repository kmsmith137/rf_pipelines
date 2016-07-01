import sys
import numpy as np

try:
    import PIL.Image
except ImportError:
    print >>sys.stderr, "rf_pipelines: import PIL.Image failed; many things will work but plotting will fail"


from rf_pipelines_c import \
   wi_stream, \
   wi_transform, \
   make_psrfits_stream, \
   make_gaussian_noise_stream, \
   make_simple_detrender, \
   make_bonsai_dedisperser


####################################################################################################


def write_png(filename, arr, weights=None, transpose=False, ytop_to_bottom=False):
    arr = np.array(arr, dtype=np.float)
    assert arr.ndim == 2

    if weights is not None:
        weights = np.array(weights, dtype=np.float)
        assert weights.shape == arr.shape
    
    if not transpose:
        arr = np.transpose(arr)
        weights = np.transpose(weights) if (weights is not None) else None
    
    if not ytop_to_bottom:
        arr = arr[::-1]
        weights = weights[::-1] if (weights is not None) else None

    (mean, rms) = (np.mean(arr), np.var(arr)**0.5)

    # color in range [0,1].
    color = 0.5 + 0.16*(arr-mean)/rms    # factor 0.16 preserves convention from some old code
    color = np.maximum(color, 0.0)
    color = np.minimum(color, 0.999999)  # 0.99999 instead of 1.0, to make roundoff-robust
    
    # rgb in range [0,1]
    red = 256. * color
    blue = 256. * (1-color)

    if weights is not None:
        (wmin, wmax) = (np.min(weights), np.max(weights))
        if wmin < 0:
            raise RuntimeError('write_png: negative weights are currently treated as an error')

        # avoid NaN
        if wmax == 0.0:
            wmax = 1.0

        weights = weights/wmax
        red *= weights
        blue *= weights

    rgb = np.zeros((arr.shape[0],arr.shape[1],3), dtype=np.uint8)
    rgb[:,:,0] = red
    rgb[:,:,2] = blue

    img = PIL.Image.fromarray(rgb)
    img.save(filename)
    print >>sys.stderr, 'wrote %s' % filename


def _downsample_2d(arr, new_nfreq, new_ntime):
    """Helper for wi_downsample."""

    assert arr.ndim == 2
    assert new_nfreq > 0
    assert new_ntime > 0

    (nfreq, ntime) = arr.shape
    assert nfreq % new_nfreq == 0
    assert ntime % new_ntime == 0

    arr = np.reshape(arr, (new_nfreq, nfreq//new_nfreq, new_ntime, ntime//new_ntime))
    arr = np.sum(arr, axis=3)
    arr = np.sum(arr, axis=1)
    return arr
    

def wi_downsample(intensity, weights, new_nfreq, new_ntime):
    """Downsamples a pair of 2D arrays (intensity, weights), returning a new pair (ds_intensity, ds_weights)."""

    wi = _downsample_2d(weights * intensity, new_nfreq, new_ntime)
    w = _downsample_2d(weights, new_nfreq, new_ntime)

    mask = (w > 0.0)
    wi = wi / np.where(mask, w, 1.0)
    wi = np.where(mask, wi, 0.0)

    return (wi, w)


class plotter_transform(wi_transform):
    def __init__(self, img_prefix, img_nfreq, img_nt, nt_chunk):
        assert nt_chunk > 0
        assert img_nfreq > 0
        assert img_nt > 0
        assert nt_chunk % img_nt == 0

        self.img_prefix = img_prefix
        self.img_nfreq = img_nfreq
        self.img_nt = img_nt
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0

        self.substream_ix = 0
        self.chunk_ix = 0


    def set_stream(self, s):
        self.nfreq = s.nfreq
        assert s.nfreq % self.img_nfreq == 0


    def start_substream(self, t0):
        pass


    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weight):
        prefix = self.img_prefix if (self.substream_ix == 0) else ('%s%d' % (self.img_prefix,self.substream_ix+1))
        filename = '%s_%d.png' % (prefix, self.chunk_ix)
        self.chunk_ix += 1

        (ds_intensity, ds_weights) = wi_downsample(intensity, weights, self.img_nfreq, self.img_nt)
        write_png(filename, intensity, weights=weights, transpose=True, ytop_to_bottom=True)


    def end_substream(self):
        self.substream_ix += 1



class frb_injector_transform(wi_transform):
    def __init__(self, snr, undispersed_arrival_time, dm, intrinsic_width=0.0, sm=0.0, spectral_index=0.0, sample_rms=1.0, nt_chunk=1024):
        self.dm = dm
        self.sm = sm
        self.snr = snr
        self.spectral_index = spectral_index
        self.intrinsic_width = intrinsic_width
        self.undispersed_arrival_time = undispersed_arrival_time
        self.sample_rms = sample_rms

        self.nfreq = 0     # to be initialized in set_stream()
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0


    def set_stream(self, stream):
        try:
            import simpulse
        except ImportError:
            # Print an installation hint...
            print >>sys.stderr, "\n*** Note: If you don't have the 'simpulse' library, it can be downloaded at https://github.com/kmsmith137/simpulse ***\n"
            # ...then print an error message and exit
            raise

        self.nfreq = stream.nfreq
        self.freq_lo_MHz = stream.freq_lo_MHz
        self.freq_hi_MHz = stream.freq_hi_MHz
        self.dt_sample = stream.dt_sample

        self.pulse = simpulse.single_pulse(1024,       # number of samples used "under the hood" in the simpulse library
                                           self.nfreq,
                                           self.freq_lo_MHz,
                                           self.freq_hi_MHz,
                                           self.dm,
                                           self.sm,
                                           self.intrinsic_width,
                                           1.0,        # fluence, placeholder value to be changed below based on S/N
                                           self.spectral_index,
                                           self.undispersed_arrival_time)

        # signal-to-noise with fluence=1
        snr0 = self.pulse.get_signal_to_noise(self.dt_sample)
        
        # adjust fluence of pulse to target S/N
        self.pulse.fluence = self.snr / snr0

    
    def start_substream(self, t0):
        # We keep track of the time range spanned by the stream, so that we can print a warning
        # in end_substream() if the substream doesn't span the pulse.
        self.substream_t0 = t0
        self.substream_t1 = t0


    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weight):
        nt_chunk = intensity.shape[1]
        t1 = t0 + nt_chunk * self.dt_sample

        # The freq_hi_to_lo flag tells the 'simpulse' library to use the rf_pipelines
        # frequency ordering (i.e. frequencies ordered from high to low).
        self.pulse.add_to_timestream(intensity, t0, t1, freq_hi_to_lo=True)
        self.substream_t1 = t1


    def end_substream(self):
        # We print warning if the pulse isn't entirely contained in the substream,
        # since this is probably unintentional.
        (pulse_t0, pulse_t1) = self.pulse.get_endpoints()

        if (pulse_t0 >= self.substream_t0) and (pulse_t1 <= self.substream_t1):
            return

        intersection_t0 = max(pulse_t0, self.substream_t0)
        intersection_t1 = min(pulse_t1, self.substream_t1)
        intersection_dt = max(intersection_t1 - intersection_t0, 0)
        missing_frac = 1.0 - intersection_dt / (pulse_t1 - pulse_t0)

        print >>sys.stderr, 'frb_injector_transform: warning: %f\% of pulse was outside stream endpoints' % (100. * missing_frac)
