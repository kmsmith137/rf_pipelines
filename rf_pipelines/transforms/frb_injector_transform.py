import sys
import rf_pipelines


class frb_injector_transform(rf_pipelines.wi_transform):
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

    
    def start_substream(self, isubstream, t0):
        # We keep track of the time range spanned by the stream, so that we can print a warning
        # in end_substream() if the substream doesn't span the pulse.
        self.substream_t0 = t0
        self.substream_t1 = t0


    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
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
