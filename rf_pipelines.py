# This file is currently just a wrapper around rf_pipelines_cython, but
# intended to be easier to read than a .pyx file.

import rf_pipelines_cython


def psrfits_stream(filename):
    return wi_stream(rf_pipelines_cython.make_psrfits_stream(filename))


def simple_detrender(nt_chunk):
    return wi_transform(rf_pipelines_cython.make_simple_detrender(nt_chunk))


####################################################################################################


class wi_stream(object):
    """
    Currently, there is no way to write a wi_stream from Python.
    Therefore, the only option is to use a wrapped C++ stream such as psrfits_stream() above.
    """

    def __init__(self, p):
        assert isinstance(p, rf_pipelines_cython.wi_stream)
        self.p = p

    @property
    def nfreq(self):
        return self.p.get_nfreq()

    @property
    def freq_lo_MHz(self):
        return self.p.get_freq_lo_MHz()

    @property
    def freq_hi_MHz(self):
        return self.p.get_freq_hi_MHz()

    @property
    def dt_sample(self):
        return self.p.get_dt_sample()

    @property
    def nt_maxwrite(self):
        return self.p.get_nt_maxwrite()

    def run(self, transforms):
        """The 'transforms' arg should be a list of objects of class wi_transform."""

        self.p.clear_transforms()

        for t in transforms:
            assert isinstance(t, wi_transform)
            self.p.add_transform(t.p)
            
        self.p.run()
        self.p.clear_transforms()


class wi_transform(object):
    def __init__(self, p):
        assert isinstance(p, rf_pipelines_cython.wi_transform)
        self.p = p
