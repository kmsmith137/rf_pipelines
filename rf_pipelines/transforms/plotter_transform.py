import sys
import numpy as np
import rf_pipelines


class plotter_transform(rf_pipelines.wi_transform):
    def __init__(self, img_prefix, img_nfreq, img_nt, downsample_nt=1, nt_chunk=0):
        assert img_nt > 0
        assert img_nfreq > 0
        assert downsample_nt > 0
        
        if nt_chunk == 0:
            nt_chunk = min(img_nt, 1024//downsample_nt+1) * downsample_nt    # default value

        assert nt_chunk > 0

        if nt_chunk % downsample_nt != 0:
            raise RuntimeError("plotter_transform: current implementation requires 'nt_chunk' to be a multiple of 'downsample_nt'")

        self.img_prefix = img_prefix
        self.img_nfreq = img_nfreq
        self.img_nt = img_nt
        self.nt_chunk_ds = nt_chunk // downsample_nt

        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0


    def set_stream(self, s):
        self.nfreq = s.nfreq

        if s.nfreq % self.img_nfreq != 0:
            raise RuntimeError("plotter_transform: current implementation requires 'img_nfreq' to be a divisor of stream nfreq")


    def start_substream(self, isubstream, t0):
        self.intensity_buf = np.zeros((self.img_nfreq,self.img_nt), dtype=np.float32)
        self.weight_buf = np.zeros((self.img_nfreq,self.img_nt), dtype=np.float32)
        self.isubstream = isubstream
        self.ifile = 0
        self.ipos = 0


    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        (intensity, weights) = rf_pipelines.wi_downsample(intensity, weights, self.img_nfreq, self.nt_chunk_ds)

        ichunk = 0
        while ichunk < self.nt_chunk_ds:
            n = min(self.nt_chunk_ds - ichunk, self.img_nt - self.ipos)
            
            self.intensity_buf[:,self.ipos:(self.ipos+n)] = intensity[:,ichunk:(ichunk+n)]
            self.weight_buf[:,self.ipos:(self.ipos+n)] = weights[:,ichunk:(ichunk+n)]
            self.ipos += n
            ichunk += n

            if self.ipos == self.img_nt:
                self._write_file()


    def end_substream(self):
        if self.ipos > 0:
            self._write_file()


    def _write_file(self):
        filename = self.img_prefix
        if self.isubstream > 0:
            filename += str(isubstream+1)
        filename += ('_%s.png' % self.ifile)

        rf_pipelines.write_png(filename, self.intensity_buf[:,:self.ipos], weights=self.weight_buf[:,:self.ipos], transpose=True, ytop_to_bottom=True)

        self.ifile += 1
        self.ipos = 0
