import rf_pipelines
import numpy as np
import h5py
from math import floor


class variance_estimator(rf_pipelines.py_wi_transform):
    """
    This is a pseudo-transform (meaning that it does not actually modify its input). 
    
    This estimates the variance of the intensity/weights arrays passed to it. First, 
    it computes the variance across v1_chunk number of timesamples. Then, it takes the 
    median of v2_chunk number of v1 estimates. Note that in the case of the majority of 
    a channel being masked, the variance will be 0.

    Thus, the final variance array produced is of dimensions
        (nfreq, total number of time samples / v1_chunk / v1_chunk)

    The variance array is written out as a .h5 file in the directory of the test script.
    The prefix for this outputted file can be determined by the fname argument. 

    Variance arrays are outputted after the v2 array accumulated at least 64 x pixels so that v2 
    does not become too large and too much data is not lost if something causes the
    pipeline run to terminate (note that files are only outputted at the end of a 
    process_chunk call). 
    """

    def __init__(self, v1_chunk=128, v2_chunk=80, nt_chunk=1024, fname=None):
        name = "variance_estimator(v1_chunk=%d, v2_chunk=%d, nt_chunk=%d)" % (v1_chunk, v2_chunk, nt_chunk)

        # Call base class constructor
        rf_pipelines.py_wi_transform.__init__(self, name)

        assert nt_chunk % v1_chunk == 0, \
            'For now, nt_chunk(=%d) must be a multiple of v1_chunk(=%d)' % (nt_chunk, v1_chunk)

        self.v1_chunk = v1_chunk
        self.v2_chunk = v2_chunk
        self.nt_chunk = nt_chunk
        self.nt_prepad = 0
        self.nt_postpad = 0 
        self.v1_t = floor(self.v2_chunk / 2) + 1  # which v1 to index for time (+1 since iv1 is len, not index)
            
        # Make and open the h5 file
        if fname is None:
            self.fname = 'var_v1_%d_v2_%d.h5' % (self.v1_chunk, self.v2_chunk)
        else:
            self.fname = '%s_v1_%d_v2_%d.h5' % (fname, self.v1_chunk, self.v2_chunk)
        self.f = h5py.File(self.fname, mode='w')
        self.f.attrs['v1_chunk'] = self.v1_chunk
        self.f.attrs['v2_chunk'] = self.v2_chunk


    def set_stream(self, s):
        # As explained in the rf_pipelines.py_wi_transform docstring, this function is
        # called once per pipeline run, when the stream (the 's' argument) has been specified.
        self.nfreq = s.nfreq


    def start_substream(self, isubstream, t0):
        # Called once per substream (a stream can be split into multiple substreams).
        self.v1 = np.zeros((self.nfreq, self.v2_chunk))
        self.iv1 = 0   # keeps track of which position we are adding v1 to 
        self.v2 = []   # stores final v2 values (reshaped when being outputted)
        self.t = []    # keeps track of the average time for each v2

        # Create h5 dsets here so we can change them later - can be extended if maxshape is None
        self.v_dset = self.f.create_dataset('variance', shape=(self.nfreq, 0), dtype=np.float32, maxshape=(self.nfreq, None))

        # I coudn't get this to work without being 2-D, so time needs to be accessed by indexing the 0th 
        # element to get a normal list in Variance_Estimates class
        self.t_dset = self.f.create_dataset('time', shape=(1,0), dtype=np.float32, maxshape=(1, None)) 


    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        # This is the main computational routine defining the transform, which is called
        # once per incoming "block" of data.  For documentation on the interface see
        # rf_pipelines.py_wi_transform docstring.
        
        for i in xrange(0, self.nt_chunk, self.v1_chunk):
            # Process the chunks for each frequency sequentially, wrap back if nt_chunk > v1_chunk
            for frequency in xrange(self.nfreq):
                # Compute v1 for each frequency
                a = self._v1(intensity[frequency, i : i+self.v1_chunk], weights[frequency, i : i+self.v1_chunk])
                self.v1[frequency, self.iv1] = a
            self.iv1 += 1

            # Check whether we need to record the time
            if self.iv1 == self.v1_t:
                # Write the current time to self.t. i tells you the number of time samples that have been gone though (ignoring frequency), 
                # so we can use that to interpolate a time. 
                if self.v2_chunk % 2 == 0:
                    self.t += [ 1. * i / (t1 - t0) + t0 ]  # take the first time index in the chunk
                else:
                    self.t += [ 1. * i / (t1 - t0) + t0 - (t1 - t0) / (self.nt_chunk / self.v1_chunk) / 2 ]  # take the middle time index of the previous chunk  

            # Once we have calculated a value for each frequency, check if we should output a v2
            if self.iv1 == self.v2_chunk:
                for frequency in xrange(self.nfreq):
                    # Empty v1[frequency] and compute median for v2 (0 if too many v1 values are 0)
                    if np.count_nonzero(self.v1[frequency]) < self.v2_chunk * 0.25: 
                        self.v2.append(0)
                    else:
                        # Here, we want to ignore elements with value 0
                        self.v2.append(np.median(self.v1[frequency][np.nonzero(self.v1[frequency])]))
                self.v1[frequency, :] = 0
                self.iv1 = 0                

        # Write a file every now and then (for now, after at least 64 x-pixels are produced)
        if len(self.v2) >= self.nfreq*64:
            self._write()


    def end_substream(self):
        """Write any remaining data"""
        if len(self.v2) > 0:
            # We need to remove the extra time that was calculated
            if 1. * len(self.v2) / self.nfreq != len(self.t):
                self.t.pop()
            self._write()
        self.f.close()


    def _v1(self, i, w):
        """Calculate weighted variance"""
        if np.count_nonzero(w) < self.v1_chunk * 0.25:
            return 0
        return np.average(i**2, weights=w) 
    

    def _write(self):
        """Writes a .h5 file"""
        # Reshape self.v2
        out = np.array(self.v2).reshape((self.nfreq, -1), order='F')

        # Write data to script directory
        shape = self.v_dset.shape
        self.v_dset.resize((self.nfreq, shape[1] + len(out[0])))
        self.t_dset.resize((1, shape[1] + len(out[0])))
        self.v_dset[:, shape[1]:] = out[:,:]
        self.t_dset[0, shape[1]:] = self.t[:]
        print 'Variance Estimator: writing variance data to', self.fname

        # Clear self.v2 and self.t
        self.v2 = []
        self.t = []
