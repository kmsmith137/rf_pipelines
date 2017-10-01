import numpy as np
import h5py
from math import floor

from rf_pipelines.rf_pipelines_c import pipeline_object, wi_transform


class variance_estimator(wi_transform):
    """
    This pseudo-transform (meaning that it does not actually modify its input)
    estimates the variance of an intensity array and writes the results into a file. 
    First, it computes the variance across 'v1_chunk' number of time samples. Then, it 
    computes the median of 'v2_chunk' number of 'v1_chunk' estimates. Thus, the final 
    variance array produced is of dimensions (nfreq, total number of time samples / 
    v1_chunk / v1_chunk). The 'var_filename' argument is used to store the variance 
    estimates in a .h5 file.
   
    Notes: 
       - If 75% of a freq channel is masked, then the variance is set to 0.
       - The variance arrays are outputted after accumulating the first 64 bins.
    """

                 
    def __init__(self, var_filename, v1_chunk=128, v2_chunk=80, nt_chunk=1024):
        name = "variance_estimator(var_filename=%s, v1_chunk=%d, v2_chunk=%d, nt_chunk=%d)" % (var_filename, v1_chunk, v2_chunk, nt_chunk)

        assert var_filename is not None
        assert isinstance(var_filename, basestring)   # if this fails, arguments are probably in the old ordering

        # Call base class constructor
        wi_transform.__init__(self, name)

        assert nt_chunk % v1_chunk == 0, \
            'For now, nt_chunk(=%d) must be a multiple of v1_chunk(=%d)' % (nt_chunk, v1_chunk)

        self.v1_chunk = v1_chunk
        self.v2_chunk = v2_chunk
        self.nt_chunk = nt_chunk
        self.v1_t = floor(self.v2_chunk / 2) + 1  # which v1 to index for time (+1 since iv1 is len, not index)
        self.var_filename = var_filename

        self.f = h5py.File(self.var_filename, mode='w')
        self.f.attrs['v1_chunk'] = self.v1_chunk
        self.f.attrs['v2_chunk'] = self.v2_chunk


    def _bind_transform(self, json_attrs):
        if not json_attrs.has_key('t_initial'):
            raise RuntimeError("rf_pipelines.variance_estimator: pipeline must contain a chime_file_stream (or another stream which defines the 't_initial' attribute)")
        
        self.t_initial = json_attrs['t_initial']
        self.dt_sample = json_attrs['dt_sample']
        
        self.v1 = np.zeros((self.nfreq, self.v2_chunk))
        self.iv1 = 0   # keeps track of which position we are adding v1 to 
        self.v2 = []   # stores final v2 values (reshaped when being outputted)
        self.t = []    # keeps track of the average time for each v2

        # Create h5 dsets here so we can change them later - can be extended if maxshape is None
        self.v_dset = self.f.create_dataset('variance', shape=(self.nfreq, 0), dtype=np.float32, maxshape=(self.nfreq, None))

        # I coudn't get this to work without being 2-D, so time needs to be accessed by indexing the 0th 
        # element to get a normal list in Variance_Estimates class
        self.t_dset = self.f.create_dataset('time', shape=(1,0), dtype=np.float32, maxshape=(1, None)) 


    def _process_chunk(self, intensity, weights, pos):
        t0 = self.t_initial + self.dt_sample * pos
        t1 = self.t_initial + self.dt_sample * (pos + self.nt_chunk)
        
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

        # Write a file every now and then (for now, after at least 64 bins are produced)
        if len(self.v2) >= self.nfreq*64:
            self._write()
                 
    def _end_pipeline(self, json_output):
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
        print 'Variance Estimator: writing variance data to', self.var_filename

        # Clear self.v2 and self.t
        self.v2 = []
        self.t = []


    def jsonize(self):
        return { 'class_name': 'variance_estimator',
                 'var_filename': self.var_filename,
                 'v1_chunk': self.v1_chunk,
                 'v2_chunk': self.v2_chunk,
                 'nt_chunk': self.nt_chunk }


    @staticmethod
    def from_json(j):
        return variance_estimator(var_filename = j['var_filename'],
                                  v1_chunk = j['v1_chunk'],
                                  v2_chunk = j['v2_chunk'],
                                  nt_chunk = j['nt_chunk'])


pipeline_object.register_json_constructor('variance_estimator', variance_estimator.from_json)
