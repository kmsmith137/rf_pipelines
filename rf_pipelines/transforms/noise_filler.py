import numpy as np
import time
import rf_pipelines


class noise_filler(rf_pipelines.py_wi_transform):
    """
    This adds simulated Gaussian noise to intensity arrays and leaves the weights arrays unchanged. 

    The current implementation is as follows:
    - the variance of the noise added is constant across all time samples in a given frequency channel
      for each call to process_chunk
    - the variance will increase or decrease accrording to the value of self.increment between each
      process_chunk call
    - the initial variance will be random for each frequency (you need to change the values of the
      random starting values generated at the beginning as well)

    Extra variables information:
    - nt_chunk: if testing with variance_estimator, choose a value close to v1_chunk * v2_chunk - default
      is for 16 and 32 for quick testing. 
    - increment: how much the variance is incremented by with each call to process_chunk - default is 
      0.1 for a reasonable value given initial variance is between 0.1 and 0.9
    - current_var: a nfreq length list containing the variance values currently being indexed by 
      process_chunk
    - var_accumulator: holds all the variance values used in a pipeline run - written out as a .npy 
      file in end_substream

    """

    def __init__(self, nt_chunk=512, increment=0.1):
        name = "noise_filler(nt_chunk=%d, increment=%f)" % (nt_chunk, increment)
        self.nt_chunk = nt_chunk
        self.increment = increment
        rf_pipelines.py_wi_transform.__init__(self, name)


    def set_stream(self, s):
         self.nfreq = s.nfreq


    def start_substream(self, isubstream, t0):
        # If increment is -ve, this may work nicely
        # self.current_var = np.random.randint(500, 550, (self.nfreq))

        # If increment is +ve, this works
        self.current_var = (np.random.ranf(size=self.nfreq) + 0.01) * 2

        self.var_accumulator = []    # Reshape at the end 


    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        # First, cycle through the chunk and replace values with simulated ones using self.current_var
        for f in range(self.nfreq):
            variance = self.current_var[f]
            intensity[f] = np.random.normal(scale=variance, size=(self.nt_chunk))                

        # Add self.current_var to accumulator
        self.var_accumulator += list(self.current_var)

        # Increment self.current_var
        self.current_var += self.increment

        # Write out if things are getting too large
        if len(self.var_accumulator) >= self.nfreq * 64:
            self._write()


    def end_substream(self):
        self._write()


    def _write(self):
       out = np.array(self.var_accumulator).reshape((self.nfreq, -1), order='F')**2
       np.save('simulated_var_%s' % (time.strftime('%y-%m-%d-%X')), out)
       print 'Noise Filler: wrote', 'simulated_var_%s' % (time.strftime('%y-%m-%d-%X')) + '.npy'
       self.var_accumulator = []
