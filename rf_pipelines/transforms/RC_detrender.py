import numpy as np
import rf_pipelines 

class RC_detrender(rf_pipelines.py_wi_transform):
    """
    In filter jargon, this detrender is composed of two single pole
    recursive infinite impulse response high-pass filters (yikes).
    One scans forwards in time and the other backwards (and is therefore
    not a true example of inifinite impulse response..)

    Naming reason: The filter responses are equivalent to analog RC circuit

    Bi-directionality: Used to mitigate the switched noise source.
    Step functions which occur on long time scales can be dealt with by 
    comparing the response of each direction and choosing the most stable 
    (the step introduces a recovery time in the filter's response, so the
    direction that hasn't experienced the jump is used instead). 

    Constructor syntax:
        
        t = RC_detrender(a=0.99, nt_chunk=1024)

        'a=0.99' determines the time constant of the filter's exponential decay.

        'nt_chunk=1024' is the buffer size
    """

    def __init__(self,nt_chunk=1024,a=0.99):
        self.nt_chunk = nt_chunk
        self.a = a

    def set_stream(self, stream):
        self.nfreq = stream.nfreq
        self.fw  = np.zeros((self.nfreq,self.nt_chunk)) 
        self.bw  = np.zeros((self.nfreq,self.nt_chunk))
        self.dfw = np.zeros(self.fw.shape)
        self.dbw = np.zeros(self.bw.shape)

    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        # forwards
        for i in xrange(self.nt_chunk):
            self.fw[:,i]  = np.where( weights[:,i], 
                                self.fw[:,i-1]*self.a + intensity[:,i]*(1-self.a),
                                self.fw[:,i-1])
            self.dfw[:,i] = np.where( weights[:,i],
                                self.dfw[:,i-1]*self.a + (self.fw[:,i]-self.fw[:,i-1])*(1-self.a),
                                self.dfw[:,i-1])
        # backwards
        self.bw[:,-1] = self.fw[:,-1] 
        for i in xrange(1,self.nt_chunk):
            self.bw[:,-i-1]  = np.where( weights[:,-i-1],
                                   self.bw[:,-i]*self.a + intensity[:,-i-1]*(1-self.a),
                                   self.bw[:,-i])
            self.dbw[:,-i-1] = np.where( weights[:,-i-1],
                                   self.dbw[:,-i]*self.a + (self.bw[:,-i-1]-self.bw[:,-i])*(1-self.a),
                                   self.dbw[:,-i])
        # choose the most stable filter response to detrend with... 
        intensity -= np.where(abs(self.dfw) < abs(self.dbw), self.fw, self.bw)
