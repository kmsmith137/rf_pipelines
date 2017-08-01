import numpy as np
import numpy.random
import random

import rf_pipelines


class adversarial_masker(rf_pipelines.py_wi_transform):
    """
    A half-finished transform, intended to stress-test the online variance estimation logic
    by masking a bunch of rectangle-shaped regions of the input array.  If the online variance
    estimation logic is working well, then the trigger plots shouldn't show any false positives!

    The constructor argument 'nt_reset' is not very well thought-out, but is vaguely intended
    to be the timescale (expressed as a number of samples) over which the variance estimation
    logic resets itself.

    TODO: right now, this isn't very "adversarial", since all it does is mask a few large
    rectangles which correspond to timestream gaps of different sizes.  Let's try to test
    as many weird cases as possible!
    """

    def __init__(self, nt_chunk=1024, nt_reset=65536):
        rf_pipelines.py_wi_transform.__init__(self, 'adversarial_masker(nt_chunk=%d, nt_reset=%d)' % (nt_chunk, nt_reset))
        self.nt_prepad = 0
        self.nt_postpad = 0
        self.nt_chunk = nt_chunk
        self.nt_reset = nt_reset


    def set_stream(self, s):
        self.nfreq = s.nfreq

        # Initialize:
        #   self.rectangles: list of (freq_lo, freq_hi, t_lo, t_hi) quadruples to be masked.
        #   self.nt_max: nominal timestream size, expressed as number of samples.
        #   self.nt_processed: samples processed so far (this counter is incremented in process_chunk())
  
        self.rectangles = [ ]
        self.nt_max = self.nt_reset 
        self.nt_processed = 0
        
        # We place timestream gaps of a few different sizes, separated by 'nt_reset' samples.
        for i in xrange(8):
            nt_rect = self.nt_reset // 2**i    # gap size
            self.rectangles += [ (0, self.nfreq, self.nt_max, self.nt_max + nt_rect) ]
            self.nt_max = self.nt_max + nt_rect + self.nt_reset


        # Let's try masking random frequency channels within each of our rectangles
        for i in xrange(8):
            nt_rect = self.nt_reset // 2**i

            # Number of channels masked (force at least 25% to be masked):
            n_masked = random.randint(self.nfreq // 4, self.nfreq-1)
            print "Adversarial masker, rectangle " + str(i) + ": Masking " + str(n_masked) + " of " + str(self.nfreq) + " channels."

            n = 0
            while n < n_masked:
                # Get a random number between [0, nfreq-1]
                masked = random.randint(0, self.nfreq-1)
                self.rectangles += [ (masked, masked+1, self.nt_max, self.nt_max + nt_rect) ]
                n += 1

            self.nt_max = self.nt_max + nt_rect + self.nt_reset


        # Try randomly masking individual samples -- oops this didn't work out; much too slow!! Maybe make a new "mask" member 
        # variable than can be more effeciently applied in process_chunk?
        # for i in xrange(8):
        #     nt_rect = self.nt_reset // 2**i
        #     # Generate random times to mask between self.nt_max and self.nt_mask + nt_rect -- mask 50% of samples
        #     masked_t = np.random.random_integers(self.nt_max, self.nt_max + nt_rect, size=(nt_rect * self.nfreq // 2))
        #     masked_f = np.random.random_integers(0, self.nfreq, size=(nt_rect * self.nfreq // 2))
        #     for j in xrange(nt_rect * self.nfreq // 2):
        #         self.rectangles += [ (masked_f[j], masked_f[j] + 1, masked_t[j], masked_t[j] + 1) ]
        #     self.nt_max = self.nt_max + nt_rect + self.nt_reset


        #  Try masking groups of adjacent frequencies
        for i in xrange(8):
            nt_rect = self.nt_reset // 2**i

            # Decide what number of frequencies to mask
            n_masked = random.randint(self.nfreq // 4, self.nfreq-1)

            # Pick a start point for masking
            start = random.randint(0, self.nfreq)
            
            # See if we need to wrap around to mask
            if start + n_masked > self.nfreq - start:
                self.rectangles += [ (0, n_masked - (self.nfreq - start), self.nt_max, self.nt_max + nt_rect) ]

            # Mask from start -> end or nfreq
            self.rectangles += [ (start, min(self.nfreq, start + n_masked), self.nt_max, self.nt_max + nt_rect) ]
            
            self.nt_max = self.nt_max + nt_rect + self.nt_reset


        # Vertical stripes!
        for i in xrange(8):
            nt_rect = self.nt_reset // 2**i
            
            # Make each stripe 1/8 of a rectangle (this is totally arbitrary and should be changed to something
            # more sensible!)
            stripe_width = nt_rect // 8
            
            # Mask! 
            stripe_end = self.nt_max
            for i in xrange(4):
                self.rectangles += [ (0, self.nfreq, stripe_end, stripe_end + stripe_width) ]
                stripe_end += 2 * stripe_width
                
            self.nt_max = self.nt_max + nt_rect + self.nt_reset
            

        # Horizontal stripes! These are kind of lame at present and can be made more interesting. 
        for i in xrange(8):
            nt_rect = self.nt_reset // 2**i
            
            # Make each stripe 1/8 of a rectangle (this is totally arbitrary and should be changed to something
            # more sensible!)
            stripe_width = self.nfreq // 8
            
            # Mask! 
            stripe_end = 0
            for i in xrange(4):
                self.rectangles += [ (stripe_end, stripe_end + stripe_width, self.nt_max, self.nt_max + nt_rect) ]
                stripe_end += 2 * stripe_width
                
            self.nt_max = self.nt_max + nt_rect + self.nt_reset
            

        # Debug
        # for rect in self.rectangles:
        #     print rect

    
    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        # Loop over rectangles
        for (ifreq0, ifreq1, it0, it1) in self.rectangles:
            # Range of indices in current chunk which overlap with the rectangle
            it0 = max(0, min(self.nt_chunk, it0 - self.nt_processed))
            it1 = max(0, min(self.nt_chunk, it1 - self.nt_processed))

            # If there is an overlap, then mask (by setting weights to zero)
            if it0 < it1:
                weights[ifreq0:ifreq1,it0:it1] = 0.0

        self.nt_processed += self.nt_chunk
