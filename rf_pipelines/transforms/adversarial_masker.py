import numpy as np
import numpy.random
import random

import rf_pipelines


def round_up(m, n):
    """Rounds up 'm' to a multiple of 'n'."""
    return ((m+n-1)//n) * n


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

    def __init__(self, nt_chunk=1024, nt_reset=65536, nt_minefield=4096*1024):
        rf_pipelines.py_wi_transform.__init__(self, 'adversarial_masker(nt_chunk=%d, nt_reset=%d, nt_minefield=%d)' % (nt_chunk, nt_reset, nt_minefield))
        self.nt_prepad = 0
        self.nt_postpad = 0
        self.nt_chunk = nt_chunk
        self.nt_reset = nt_reset
        self.nt_minefield = nt_minefield


    def set_stream(self, s):
        self.nfreq = s.nfreq

        # Initialize:
        #   self.rectangles: list of (freq_lo, freq_hi, t_lo, t_hi) quadruples to be masked.
        #   self.mask: list of (t_lo, t_hi) to be masked
        #   self.nt_max: nominal timestream size, expressed as number of samples.
        #   self.nt_processed: samples processed so far (this counter is incremented in process_chunk())
  
        self.rectangles = [ ]
        self.mask = [ ]
        self.nt_max = self.nt_reset 
        self.nt_processed = 0
        
        print 'adversarial_masker (nt=%d): timestream gaps of different sizes, with all frequencies masked' % self.nt_max

        for i in xrange(8):
            nt_rect = self.nt_reset // 2**i    # gap size
            self.rectangles += [ (0, self.nfreq, self.nt_max, self.nt_max + nt_rect) ]
            self.nt_max += nt_rect + self.nt_reset

        print 'adversarial_masker (nt=%d): long stretches with (3/4) of the frequency band masked' % self.nt_max

        self.nt_max = round_up(self.nt_max, 2048)
        self.rectangles += [ (self.nfreq//4, self.nfreq, self.nt_max, self.nt_max + self.nt_reset) ]
        self.nt_max += 2 * self.nt_reset

        self.nt_max = round_up(self.nt_max, 2048)
        self.rectangles += [ (0, (3*self.nfreq)//4, self.nt_max, self.nt_max + self.nt_reset) ]
        self.nt_max += 2 * self.nt_reset
        
        print 'adversarial_masker (nt=%d): masking random frequency channels within each of our rectangles' % self.nt_max

        for i in xrange(8):
            nt_rect = self.nt_reset // 2**i

            # Number of channels masked (force at least 25% to be masked):
            n_masked = random.randint(self.nfreq // 4, self.nfreq-1)
            print "    Adversarial masker, rand. rectangle " + str(i) + ": Masking " + str(n_masked) + " of " + str(self.nfreq) + " channels."

            n = 0
            while n < n_masked:
                # Get a random number between [0, nfreq-1]
                masked = random.randint(0, self.nfreq-1)
                self.rectangles += [ (masked, masked+1, self.nt_max, self.nt_max + nt_rect) ]
                n += 1

            self.nt_max += nt_rect + self.nt_reset

        print 'adversarial_masker (nt=%d): masking groups of adjacent frequencies' % self.nt_max

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
            
            self.nt_max += nt_rect + self.nt_reset

        print 'adversarial_masker (nt=%d): vertical stripes!' % self.nt_max

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
                
            self.nt_max += nt_rect + self.nt_reset
            
        print 'adversarial_masker (nt=%d): horizontal stripes!' % self.nt_max

        n_stripes = 40
        for i in xrange(8):
            nt_rect = self.nt_reset // 2**i
            
            # Make each stripe 1/40 of a rectangle (this is totally arbitrary and should be changed to something
            # more sensible!)
            stripe_width = self.nfreq // (2 * n_stripes)
            
            # Mask! 
            stripe_end = 0
            for i in xrange(n_stripes):
                self.rectangles += [ (stripe_end, stripe_end + stripe_width, self.nt_max, self.nt_max + nt_rect) ]
                stripe_end += 2 * stripe_width
                
            self.nt_max += nt_rect + self.nt_reset
            
        print 'adversarial_masker (nt=%d): rectangles with randomly masked samples' % self.nt_max

        for i in xrange(8):
            nt_rect = self.nt_reset // 2**i
            self.mask += [ (self.nt_max, self.nt_max + nt_rect) ]
            self.nt_max += nt_rect + self.nt_reset

        print 'adversarial_masker (nt=%d): "minefield" of randomly placed rectangles' % self.nt_max

        for i in xrange(self.nt_minefield // 3000):
            freq0 = random.randint(-int(0.2*self.nfreq), self.nfreq-1)
            freq1 = random.randint(0, int(1.2*self.nfreq))
            freq_lo = min(freq0,freq1)
            freq_hi = max(freq0,freq1) + 1
            freq_lo = max(freq_lo, 0)
            freq_hi = min(freq_hi, self.nfreq)
            assert 0 <= freq_lo < freq_hi <= self.nfreq
            
            # cube root
            t = random.uniform(-1.,1.)
            t = t**(1/3.) if (t >= 0.) else -(-t)**(1/3.)
            it0 = self.nt_max + int((t+1)/2. * self.nt_minefield)

            nt = int(np.exp(random.uniform(np.log(100.), np.log(10000.))))

            self.rectangles += [ (freq_lo, freq_hi, it0, it0+nt) ]

        self.nt_max += self.nt_minefield

        print 'adversarial_masker (nt=%d): all done!  (nrect=%d)' % (self.nt_max, len(self.rectangles))


    
    def process_chunk(self, t0, t1, intensity, weights, pp_intensity, pp_weights):
        # Loop over rectangles
        for (ifreq0, ifreq1, it0, it1) in self.rectangles:
            # Range of indices in current chunk which overlap with the rectangle
            it0 = max(0, min(self.nt_chunk, it0 - self.nt_processed))
            it1 = max(0, min(self.nt_chunk, it1 - self.nt_processed))

            # If there is an overlap, then mask (by setting weights to zero)
            if it0 < it1:
                weights[ifreq0:ifreq1,it0:it1] = 0.0

        # Loop over mask
        for (it0, it1) in self.mask:
            # Range of indices in current chunk which overlap with the rectangle
            it0 = max(0, min(self.nt_chunk, it0 - self.nt_processed))
            it1 = max(0, min(self.nt_chunk, it1 - self.nt_processed))

            # If there is an overlap, then mask randomly
            if it0 < it1:
                p = 0.50
                a = np.random.choice(a=[False, True], size=weights[:, it0:it1].shape, p=[p, 1-p])
                weights[a] = 0

        self.nt_processed += self.nt_chunk
