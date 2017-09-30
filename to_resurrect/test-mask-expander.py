#!/usr/bin/env python

import numpy as np
from rf_pipelines.transforms.mask_expander import expand_mask


def test_expand_mask(nfreq, nt_chunk, thr, axis):
    num = 0.    # count of unmasked array elements
    den = 0.    # count of all array elements

    for i in xrange(1000):
        weights = np.random.normal(thr, 0.01, size=(nfreq,nt_chunk))
        expand_mask(weights, thr, axis)
        num += np.count_nonzero(weights)
        den += weights.size

    print 'unmasked_fraction =', (num/den)


for axis in [ None, 0, 1 ]:
    test_expand_mask(256, 512, 0.2, axis)
