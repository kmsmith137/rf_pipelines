#!/usr/bin/env python

import rf_pipelines
import numpy as np
import h5py
import sys


# A command line tool for plotting variances h5 files and making comparison plots 
# (dividing one h5 file by another). If called with one argument, it will plot 
# the h5 file and if called with three, it will produce a comparison plot. 
# This is not very friendly and is written assuming the user uses it correctly :)


def read_h5(fname):
    """Helper function to extract variances + interpolate from h5 files"""
    with h5py.File(fname, 'r') as hf:
        var = hf['variance'][:]
    # Interpolate zeros
    x = len(var[0])
    indices = np.arange(x)
    for frequency in xrange(len(var)):
        nonzero = np.nonzero(var[frequency])[0]
        if len(nonzero) < 0.25 * x:
            var[frequency] = np.zeros((x))
        else:
            var[frequency] = np.interp(indices, nonzero, var[frequency, nonzero])
    return var


def plot(h5_file):
    """Plot a single h5 file"""
    name = h5_file[:-3] + '.png'
    var = read_h5(h5_file)
    print name, 'shape:', var.shape

    rf_pipelines.utils.var_to_png(name, var, comparison=False)


def compare(h5_file1, h5_file2, factor):
    """Make a comparison plot for two h5 files"""
    # Read files and remove zeros in the denominator
    var1 = read_h5(h5_file1)
    var2 = read_h5(h5_file2)
    var2[var2 < 0.00001] = 0.00001

    # Start with smallest one, upsample to nt * t, and make the same size
    if var1.shape[1] > var2.shape[1]:
        var2 = rf_pipelines.utils.upsample(var2, var2.shape[0], var2.shape[1] * int(factor)) 
        var1 = var1[:,:var2.shape[1]]
    else:
        var1 = rf_pipelines.utils.upsample(var1, var1.shape[0], var1.shape[1] * int(factor))
        var2 = var2[:,:var1.shape[1]]

    comparison = var1/var2
    name = 'comp_' + h5_file1[:-3] + '_' + h5_file2[:-3] + '.png'

    rf_pipelines.utils.var_to_png(name, comparison, comparison=True)


if __name__ == "__main__":
    if len(sys.argv) == 2:
        # Plot a single variance file
        plot(sys.argv[1])
    else:
        # Three args - var1, var2, conversion factor
        compare(sys.argv[1], sys.argv[2], sys.argv[3])
