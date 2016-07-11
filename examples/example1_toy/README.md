### DESCRIPTION

This is a toy example illustrating some of the builtin transforms, plus
writing a new transform in Python.

  - We start with a stream which emits Gaussian random noise

  - We apply a 'toy_transform' written in python which injects a periodic
    sinusoidal signal plus a sinusoidal mask.  (The details here are arbitrary
    and just an excuse to demonstrate the python interface.)

  - We add a simulated FRB.

  - We run a plotter which plots the input timestream consisting of the sum
    of the signals above

  - We apply the bonsai dedisperser and generate coarse-grained triggers.


### INSTRUCTIONS FOR RUNNING

First you'll need to generate the bonsai config hdf5 file from the bonsai text file.
(This is a temporary workaround for a currently-unimplemented feature in bonsai: on-the-fly
estimation of trigger variances.)

```
    bonsai-mkweight bonsai_config.txt bonsai_config.hdf5
```

Then run the example:

```
   ./example1.py
```

This will generate a bunch of waterfall plots plus a file 'triggers.hdf5' containing
coarse-grained triggers.  The trigger file can be plotted with:

```
   bonsai-plot-triggers.py triggers.hdf5
```

If everything worked then the file 'waterfall_5.png' should contain the simulated FRB, 
gaussian noise, plus periodic artifacts from the toy_transform.  The file 'triggers.png'
should contain some wide diagonal stripes due to the toy_transform artifacts, plus a
sharply peaked bowtie from the FRB.  The github repo contains reference versions of these
plots for comparison:

![alt](https://github.com/kmsmith137/rf_pipelines/blob/master/examples/example1_toy/reference_waterfall_5.png "Text")

!(https://github.com/kmsmith137/rf_pipelines/blob/master/examples/example1_toy/reference_triggers.png)
