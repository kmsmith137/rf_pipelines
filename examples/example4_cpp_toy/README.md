### DESCRIPTION

This is a minimal rf_pipelines run which demonstrates running an rf_pipeline through
its C++ interface.  It makes a stream of Gaussian random numbers, applies a (superfluous)
detrender, and runs the output through the bonsai dedispersion code.

### INSTRUCTIONS FOR RUNNING

You'll need to compile the bonsai code, and make sure that rf_pipelines has been compiled
and installed with HAVE_BONSAI=y in Makefile.local.

Then you should be able to run the example as follows.
For more details, see comments in example4.cpp.
```
# Compile example4.cpp, either by hand or using the micro-Makefile in the example4 directory
make

# You'll need to generate the bonsai config hdf5 file from the bonsai text file.
# (This is a temporary workaround for a currently-unimplemented feature in bonsai: on-the-fly
# estimation of trigger variances.)
bonsai-mkweight bonsai_config.txt bonsai_config.hdf5

# Run the pipeline.  We now run bonsai in a mode where it splits the output across multiple hdf5 files.
./example4

# Plot the output of the dedispersion transform.
# This actually generates six plots, since the triggers have been split into two hdf5 files, and the bonsai
# config file defines three dedispersion trees correpsonding to different DM and pulse width ranges.
bonsai-plot-triggers.py bonsai_outputs.hdf5
```
