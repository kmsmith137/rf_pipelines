# Makefile.local for cedar.computecanada.ca

HAVE_PSRFITS=n
HAVE_SIMPULSE=y
HAVE_CH_FRB_IO=y
HAVE_BONSAI=y
HAVE_HDF5=y
HAVE_PNG=y

# Directory where C++ libraries will be installed
LIBDIR=$(HOME)/lib

# Directory where C++ header files will be installed
INCDIR=$(HOME)/include

# Directory where executables will be installed
BINDIR=$(HOME)/bin

# Directory where Python modules will be installed
PYDIR=$(HOME)/lib/python2.7/site-packages

# Read-only location of python header files
PYTHON_INCDIR=/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/python/2.7.14/include/python2.7

# Hint for NUMPY_INCDIR: import numpy; numpy.get_include()
# On cedar, I'm assuming numpy has been installed with 'pip install --user'.
NUMPY_INCDIR=$(HOME)/.local/lib/python2.7/site-packages/numpy/core/include

#
# C++ command line
# Must support c++11
# Don't forget to put -L. and -L$(LIBDIR) on the command line (in this order)
# Don't forget to add . and $(LIBDIR) in your LD_LIBRARY_PATH environment variable (in this order)
# Don't forget -pthread and -fPIC
# Don't forget -march=native
#

CPP=g++ -std=c++11 -pthread -fPIC -Wall -Wno-literal-suffix -O3 -march=native -ffast-math -funroll-loops -I. -I$(INCDIR) -I/usr/include/jsoncpp -I$(PYTHON_INCDIR) -I$(NUMPY_INCDIR) -L. -L$(LIBDIR)

CPP_LFLAGS= -L. -L$(LIBDIR)
