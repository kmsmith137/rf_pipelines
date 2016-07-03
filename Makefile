# Makefile.local must define the following variables
#   LIBDIR      install dir for C++ libraries
#   INCDIR      install dir for C++ headers
#   PYDIR       install dir for python modules
#   CPP         C++ compiler command line
#
# Some optional variables which I only use for osx/clang:
#   CPP_LFLAGS      extra linker flags when creating a .so or executable file from .o files
#   LIBS_PYMODULE   any extra libraries needed to link a python extension module (osx needs -lPython)
#
# See site/Makefile.local.* for examples.

include Makefile.local

INCFILES=rf_pipelines.hpp rf_pipelines_internals.hpp

# Source files for the core C++ library 'librf_pipelines.so'
OFILES=bonsai_dedisperser.o \
	chime_file_stream.o \
	gaussian_noise_stream.o \
	misc.o \
	psrfits_stream.o \
	simple_detrender.o \
	wi_run.o \
	wraparound_buf.o

# Files that get installed in $(PYDIR)
# Includes both Python source files and the extension module rf_pipelines_c.so (written in C++)
PYFILES=rf_pipelines/rf_pipelines_c.so \
	rf_pipelines/__init__.py \
	rf_pipelines/utils.py \
	rf_pipelines/transforms/__init__.py \
	rf_pipelines/transforms/frb_injector_transform.py \
	rf_pipelines/transforms/plotter_transform.py

# FIXME generate this from PYFILES using Makefile rule
PYCFILES=rf_pipelines/__init__.pyc \
	rf_pipelines/transforms/__init__.pyc \
	rf_pipelines/utils.pyc \
	rf_pipelines/transforms/frb_injector_transform.pyc \
	rf_pipelines/transforms/plotter_transform.pyc

CLEANDIRS=. rf_pipelines rf_pipelines/transforms

LIBS=


####################################################################################################


ifndef CPP
$(error Fatal: Makefile.local must define CPP variable)
endif

ifndef INCDIR
$(error Fatal: Makefile.local must define INCDIR variable)
endif

ifndef LIBDIR
$(error Fatal: Makefile.local must define LIBDIR variable)
endif

ifndef PYDIR
$(error Fatal: Makefile.local must define PYDIR variable)
endif

ifeq ($(HAVE_BONSAI),y)
	CPP += -DHAVE_BONSAI
	LIBS += -lbonsai -lhdf5
endif

ifeq ($(HAVE_PSRFITS),y)
	CPP += -DHAVE_PSRFITS
	LIBS += -lpsrfits_utils -lcfitsio
endif

ifeq ($(HAVE_CH_FRB_IO),y)
	CPP += -DHAVE_CH_FRB_IO
	LIBS += -lch_frb_io -lhdf5
endif


####################################################################################################


all: librf_pipelines.so rf_pipelines/rf_pipelines_c.so run-unit-tests

install: librf_pipelines.so rf_pipelines/rf_pipelines_c.so
	mkdir -p $(INCDIR)/ $(LIBDIR)/ $(PYDIR)/rf_pipelines/transforms
	cp -f $(INCFILES) $(INCDIR)/
	for f in $(PYFILES); do cp $$f $(PYDIR)/$$f; done
	cp -f librf_pipelines.so $(LIBDIR)/

uninstall:
	for f in $(INCFILES); do rm -f $(INCDIR)/$$f; done
	for f in $(PYFILES) $(PYCFILES); do rm -f $(PYDIR)/$$f; done
	rmdir $(PYDIR)/rf_pipelines/transforms
	rmdir $(PYDIR)/rf_pipelines
	rm -f $(LIBDIR)/librf_pipelines.so

clean:
	rm -f run-unit-tests
	for d in $(CLEANDIRS); do rm -f $$d/*~ $$d/*.o $$d/*.so $$d/*.pyc; done

%.o: %.cpp $(INCFILES)
	$(CPP) -c -o $@ $<

librf_pipelines.so: $(OFILES)
	$(CPP) $(CPP_LFLAGS) -shared -o $@ $^ $(LIBS)

rf_pipelines/rf_pipelines_c.so: rf_pipelines/rf_pipelines_c.cpp $(INCFILES) rf_pipelines/python_extension_helpers.hpp librf_pipelines.so
	$(CPP) $(CPP_LFLAGS) -Wno-strict-aliasing -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION -shared -o $@ $< -lrf_pipelines $(LIBS) $(LIBS_PYMODULE)

run-unit-tests: run-unit-tests.o librf_pipelines.so
	$(CPP) $(CPP_LFLAGS) -o $@ $< -lrf_pipelines
