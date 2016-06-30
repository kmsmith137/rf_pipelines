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

OFILES=bonsai_dedisperser.o \
	chime_file_stream.o \
	detrenders.o \
	gaussian_noise_stream.o \
	misc.o \
	psrfits_stream.o \
	wi_run.o \
	wraparound_buf.o

PYFILES=rf_pipelines_c.so rf_pipelines.py

LIBS=

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

all: librf_pipelines.so rf_pipelines_c.so run-unit-tests

install: librf_pipelines.so rf_pipelines_c.so
	cp -f $(INCFILES) $(INCDIR)/
	cp -f $(PYFILES) $(PYDIR)/
	cp -f librf_pipelines.so $(LIBDIR)/

uninstall:
	for f in $(INCFILES); do rm -f $(INCDIR)/$$f; done
	for f in $(PYFILES); do rm -f $(PYDIR)/$$f; done
	rm -f $(LIBDIR)/librf_pipelines.so

clean:
	rm -f *~ *.o *.so run-unit-tests

%.o: %.cpp $(INCFILES)
	$(CPP) -c -o $@ $<

librf_pipelines.so: $(OFILES)
	$(CPP) $(CPP_LFLAGS) -shared -o $@ $^ $(LIBS)

rf_pipelines_c.so: rf_pipelines_c.cpp $(INCFILES) python_extension_helpers.hpp librf_pipelines.so
	$(CPP) $(CPP_LFLAGS) -Wno-strict-aliasing -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION -shared -o $@ $< -lrf_pipelines $(LIBS) $(LIBS_PYMODULE)

run-unit-tests: run-unit-tests.o librf_pipelines.so
	$(CPP) $(CPP_LFLAGS) -o $@ $< -lrf_pipelines
