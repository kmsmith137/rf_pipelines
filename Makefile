# Makefile.local must define the following variables
#   LIBDIR      install dir for C++ libraries
#   INCDIR      install dir for C++ headers
#   PYDIR       install dir for python/cython modules
#   CPP         C++ compiler command line
#   CPP_LFLAGS  extra linker flags when creating a .so or executable file from .o files
#
# See site/Makefile.local.* for examples.

include Makefile.local

INCFILES=rf_pipelines.hpp rf_pipelines_internals.hpp
PYFILES=rf_pipelines.py rf_pipelines_cython.so
LIBFILES=librf_pipelines.so

OFILES=chime_file_stream.o \
	detrenders.o \
	misc.o \
	psrfits_stream.o \
	wi_run.o \
	wraparound_buf.o

LIBS=

ifeq ($(HAVE_PSRFITS),y)
	CPP += -DHAVE_PSRFITS
	LIBS += -lpsrfits_utils -lcfitsio
endif

ifeq ($(HAVE_CH_FRB_IO),y)
	CPP += -DHAVE_CH_FRB_IO
	LIBS += -lch_frb_io -lhdf5
endif


####################################################################################################


all: librf_pipelines.so rf_pipelines_cython.so run-unit-tests

install: librf_pipelines.so rf_pipelines_cython.so
	cp -f $(INCFILES) $(INCDIR)/
	cp -f $(LIBFILES) $(LIBDIR)/
	cp -f $(PYFILES) $(PYDIR)/

uninstall:
	for f in $(INCFILES); do rm -f $(INCDIR)/$$f; done
	for f in $(LIBFILES); do rm -f $(LIBDIR)/$$f; done
	for f in $(PYFILES); do rm -f $(PYDIR)/$$f; done

clean:
	rm -f *~ *.o *.so rf_pipelines_cython.cpp run-unit-tests

%.o: %.cpp $(INCFILES)
	$(CPP) -c -o $@ $<

librf_pipelines.so: $(OFILES)
	$(CPP) $(CPP_LFLAGS) -shared -o $@ $^ $(LIBS)

rf_pipelines_cython.cpp: rf_pipelines_cython.pyx rf_pipelines_pxd.pxd
	cython --cplus $<

rf_pipelines_cython.so: rf_pipelines_cython.cpp rf_pipelines_cython.hpp librf_pipelines.so
	$(CPP) $(CPP_LFLAGS) -shared -o $@ $< -lrf_pipelines $(LIBS)

run-unit-tests: run-unit-tests.o librf_pipelines.so
	$(CPP) $(CPP_LFLAGS) -o $@ $< -lrf_pipelines
