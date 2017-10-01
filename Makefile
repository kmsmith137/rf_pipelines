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

INCFILES = rf_pipelines.hpp \
	rf_pipelines_base_classes.hpp \
	rf_pipelines_inventory.hpp \
	rf_pipelines_internals.hpp

# Source files for the core C++ library 'librf_pipelines.so'

OFILES = badchannel_mask.o \
	bonsai_dedisperser.o \
	chime_file_stream.o \
	chime_file_stream_base.o \
	chime_file_writer.o \
        chime_frb_file_stream.o \
	chime_network_stream.o \
	chime_packetizer.o \
	chunked_pipeline_object.o \
	file_utils.o \
	gaussian_noise_stream.o \
	intensity_clippers.o \
	json_utils.o \
	lexical_cast.o \
	outdir_manager.o \
	pipeline.o \
	pipeline_object.o \
	polynomial_detrenders.o \
	ring_buffer.o \
	spline_detrenders.o \
	std_dev_clippers.o \
	wi_sub_pipeline.o \
	wi_stream.o \
	wi_transform.o

# Files that get installed in $(PYDIR)
# Includes both Python source files and the extension module rf_pipelines_c.so (written in C++)
PYFILES=rf_pipelines/rf_pipelines_c.so \
	rf_pipelines/__init__.py \
	rf_pipelines/utils.py \
	rf_pipelines/L1b.py \
	rf_pipelines/L1_event.py \
	rf_pipelines/streams/__init__.py \
	rf_pipelines/streams/chime_streams.py \
	rf_pipelines/transforms/__init__.py \
	rf_pipelines/transforms/adversarial_masker.py \
	rf_pipelines/transforms/bonsai_dedisperser.py \
	rf_pipelines/transforms/frb_injector_transform.py \
	rf_pipelines/transforms/mask_filler.py \
	rf_pipelines/transforms/noise_filler.py \
	rf_pipelines/transforms/plotter_transform.py \
	rf_pipelines/transforms/variance_estimator.py \
	rf_pipelines/retirement_home/__init__.py \
	rf_pipelines/retirement_home/intensity_clipper.py \
	rf_pipelines/retirement_home/polynomial_detrender.py \
	rf_pipelines/retirement_home/std_dev_clipper.py


TESTBINFILES = test-ring-buffer test-core-pipeline-logic test-file-stream-base

# Used in 'make clean'
CLEANDIRS=. site rf_pipelines rf_pipelines/streams rf_pipelines/transforms rf_pipelines/retirement_home

# Used in 'make uninstall': header files which no longer exist, but did exist in previous versions of rf_pipelines
DUMMY_INCFILES=chime_packetizer.hpp chime_file_stream_base.hpp reverter.hpp


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


####################################################################################################


LIBS=

ifeq ($(HAVE_BONSAI),y)
	CPP += -DHAVE_BONSAI
	LIBS += -lbonsai -lhdf5
endif

#ifeq ($(HAVE_PSRFITS),y)
#	CPP += -DHAVE_PSRFITS
#	LIBS += -lpsrfits_utils -lcfitsio
#endif

ifeq ($(HAVE_CH_FRB_IO),y)
	CPP += -DHAVE_CH_FRB_IO
	LIBS += -lch_frb_io -lhdf5
endif

#ifeq ($(HAVE_CH_FRB_IO),y)
#	CPP += -DHAVE_CH_FRB_IO
	LIBS += -lsimpulse
#endif

LIBS += -lrf_kernels -ljsoncpp


####################################################################################################


all: librf_pipelines.so rf_pipelines/rf_pipelines_c.so $(TESTBINFILES)

install: librf_pipelines.so rf_pipelines/rf_pipelines_c.so
	mkdir -p $(INCDIR)/ $(LIBDIR)/ $(PYDIR)/rf_pipelines/streams $(PYDIR)/rf_pipelines/transforms $(PYDIR)/rf_pipelines/retirement_home
	cp -f $(INCFILES) $(INCDIR)/
	for f in $(PYFILES); do cp $$f $(PYDIR)/$$f; done
	cp -f librf_pipelines.so $(LIBDIR)/

uninstall:
	for f in $(INCFILES) $(DUMMY_INCFILES); do rm -f $(INCDIR)/$$f; done
	rm -f $(LIBDIR)/librf_pipelines.so
	rm -rf  $(PYDIR)/rf_pipelines_c.so $(PYDIR)/rf_pipelines/

clean:
	rm -f $(TESTBINFILES) rf_pipeline_0.json
	for d in $(CLEANDIRS); do rm -f $$d/*~ $$d/*.o $$d/*.so $$d/*.pyc; done

%.o: %.cpp $(INCFILES)
	$(CPP) -c -o $@ $<

librf_pipelines.so: $(OFILES)
	$(CPP) $(CPP_LFLAGS) -shared -o $@ $^ $(LIBS)

rf_pipelines/rf_pipelines_c.o: rf_pipelines/rf_pipelines_c.cpp $(INCFILES)
	$(CPP) -Wno-strict-aliasing -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION -c -o $@ $<

rf_pipelines/rf_pipelines_c.so: rf_pipelines/rf_pipelines_c.o librf_pipelines.so
	$(CPP) $(CPP_LFLAGS) -shared -o $@ $< -lrf_pipelines -lpyclops $(LIBS) $(LIBS_PYMODULE)


test-core-pipeline-logic: test-core-pipeline-logic.o $(OFILES)
	$(CPP) $(CPP_LFLAGS) -o $@ $^ $(LIBS)

test-ring-buffer: test-ring-buffer.o $(OFILES)
	$(CPP) $(CPP_LFLAGS) -o $@ $^ $(LIBS)

test-file-stream-base: test-file-stream-base.o $(OFILES)
	$(CPP) $(CPP_LFLAGS) -o $@ $^ $(LIBS)
