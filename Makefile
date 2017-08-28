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

INCFILES=rf_pipelines.hpp rf_pipelines_internals.hpp chime_packetizer.hpp chime_file_stream_base.hpp reverter.hpp

KERNEL_INCFILES=kernels/downsample.hpp \
	kernels/intensity_clippers.hpp \
	kernels/mask.hpp \
	kernels/mean_variance.hpp \
	kernels/polyfit.hpp \
	kernels/std_dev_clippers.hpp

# Source files for the core C++ library 'librf_pipelines.so'
OFILES=badchannel_mask.o \
	bitmask_maker.o \
	bonsai_dedisperser.o \
	chime_file_stream.o \
	chime_frb_file_stream.o \
	chime_file_stream_base.o \
	chime_file_writer.o \
	chime_network_stream.o \
	chime_packetizer.o \
	gaussian_noise_stream.o \
	intensity_clippers.o \
	lexical_cast.o \
	misc.o \
	online_mask_filler.o \
	outdir_manager.o \
	polynomial_detrenders.o \
	reverter.o \
	pulse_adder.o \
	psrfits_stream.o \
	spline_detrenders.o \
	std_dev_clippers.o \
	timing_thread.o \
	udsample.o \
	wi_run_state.o \
	wi_stream.o \
	wi_transform.o \
	wraparound_buf.o

# Files that get installed in $(PYDIR)
# Includes both Python source files and the extension module rf_pipelines_c.so (written in C++)
PYFILES=rf_pipelines/rf_pipelines_c.so \
	rf_pipelines/__init__.py \
	rf_pipelines/utils.py \
	rf_pipelines/L1b.py \
	rf_pipelines/L1_event.py \
	rf_pipelines/streams/__init__.py \
	rf_pipelines/streams/chime_streams.py \
	rf_pipelines/streams/psrfits_stream.py \
	rf_pipelines/streams/gaussian_noise_stream.py \
	rf_pipelines/transforms/__init__.py \
	rf_pipelines/transforms/adversarial_masker.py \
	rf_pipelines/transforms/chime_packetizer.py \
	rf_pipelines/transforms/chime_transforms.py \
	rf_pipelines/transforms/bonsai_dedisperser.py \
	rf_pipelines/transforms/frb_injector_transform.py \
	rf_pipelines/transforms/plotter_transform.py \
	rf_pipelines/transforms/badchannel_mask.py \
	rf_pipelines/transforms/polynomial_detrender.py \
	rf_pipelines/transforms/mask_expander.py \
	rf_pipelines/transforms/kurtosis_filter.py \
	rf_pipelines/transforms/std_dev_clipper.py \
	rf_pipelines/transforms/thermal_noise_weight.py \
	rf_pipelines/transforms/RC_detrender.py \
	rf_pipelines/transforms/intensity_clipper.py \
	rf_pipelines/transforms/variance_estimator.py \
	rf_pipelines/transforms/mask_filler.py \
	rf_pipelines/transforms/noise_filler.py \
	rf_pipelines/transforms/online_mask_filler.py

TESTBINFILES=run-unit-tests \
	test-kernels \
	test-file-stream-base \
	time-clippers

# Used in 'make clean'
CLEANDIRS=. site rf_pipelines rf_pipelines/streams rf_pipelines/transforms \
	examples/example1_toy examples/example2_gbncc examples/example3_chime


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

ifeq ($(HAVE_PSRFITS),y)
	CPP += -DHAVE_PSRFITS
	LIBS += -lpsrfits_utils -lcfitsio
endif

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
	mkdir -p $(INCDIR)/ $(LIBDIR)/ $(PYDIR)/rf_pipelines/streams $(PYDIR)/rf_pipelines/transforms
	cp -f $(INCFILES) $(INCDIR)/
	for f in $(PYFILES); do cp $$f $(PYDIR)/$$f; done
	cp -f librf_pipelines.so $(LIBDIR)/

uninstall:
	for f in $(INCFILES); do rm -f $(INCDIR)/$$f; done
	rm -f $(LIBDIR)/librf_pipelines.so
	rm -rf $(PYDIR)/rf_pipelines/

clean:
	rm -f $(TESTBINFILES)
	for d in $(CLEANDIRS); do rm -f $$d/*~ $$d/*.o $$d/*.so $$d/*.pyc; done

%.o: %.cpp $(INCFILES) $(KERNEL_INCFILES)
	$(CPP) -c -o $@ $<

librf_pipelines.so: $(OFILES)
	$(CPP) $(CPP_LFLAGS) -shared -o $@ $^ $(LIBS)

rf_pipelines/rf_pipelines_c.so: rf_pipelines/rf_pipelines_c.cpp $(INCFILES) rf_pipelines/python_extension_helpers.hpp librf_pipelines.so
	$(CPP) $(CPP_LFLAGS) -Wno-strict-aliasing -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION -shared -o $@ $< -lrf_pipelines $(LIBS) $(LIBS_PYMODULE)

run-unit-tests: run-unit-tests.o librf_pipelines.so
	$(CPP) $(CPP_LFLAGS) -o $@ $< -lrf_pipelines $(LIBS)

test-kernels: test-kernels.cpp $(INCFILES) $(KERNEL_INCFILES)
	$(CPP) -o $@ $<

test-file-stream-base: test-file-stream-base.cpp $(INCFILES) $(KERNEL_INCFILES) librf_pipelines.so
	$(CPP) $(CPP_LFLAGS) -o $@ $< -lrf_pipelines $(LIBS)

time-clippers: time-clippers.cpp $(INCFILES) $(KERNEL_INCFILES) librf_pipelines.so
	$(CPP) $(CPP_LFLAGS) -o $@ $< -lrf_pipelines $(LIBS)
