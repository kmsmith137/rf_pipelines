### RF_PIPELINES

rf_pipelines: A plugin-based framework for processing channelized intensity data in radio astronomy.

Note: this repo now includes the web viewer code which was previously
in [mburhanpurkar/web_viewer](https://github.com/mburhanpurkar/web_viewer).

Warning: At the moment, rf_pipelines is **undocumented**.
This is a big problem, and fixing it is a high priority!
(There used to be some minimal documentation and examples, but I temporarily removed them,
since they were so out-of-date that they were more misleading than helpful.)

In the meantime, if you're interested in running the real-time CHIME FRB pipeline,
please refer to the [kmsmith137/ch_frb_l1](https://github.com/kmsmith137/ch_frb_l1/)
repository, which does have reasonable documentation.  If you're interested in
running the "offline" CHIME FRB pipeline (including the web viewer), you may
find the example scripts in the [mrafieir/ch_frb_rfi](https://github.com/mrafieir/ch_frb_rfi/)
repository useful.

Installation instructions follow.

### REQUIRED DEPENDENCIES

  - First note that if you build rf_pipelines with only its required dependencies
    (i.e. no optional dependencies) then it will not be very useful!  For
    example, you can't dedisperse data without the 'bonsai' optional dependency,
    you can't read/write common file formats without optional dependencies
    such as 'hdf5', 'psrfits', etc.

    For CHIME, I recommend building rf_pipelines with the following optional
    dependencies: ch_frb_io, bonsai, simpulse, hdf5, png.  See the next section
    [Optional Dependencies](#user-content-optional-dependencies) for more info.
    
  - A gcc which is recent enough that C++11 is supported.  
    I know that gcc 4.8.1 works, and that 4.4.7 is too old.

  - A recent 2.x python (I know that 2.7.5 and 2.7.12 work).

  - The following helper libraries (see their respective README's for compilation instructions):
      - [kmsmith137/simd_helpers](https://github.com/kmsmith137/simd_helpers):
        header-only library for writing x86 assembly language kernels.
      - [kmsmith137/pyclops](https://github.com/kmsmith137/pyclops):
        some hacks for writing hybrid C++/python code.
      - [kmsmith137/rf_kernels](https://github.com/kmsmith137/rf_kernels):
        fast C++/assembly kernels for RFI removal and related tasks.

  - The following python modules: numpy, scipy, h5py, Pillow.  
    (Installation hint: `sudo pip install <module_name>`)

  - jsoncpp (https://github.com/open-source-parsers/jsoncpp)
      - osx one-liner: `brew install jsoncpp`
      - centos one-liner: `yum install jsoncpp-devel`
      - ubuntu one-liner: `apt-get install libjsoncpp-dev`
      - Building from scratch is a pain, but the following procedure worked for me:
        ```
        git clone https://github.com/open-source-parsers/jsoncpp
        mkdir -p build/debug
        cd build/debug
        cmake -DCMAKE_INSTALL_PREFIX:PATH=$HOME -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_C_FLAGS=-fPIC -DCMAKE_BUILD_TYPE=debug -G "Unix Makefiles" ../..
        make install
        ```

<a name="optional-dependencies"></a> 
### OPTIONAL DEPENDENCIES

For each optional dependency below, there is a corresponding Makefile variable
(e.g. `HAVE_BONSAI`).  In our build system, these variables are declared in Makefile.local
(see the next section for more info).

  - **bonsai:** very optimized x86 incoherent dedispersion using a tree algorithm.

    See installation instructions in its github repo:
    [kmsmith137/bonsai](https://github.com/kmsmith137/bonsai).
    Note that bonsai has its own optional dependencies, which
    you may want to include when building it.

  - **ch_frb_io:** CHIMEFRB file formats and network I/O.

    See installation instructions in its github repo:
    [CHIMEFRB/ch_frb_io](https:://github.com/CHIMEFRB/ch_frb_io).

    Note: ch_frb_io is really two libraries, a python library and a C++ library, 
    and you only need the C++ part.  The C++ part has its own dependencies
    (lz4, msgpack, zeromq, cppzmq, hdf5, bitshuffle).

  - **simpulse:** C++/python library for simulating pulses in radio astronomy

    See installation instructions in its github repo:
    [kmsmith137/simpulse](https://github.com/kmsmith137/simpulse)

    Note that simpulse requires cython and fftw3.

  - **libpng:** Strictly speaking, this is an optional dependency, but you'll
    almost certainly want to enable it, since making png images for visualization
    purposes is a really useful feature in rf_pipelines.
      - osx one-liner: `brew install libpng`
      - centos one-liner: `yum install libpng-devel`
      - ubuntu one-liner: `apt-get install libpng-dev`

  - **hdf5:** Widely-used file format for scientific data.
  
    You'll need to install HDF5 (including C++ support), if it's not already installed.
    You'll also need the sp_hdf5 header-only C++ wrapper library (https://github.com/kmsmith137/sp_hdf5).
    See the README.md file there for installation instructions.

    **Note that sp_hdf5 requires HDF5 version 1.8.12 or later,
    but does not work with version 1.10.x.  This will be fixed eventually!**
    In the meantime, the sp_hdf5 README includes instructions for installing a version of HDF5
    which is neither too old nor too new.
      
  - **psrfits_utils:** Utility library for working with search- and fold-mode PSRFITS pulsar data files.
    Needed for e.g. gbncc.  (Not needed for CHIME.)

    See installation instructions in its github repo:
    [scottransom/psrfits_utils](https://github.com/scottransom/psrfits_utils)
    
    Note: if psrfits_utils installation fails due to autoconf warnings which don't look 
          serious, it may help to edit configure.ac and remove the "-Werror" in the line 
          AM_INIT_AUTOMAKE([-Wall -Werror foreign])


### INSTALLATION

  - The rf_pipelines Makefile assumes the existence of a file `Makefile.local` which defines
    a few machine-dependent Makefile variables:
    ```
      INCDIR     Installation directory for C++ header files
      LIBDIR     Installation directory for libraries
      CPP        C++ compiler executable + flags, see below for tips!
        etc.
    ```

    For a complete list of variables which must be defined, see comments at the top of ./Makefile.

    Rather than write a Makefile.local from scratch, I recommend that you start with one of the
    examples in the site/ directory, which contains Makefile.locals for a few frequently-used
    CHIME machines.  In particular, site/Makefile.local.kms_laptop16 is a recent osx machine,
    and site/Makefile.local.frb1 is a recent CentOS Linux machine.  (If you're a member of
    CHIME and you're using one of these machines, you can just symlink the appropriate file in
    site/ to ./Makefile.local)

  - Do `make all install` to build.

  - If you have trouble getting rf_pipelines to build/work, then the problem probably has
    something to do with your compiler flags (specified as part of CPP) or environment 
    variables.  Here are a few hints:

      - You probably need `-std=c++11` in your compiler flags, for C++11 support
      - You probably need `-pthread` in your compiler flags, in order to compile
        some multithreaded timing tests.
      - You probably need `-march=native` in your compiler flags, to get AVX/AVX2
        intrinsics.  (I usually use optimization flags `-O3 -march=native -ffast-math -funroll-loops`.)
      - You probably want `-Wall -fPIC` in your compiler flags on general principle.
      - The rf_pipelines build procedure assumes that the current directory is searched for header
        files and libraries, i.e. you should have `-I. -L.` in your compiler flags.
      - You also probably want `-I$(INCDIR) -L$(LIBDIR)` in your compiler flags, so that
        these install dirs are also searched for headers/libraries (e.g. simpulse)
      - You may need more -I and -L flags to find all necessary headers/libraries.
      - In particular, if you get the error message "Python.h not found", then you
        probably need something like -I/usr/include/python2.7.  You can get the header
	directory for your version of python with `distutils.sysconfig.get_python_inc()`
      - If you get the error message "numpy/arrayobject.h not found", then you probably 
        need something like -I/usr/lib64/python2.7/site-packages/numpy/core/include.
        You can get the header directory for your numpy installation with 
	`numpy.get_include()`.
      - If everything compiles but libraries are not being found at runtime, then you
        probably need to add `.` or LIBDIR to the appropriate environment variable
        ($LD_LIBRARY_PATH in Linux, or $DYLD_LIBRARY_PATH in osx)
      - Likewise, if you're not able to load `rf_pipelines` python module, then you probably
        need to add your PYDIR to the $PYTHONPATH environment variable.  If you're
	not able to run the rf_pipelines command-line utils (e.g. `rfp-time`)
	then you probably need to add your BINDIR to the $PATH environment variable.

    Feel free to email me if you have trouble!
