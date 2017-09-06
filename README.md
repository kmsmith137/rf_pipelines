rf_pipelines: A plugin-based framework for processing channelized intensity data in radio astronomy.

Some day, there will be documentation in [MANUAL.md](./MANUAL.md).

For a high-level overview, here are some slides from CHIME telecons:

  - June 20 proposal:
      [ docs/16-06-20-rf_pipelines_proposal.pdf ] 
      ( docs/16-06-20-rf_pipelines_proposal.pdf ) 

  - July 5 update:
      [ docs/16-07-05-rf_pipelines_update.pdf ] 
      ( docs/16-07-05-rf_pipelines_update.pdf ) 

This code is under construction!  A current to do list is in [TODO.md](./TODO.md).

### DEPENDENCIES

  - A gcc which is recent enough that C++11 is supported.  
    I know that gcc 4.8.1 works, and that 4.4.7 is too old.

  - The simd_helpers library (https://github.com/kmsmith137/simd_helpers)
    This is new and experimental, so let me know if you have any trouble compiling it!

  - libhdf5 (https://www.hdfgroup.org/HDF5/release/obtain5.html)
    Note that this is a link to HDF5 v1.8.  I imagine v1.10 also works but haven't tested it yet.

  - The following python modules: numpy, scipy, h5py.  (Installation hint: `sudo pip install <module_name>`)

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

  - Optional but recommended: The 'PIL' python imaging library (you can test whether you have 
    it with 'import PIL' from python).  If you need to install it, I recommend the 'Pillow' 
    variant (pip install Pillow)

  - Optional: simpulse (https://github.com/kmsmith137/simpulse)
    You'll need this if you want to inject simulated FRB's.
    Note that simpulse requires cython and fftw3 (see its README).

  - Optional: bonsai (https://github.com/CHIMEFRB/bonsai)
    This is an incremental incoherent dedisperser.
    Note that bonsai v3 or later is required!
    Note that bonsai requires a very recent cython (see its README).

  - Optional: ch_frb_io (https://github.com/CHIMEFRB/ch_frb_io)
    You'll need this if you want to analyze CHIME data.  

    Note: ch_frb_io is really two libraries, a python library and a C++ library, 
    and you only need the C++ part (see installation instructions in the ch_frb_io README).

  - Optional: psrfits_utils (https://github.com/scottransom/psrfits_utils)
    You'll need this if you want to analyze data in PSRFITS format (e.g. gbncc).

    Note: if psrfits_utils installation fails due to autoconf warnings which don't look 
          serious, it may help to edit configure.ac and remove the "-Werror" in the line 
          AM_INIT_AUTOMAKE([-Wall -Werror foreign])


### INSTALLATION

  - Create a file ./Makefile.local containing compiler flags, library locations, etc.
    The format is defined in Makefile, but it will probably be easiest to copy one of
    the examples in site/ and customize.  For each of the optional dependencies above,
    there is a y/n flag in Makefile.local to indicate whether you have it.

    Note that your C++ compiler flags should probably include at minimum
       -std=c++11 -O3 -march=native -fPIC -I. -I$(INCDIR) -L. -L$(LIBDIR)
    and additionally -pthread if using g++.  If -march=native is omitted, you'll
    get a lot of compiler errors!

  - make all install

  - For some quick unit tests, do 
```
      ./run-unit-tests
      ./test-kernels
      ./test-file-stream-base
      ./test-cpp-python-equivalence.py
```

### QUICK START

  - This code is best "documented by example", so defintely start by checking out the following:

    - [examples/example1_toy] (examples/example1_toy):
      Illustrates basic interface + writing a new transform in python

    - [examples/example2_gbncc] (examples/example2_gbncc):
      Detrend and dedisperse GBNCC data

    - [examples/example3_chime] (examples/example3_chime):
      Detrend and dedisperse CHIME pathfinder data

    - [examples/example4_cpp_toy] (examples/example4_cpp_toy):
      Illustrates running rf_pipelines through its C++ interface

  - If you want to run rf_pipelines from python, I recommend browsing python docstrings next.  Some useful docstrings:
    ```
    rf_pipelines                    [ 'help(rf_pipelines)' from python interpreter. ]
    rf_pipelines.py_wi_transform    [ 'help(rf_pipelines.py_wi_transform)' ]
 
    + docstrings for stream/transform objects, e.g. rf_pipelines.chime_stream_from_acqdir
      or rf_pipelines.plotter_transform
    ```

  - If you want to write a new stream in C++, I recommend reading comments in rf_pipelines.hpp, especially the
    comments in 'class wi_stream' and 'class wi_run_state'.

  - For a specification of the pipeline json output, see [./docs/json_spec.txt](docs/json_spec.txt).
