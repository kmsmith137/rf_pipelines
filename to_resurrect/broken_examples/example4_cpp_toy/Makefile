# You'll probably need to do 'make install' from the toplevel rf_pipelines directory
# before this will compile.

include ../../Makefile.local

all: example4

example4: example4.cpp
	$(CPP) $(CPP_LFLAGS) -o $@ $< -lrf_pipelines

clean:
	rm -f *~ example4
