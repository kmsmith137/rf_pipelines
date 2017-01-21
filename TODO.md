### High-priority missing features

- Make the bonsai trigger array inspectable from python.
  (This is next on my todo list.)

- Normalize the bonsai triggers.  This may prove to be a long-term effort which needs a lot of tweaking!
  (This is #2 on my todo list.)

- Write a transform which makes a running estimate of the variance in each channel, and writes it to a file,
  which can be read and used in other transforms.  This is a prerequisite for correctly normalizing the
  frb_injector transform.  It may also help with normalizing the bonsai triggers.

- Low-level pipeline logic for saving the intensity/weights arrays for later use.

- Low-level pipeline logic for downsampling/upsampling the pipeline.

- Correctly normalize the frb_injector_transform (needed eventually for end-to-end testing.)

- There should be a way for python transforms to set arbitrary json output.
  (Right now, they can add plot_groups and plots to the json, but that's all.)

### Low-priority missing features

- Low-level pipeline logic for skipping ahead in the pipeline

- Low-level pipeline logic for "forking" the pipeline

- Low-level pipeline logic for multithreading the pipeline

- Can we find a systematic way of looking for memory/refcount leaks?

### Cleanups

- The pipeline needs a lot more documentation!!

- Eventually we'll get rid of the 'imitate_cpp' optional argument on some transforms.

- Some of the examples in the examples/ directory are currently broken.

- Allowing 'axis' arugments to be strings { 'freq', 'time' } in python would be a nice cleanup.

### Performance

- There are still some optimizations/improvements which could be made to the assembly
  language kernels.  For a complete list see [TODO_kernels.md](./TODO_kernels.md).
  My plan for now is to wait until the RFI pipeline converges, see which transforms 
  end up dominating the running time, and making a special effort to optimize those.
  
- Stride-tweaking!  Due to cache associativity issues, the transforms can run faster when
  the stride is not a large power of two.  Here is an example (on frb1):

  ```
  [kendrick@psrcontrol rf_pipelines]$ ./time-detrenders 1024 1024 4096 2 6
  nthreads = 6
  time-detrenders: nfreq=1024, nt_chunk=1024, stride=4096, polydeg=2, niter=16
  polynomial_detrender_cpp(nt_chunk=1024, axis=1, polydeg=2, epsilon=0.01): 0.0343157 sec
  polynomial_detrender_cpp(nt_chunk=1024, axis=0, polydeg=2, epsilon=0.01): 0.106008 sec

  [kendrick@psrcontrol rf_pipelines]$ ./time-detrenders 1024 1024 4112 2 6
  nthreads = 6
  time-detrenders: nfreq=1024, nt_chunk=1024, stride=4112, polydeg=2, niter=16
  polynomial_detrender_cpp(nt_chunk=1024, axis=1, polydeg=2, epsilon=0.01): 0.0362637 sec
  polynomial_detrender_cpp(nt_chunk=1024, axis=0, polydeg=2, epsilon=0.01): 0.0428308 sec
  ```

  I suspect that when the multithreaded pipeline is implemented, and we start doing
  test node runs, we'll want to experiment with implementing stride-tweaking in the
  core rf_pipelines logic.  I think a good place to start is to compare the actual
  pipeline timings with the result of 'time-detrenders', 'time-clippers, etc.

  Reminder: playing with vtune might be interesting.
