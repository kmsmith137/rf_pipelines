### RF_PIPELINES

rf_pipelines: A plugin-based framework for processing channelized intensity data in radio astronomy.

Note: this repo now includes the web viewer code which was previously
in [mburhanpurkar/web_viewer](https://github.com/mburhanpurkar/web_viewer).

This manual is mostly a placeholder, and just contains a few random sections,
rather than systematic documentation.  Fixing this is a high priority!


### CONTENTS

  - [Summary of changes since v15](#user-content-summary-of-changes-since-v15)
  - [Command-line utilities](#user-content-command-line-utilities)
  - [Current status](#user-content-current-status)
  - [Class hierarchy](#user-content-class-hierarchy)


<a name="summary-of-changes-since-v15"></a>
### SUMMARY OF CHANGES SINCE V15

In v16, the low-level logic in rf_pipelines was largely rewritten!  The
low-level building blocks (streams, transforms, etc.) have a different API now,
and any "pre-v16" code will need substantial changes to work "post-v16".
This is part of the "2017 Mega Merge" affecting many parts of the CHIMEFRB pipeline.

To make matters worse, due to the CHIME pre-deployment scramble, the new API
is fragile and incomplete, and not all pre-v16 features are working.  We hope
that everything will be finished and well-tested in a few more revisions.

For now, here is a bullet-point summary of changes made since v15.
This will be expanded into more complete documentation later!

  - One important thing to say up-front is that a lot of functionality has not been python-wrapped!
    As a result, the ability to define new pipeline_objects in python is pretty limited, compared to
    the C++ part of the code.  Fixing this is a high priority!

  - Previously, the pipeline operated on a fixed pair of ring buffers, 'intensity' and 'weights',
    and there were two types of objects: "wi_streams", which output the (intensity, weights) buffers,
    and "wi_transforms" which operate on (intensity, weights) buffers in an input/output sense.  A
    pipeline is always a wi_stream followed by N wi_transforms.

  - Now, a pipeline can contain any number of ring buffers (labeled by a 'name' string),
    and pipeline_objects can have M inputs and N outputs.  The 'wi_stream' class corresponds
    to the special case (M,N)=(0,2) and the 'wi_transform' class corresponds to (M,N)=(2,2).
    
    (In more detail, there is a class hierarchy with a single abstract base class `pipeline_object`
    at the root.  The `wi_stream` and `wi_transform` classes are subclasses of pipeline_object,
    which are "semi-abstract" in the sense that they are more constrained than `pipeline_object`,
    but leave most of their behavior undefined, to be specified by a further level of subclassing.)

  - There are also container classes, most importantly the ubiquitous `pipeline` class, which contains
    a list of pipeline_objects, to be run sequentially (i.e. the output from each is the input
    to the next).  The `pipeline` class is a subclass of pipeline_object, so it can be used anywhere
    that a pipeline_object can (in particular, pipelines can contain pipelines).

  - In the simplest case where the pipeline is a single wi_stream followed by N wi_transforms, the
    syntax might look something like this:
    ```
    # An example of a wi_stream, see docstring for more info
    s = rf_pipelines.chime_dummy_network_stream(...)

    # Examples of wi_transforms, see docstrings for more info
    t1 = rf_pipelines.spline_detrender(...)
        ...
    tN = rf_pipelines.intensity_clipper(...)

    # Chain stream and transforms together to make a pipeline.
    p = rf_pipelines.pipeline([s,t1,...,tN])

    # Simplest way to run a pipeline, but we usually use rf_pipelines.run_for_web_viewer().
    p.run()
    ```

  - In general, a pipeline_object is runnable if the number of input ring buffers is zero.
    In the example above, this is the case because we started the pipeline `p` with a `wi_stream`,
    but it would not be the case if it started with a `wi_transform`.

  - Ring buffers may be downsampled (in time, or along other axes such as frequency) relative to the
    "top-level" pipeline sampling rate.

  - In particular, there is a 'wi_sub_pipeline' container class which is crucial for our fast RFI removal code.
    Like the `pipeline` container class, it contains an arbitrary sequence of pipeline_objects (probably wi_transforms),
    and the `wi_sub_pipeline` appears as a single pipeline_object which runs them sequentially.  However, the
    `wi_sub_pipeline` runs its "inner" sequence of transforms at lower (frequency, time) resolution.  This is
    done as follows:

       - The input to the `wi_sub_pipeline` is a pair of high-resolution (intensity, weights) arrays.
       - These arrays are downsampled to produce a pair of low-resolution (intensity, weights) arrays.
       - The "inner" sequence of transforms is run on the low-resolution arrays.  In particular, the
         low-resolution weights array now has extra masking applied.  (In general, we represent the
         RFI mask by zeroing elements of the weights array.)
       - The resulting low-resolution **mask** is upsampled back to high-resolution, and applied to
         the high-resolution weights array.  The low-resolution intensity array doesn't get upsampled
         and is "thrown away".

    Note that from the perspective of the high-resolution pipeline, the `wi_sub_pipeline` is a "clipper":
    it modifies the weights array by applying extra masking, and does not modify the intensity array.
    Under the hood, this is implemented by downsampling the data, running a full pipeline at low resolution
    to produce a low-res RFI mask, then upsampling this mask and applying it at high resolution.

  - For example, our CHIME RFI removal code currently works as follows.  We start with data with 16K
    frequency channels and 1ms sampling, i.e. (nfreq, dt) = (16K, 1ms).  We downsample by a factor of 16
    in frequency, i.e. to (nfreq, dt) = (1K, 1ms).  We then apply a long, complex sequence of transforms
    (around 100) and produce a 1K-channel RFI mask.  This mask is upsampled to the original 16K resolution 
    and applied.  At this point, the 16K-channel data has a complete RFI mask, but the intensities have
    not been detrended (or otherwise modified).  We apply two 16K-channel detrenders (one in the time
    direction and one in the frequency direction) to complete processing.

    This scheme is implemented as a "top-level" pipeline with three transforms: a wi_sub_pipeline and
    two detrenders.  Inside the wi_sub_pipeline, there is a long chain of around 100 transforms (detrenders
    and clippers) which operate at 1K-resolution.  This means that almost
    all of the transforms operate at 1K-resolution where they are 16 times faster.

    We plan to experiment with speeding things up further, by applying the "sub-pipeline" idea
    recursively.  Our (16K, 1ms) RFI removal is fast because we downsample to (1K, 1ms) and only
    have a few transforms which operate at the full (16K, 1ms) resolution.  Maybe we can also
    speed up the (1K, 1ms) RFI removal, by downsampling to an even lower resolution?

  - With the increased generality of the "new" v16 rf_pipelines, there are some new ideas for
    RFI removal that we plan to pursue soon.

      - A `mask_expander` which identifies regions in the (freq, time) plane which are mostly 
        RFI-masked, or where the RFI mask is slow to converge, and fully masks these regions.

        (It probably makes more sense to mask regions where the "delta-mask" is large, where
        the "delta-mask" is the difference between the current RFI mask, and the RFI mask at
        some previous point in the pipeline.  This means we need to "fork", or save a copy
        of the RFI mask for later use, which requires the v16 generality.)

      - A `kurtosis_filter` which compares the variance of the intensity to the square of its mean.
        If the electric field samples are Gaussian distributed (before squaring and averaging to get 
	intensity samples), then the variance is proportional to the square of the mean.  This
	gives a simple RFI-masking criterion which can be applied independently to tiny subsets 
	of the data.

        In the literature, this is sometimes called "kurtosis filtering" because it is often
        applied directly to the electric field samples (before squaring), in which case the
        comparison is between the kurtosis and variance (rather than variance and mean).
	
      - There is a lot of scope for internal simplifications by defining more ring buffers.
      
        For example, dedispersion should be implemented as a pipeline_object which takes
	(intensity, weights) buffers as inputs, and outputs one or more buffers containing
	coarse-grained triggers.  Further postprocessing of the triggers (e.g. plotting) can
	be done by defining new pipeline_objects which are placed later in the pipeline and
	operate on the coarse-grained trigger ring_buffers.

        In the old API, there was no way to define pipeline ring buffers for the coarse-grained
        trigger arrays.  As a result, our `bonsai_dedisperser` transform contains an
        ad hoc mini-pipeline which implements all possible postprocessing actions, since
        there is no way to pass them to the next transform.

        Similarly, the parts of our code related to variance estimation can be improved
        by defining ring buffers for the variance estimates.

  - Any object of type `rf_pipelines.pipeline_object` can now be serialized to a json file.
    This is a necessary step for running the pipeline with the `rfp-run` utility (see below),
    timing the pipeline with `rfp-time`, or using an RFI transform chain in the real-time
    L1 server ([https://github.com/kmsmith137/ch_frb_l1](https://github.com/kmsmith137/ch_frb_l1)).

    To write an rf_pipelines json file, the syntax is
    ```
    j = p.jsonize()    # where p is an object of type rf_pipelines.pipeline_object
    rf_pipelines.json_write(filename, j)
    ```
    To read an rf_pipelines json file, the syntax is
    ```
    j = rf_pipelines.json_read(filename)
    p = rf_pipelines.pipeline_object.from_json(j)   # returns an object of type rf_pipelines.pipeline_object
    ```

  - There are new command-line utilities which operate on json-seralized pipeline_objects, which we now describe in the next section.

<a name="command-line-utilities"></a>
### COMMAND-LINE UTILITIES

- **rfp-run:**  Runs a pipeline from the command line.  Can run in "batch mode", and use multiple cores in parallel.
  ```
  Usage: rfp-run [-nosh] [-w run_name] [-v verbosity] [-t nthreads] file1.json [file2.json file3.json ...]
  Each json file should contain either a jsonized pipeline_object, or a "run-list" of [suffix, json_filename] pairs
      -n: runs the pipeline with no output directory
      -o: show stdout during pipeline run (by default, stdout is suppressed, but stderr is shown)
      -s: throws exception unless single pipeline run (i.e. no run-lists allowed), infrequently used
      -w: runs the pipeline in a directory which is indexed by the web viewer (frb1 only)
      -v: specifies pipeline verbosity (integer, default 2)
      -t: number of threads (default 1, note that multiple threads are only useful if at least one json file is a run-list)
      -h: show longer help message and exit
  ```
  The rfp-run utility runs a pipeline from the command line.  The pipeline is constructed from a sequence of
  pipeline_objects which have previously been serialized with jsonize(), and specified on the command line.

  Alternatively, each command-line json file can be a "run-list" which points to a list of pipeline_object json files.
  In this case, multiple pipeline runs are performed in batch processing mode.

  A run-list file is just a json file containing a list of [suffix, json_filename] string pairs.  The run-list 'suffix'
  is a short string such as 'run1' which is appended to the web viewer run_name, so that we get a unique run_name for
  each run.  The run-list 'json_filename' is interpreted relative to the directory containing the run-list (not the
  current working directory).  The run-list file format is intended to be minimal enough that run-lists are easy to
  make by hand.  There are lots of examples of run-lists in [https://github.com/mrafieir/ch_frb_rfi](https://github.com/mrafieir/ch_frb_rfi).

  If the -w flag is specified, then the run will be viewable in the web viewer after it completes.  Exactly one of the -w,-n
  flags must be specified.

  If multiple threads are specified with -t NTHREADS, then multiple pipeline runs will be performed in parallel.  This only helps
  if run-lists are being used (so that there is more than one pipeline run to perform).  It often makes sense to set NTHREADS equal
  to the number of cores in the node.

  By default, the pipeline's stdout is not displayed to the screen, whereas stderr is.  Generally speaking, in rf_pipelines, we
  try to observe a convention where stderr is used to report warnings and unusual events, and stdout is used to report routine
  events.  If the -o flag is specified, then both stdout and stderr will be displayed.  Note that for a web viewer run,
  stdout always gets written to a log file, which can be viewed afterwards in the web viewer.

  The optional environment variable RF_PIPELINE_ATTRS contains additional pipeline attributes (specified as a json object of
  (key,value) pairs, serialized to a single string).  These attributes will be passed to the _bind() and _start_pipeline()
  methods of all pipeline_objects, and are also written to the pipeline's json output.

- **rfp-time.**  Time a pipeline from the command line.  Can run multiple pipeline instances in parallel on different cores,
  to emulate a "production" environment.
  ```
  Usage: rfp-time [-rP] [-t NTHREADS] [-j JSON_OUTFILE] file.json [file2.json file3.json ...]
     -t: change number of worker threads (default 1)
     -P: don't pin threads to cores (default is to pin threads, this should be done on an otherwise idle machine)
     -j: write json output from thread 0 to specified file (must not already exist)
  ```

- **rfp-analyze.**  Shows auxiliary info for a pipeline: ring buffer latencies, memory footprint.
  It could use more documentation, so the output may be cryptic, but we mention it here for completeness!
  ```
  Usage: rfp-analyze [-r] [-d DEPTH] [-j JSON_OUTFILE] file1.json [file2.json file3.json ...]
      -r: runs pipeline, and computes some extra information
      -d: limits depth of latency analysis (integer)
      -j: dumps result of pipeline_object.get_info() to a json file (usually for debugging)
  ```

- **rfp-json-show.**  Pretty-prints the contents of a json file (the script is wrapper around `rf_pipelines.json_show()`).
  ```
  Usage: rfp-json-show [-d DEPTH] file.json [key1 key2 ...]
     The -d flag expands output to specified depth (default 1)

     The "keys" are applied sequentially to the json object, before printing it.
     For example, if file.json contains an Object x whose "f" member is a list, then
        rfp-json-show.py file.json f 3
     will print x["f"][3] instead of printing x.
  ```

<a name="current-status"></a>
### CURRENT STATUS

The following table shows the current status, in the aftermath of the Mega Merge.
Some v15 features have not been ported yet to the v16 API, and there are some new
features which were not in v15 (for example, spline_detrender).

<table>
  <tr><th>Name</th> <th>Language</th> <th>Level of testing</th></tr>
  <tr> <th colspan="3" align="center">Containers</td> </tr>  
  <tr> <td>pipeline</td> <td>C++</td> <td>Fully tested</td> </tr>
  <tr> <td>wi_sub_pipeline</td> <td>C++/assembly</td> <td>Fully tested</td> </tr>
  <tr> <th colspan="3" align="center">Streams</td> </tr>
  <tr> <td>chime_stream_from_acqdir</td> <td>C++</td> <td>Fully tested</td>
  <tr> <td>chime_stream_from_filename</td> <td>C++</td> <td>Fully tested</td>
  <tr> <td>chime_stream_from_filename_list</td> <td>C++</td> <td>Fully tested</td>
  <tr> <td>chime_stream_from_times</td> <td>C++/python</td> <td>Fully tested</td>
  <tr> <td>chime_frb_stream_from_filename</td> <td>C++</td> <td>Not sure</td>
  <tr> <td>chime_frb_stream_from_filename_list</td> <td>C++</td> <td>Not sure</td>
  <tr> <td>chime_frb_stream_from_glob</td> <td>C++</td> <td>Not sure</td>
  <tr> <td>chime_network_stream</td> <td>C++</td> <td>Fully tested</td>
  <tr> <td>chime_dummy_network_stream</td> <td>C++</td> <td>Fully tested</td>
  <tr> <td>gaussian_noise_stream</td> <td>C++</td> <td>Needs unit test</td>
  <tr> <td>psrfits_stream</td> <td> -- </td> <td>Not ported from v15 yet</td>
  <tr> <th colspan="3" align="center">Detrenders</td> </tr>
  <tr> <td>spline_detrender</td> <td>C++/assembly</td> <td>Fully tested, but only AXIS_FREQ is implemented</td>
  <tr> <td>polynomial_detrender</td> <td>C++/assembly</td> <td>Fully tested</td>
  <tr> <th colspan="3" align="center">Clippers</td> </tr>
  <tr> <td>intensity_clipper</td> <td>C++/assembly</td> <td>Fully tested</td>
  <tr> <td>std_dev_clipper</td> <td>C++/assembly</td> <td>Fully tested</td>
  <tr> <td>mask_expander</td> <td>C++ (poorly optimized)</td> <td>Half-finished, untested</td>
  <tr> <th colspan="3" align="center">CHIME-specific</td> </tr>
  <tr> <td>chime_file_writer</td> <td>C++</td> <td>Fully tested</td>
  <tr> <td>chime_packetizer</td> <td>C++</td> <td>Untested since porting from v15</td>
  <tr> <td>chime_16k_spike_mask</td> <td>C++</td> <td>Fully tested</td>
  <tr> <td>chime_16k_derippler</td> <td>C++</td> <td>Fully tested</td>
  <tr> <td>chime_16k_stripe_analyzer</td> <td>C++</td> <td>Fully tested</td>
  <tr> <th colspan="3" align="center">Miscellaneous</td> </tr>
  <tr> <td>adversarial_masker</td> <td>Python</td> <td>Untested since porting from v15</td>
  <tr> <td>badchannel_mask</td> <td>C++</td> <td>Fully tested</td>
  <tr> <td>bonsai_dedisperser</td> <td>Python and C++ (*)</td> <td>Partially tested since porting from v15</td>
  <tr> <td>frb_injector_transform</td> <td>Python (**)</td> <td>Untested since porting from v15</td>
  <tr> <td>mask_filler</td> <td>Python</td> <td>Untested since porting from v15</td>
  <tr> <td>noise_filler</td> <td>Python</td> <td>Untested since porting from v15</td>
  <tr> <td>online_mask_filler</td> <td> -- </td> <td>Not ported from v15 yet</td>
  <tr> <td>pipeline_fork</td> <td>C++</td> <td>Untested</td>
  <tr> <td>plotter_transform</td> <td>Python</td> <td>Fully tested</td>
  <tr> <td>reverter</td> <td> -- </td> <td>Not ported from v15 yet</td>
  <tr> <td>rfi_bitmask</td> <td> -- </td> <td>Not ported from v15 yet</td>
  <tr> <td>spectrum_analyzer</td> <td>C++</td> <td>Fully tested</td>
  <tr> <td>variance_estimator</td> <td>Python</td> <td>Untested since porting from v15</td>
</table>

(*) We currently have two versions of the bonsai_dedisperser, one written in C++ and one written in python,
    with tradeoffs between them (each one has features that the other is missing).  The long-term plan is
    to consolidate into a "grand unified bonsai_dedisperser" written in C++.

(**) Similarly, the frb_injector transform used to have two versions, one written in C++ and one written in python,
     but only the python version has been ported from v15.

<a name="class-hierarchy"></a>
### CLASS HIERARCHY

- pipeline_object
  - pipeline
     - wi_sub_pipeline
  - pipeline_fork
  - chunked_pipeline_object
      - chime_16k_spike_mask
      - chime_16k_derippler
      - chime_16k_stripe_analyzer
      - mask_expander
      - wi_stream
          - chime_file_stream_base
	      - chime_file_stream
	      - chime_frb_file_stream
          - chime_network_stream
	  - chime_dummy_network_stream
	  - gaussian_noise_stream
      - wi_transform
          - badchannel_mask
	  - bonsai_dedisperser (has both C++ and python versions)
	  - chime_file_writer
	  - chime_packetizer
	  - intensity_clipper
	  - polynomial_detrender
	  - spectrum_analyzer
	  - spline_detrender
	  - std_dev_clipper
	  - adversarial_masker (python)
	  - bonsai_dedisperser (python)
	  - frb_injector_transform (python)
	  - mask_filler (python)
	  - plotter_transform (python)
	  - variance_estimator (python)

