### RF_PIPELINES

rf_pipelines: A plugin-based framework for processing channelized intensity data in radio astronomy.

Note: this repo now includes the web viewer code which was previously
in [mburhanpurkar/web_viewer](https://github.com/mburhanpurkar/web_viewer).

This manual is mostly a placeholder, and just contains a few random sections,
rather than systematic documentation.  Fixing this is a high priority!

### CONTENTS

  - [Current status](#user-content-current-status)
  - [Class hierarchy](#user-content-class-hierarchy)

<a name="current-status"></a>
### CURRENT STATUS

In v16, the low-level logic in rf_pipelines was largely rewritten!  The
low-level building blocks (streams, transforms, etc.) have a different API now,
and any "pre-v16" code will need substantial changes to work "post-v16".
This is part of the "2017 Mega Merge" affecting many parts of the CHIMEFRB pipeline.

To make matters worse, due to the CHIME pre-deployment scramble, the new API
is fragile and incomplete, and not all pre-v16 features are working.  We hope
that everything will be finished and well-tested in a few more revisions.

In the meantime, the following table shows the current status.  Some v15
features have not been ported yet to the new API, and there are some new
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

