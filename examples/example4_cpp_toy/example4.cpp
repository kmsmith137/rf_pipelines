#include <rf_pipelines.hpp>

using namespace std;
using namespace rf_pipelines;


int main(int argc, char **argv)
{
    // An rf_pipeline consists of a "stream" object which generates chunks of weighted intensity data,
    // followed by a sequence of "transforms" which operate on chunks of weighted intensity data.

    // First let's make the stream.  There are several experiment-specific streams defined in rf_pipelines.hpp, 
    // but we'll just use a simple stream which outputs Gaussian random numbers.

    int nfreq = 1024;               // Number of frequency channels
    int nt_tot = 32000;             // Number of time samples before end of stream
    double freq_lo_MHz = 400.;      // CHIME
    double freq_hi_MHz = 800.;      // CHIME
    double dt_sample = 1.31072e-3;  // CHIME pathfinder

    shared_ptr<wi_stream> stream = make_gaussian_noise_stream(nfreq, nt_tot, freq_lo_MHz, freq_hi_MHz, dt_sample);

    // Next let's make some transforms.  See rf_pipelines.hpp for a complete list.
    vector< shared_ptr<wi_transform> > transform_list;
    
    // First, a simple detrending transform which is probably too simple for the real world!
    transform_list.push_back( make_simple_detrender(1024) );
    
    // In a real pipeline, we'd put some RFI-removing transforms next, but these aren't implemented in C++ yet
    // (coming soon, let me know if you'd like to help) so this comment is a placeholder.

    // For debugging it's sometimes useful to "see" what's going on by plotting the intensity stream
    // somewhere in the pipeline.  Currently the best way of doing this from C++ is to insert a transform
    // which writes data to disk in CHIME hdf5 format (this is a "pseudo-transform", meaning that it doesn't
    // actually modify its input).  You can then use the utility 'ch-plot-intensity-file' in the ch_frb_io
    // github repo to make a waterfall plot.
    
#if 1  // uncomment to enable
    string filename = "example4_data.hdf5";
    bool clobber = true;
    transform_list.push_back( make_chime_file_writer(filename,clobber) );
#endif

    //
    // Now add a 'bonsai_transform', a pseudo-transform which dedisperses the data.
    //
    // Before this will work, you'll need to create the file 'bonsai_config.hdf5' from the
    // configuration file 'bonsai_config.txt', using the command:
    //    bonsai-mkweight bonsai_config.txt bonsai_config.hdf5
    //
    string bonsai_config_filename = "bonsai_config.hdf5";
    string bonsai_output_filename = "bonsai_outputs.hdf5";
    transform_list.push_back( make_bonsai_dedisperser(bonsai_config_filename, bonsai_output_filename) );
    
    //
    // Now that the stream and transforms have been set up, this line of code runs the pipeline!
    //
    // To plot the output of the dedispersion transform, use the command:
    //    bonsai-plot-triggers.py bonsai_outputs.hdf5
    //
    // This actually generates three plots (bonsai_outputs_treeN.png, where N=0,1,2) since the bonsai
    // config file defines three dedispersion trees correpsonding to different DM and pulse width ranges.

    stream->run(transform_list);

    return 0;
}
