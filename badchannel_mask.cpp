#include "rf_pipelines_internals.hpp"
#include <cassert>
#include <fstream>

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


// This file is a placeholder for a C++ implementation of the 'badchannel_mask' class.
//
// Right now, it just contains some general boilerplate needed to make a C++ badchannel_mask
// class available throughout the pipeline.  In particular, the C++ badchannel_mask can be
// obtained from python as 'rf_pipelines.rf_pipelines_c.make_badchannel_mask()'.
//
// The actual implementation of the badchannel_mask class is still missing -- right now,
// all member functions throw exceptions.

struct badchannel_mask : public wi_transform {
    // Note: inherits { nfreq, nt_chunk, nt_prepad, nt_postpad } from base class wi_transform

    badchannel_mask(const string &maskpath, int nt_chunk)
    {
        stringstream ss;
	ss << "badchannel_mask_cpp(maskpath=" << maskpath << ", nt_chunk=" << nt_chunk << ")";
        this->name = ss.str();

        // Extract the channels to be removed into bad_channels
        vector<float> bad_channels; 
	this->_get_bad_channels(maskpath, bad_channels);
    }

    void _get_bad_channels(const string &maskpath, vector<float> &bad_channels)
    {
        ifstream inf(maskpath);
        string line;
        string freq;
        int length = 0;
        float low, high;

        while (getline(inf,line))
        {
            stringstream linestream(line);

            if (line[0] != '#')
            {
	        if (length == 0)
	        {
	            // Get the number of frequency pairs that will be masked 
	            length = 2*atoi(line.c_str());
	            continue;
	        }

	        getline(linestream, freq, ',');
	        low = atof(freq.c_str());
	        getline(linestream, freq, ',');
	        high = atof(freq.c_str());
	        assert (low < high && "The lower frequency must be less than the high frequency!");

	        bad_channels.push_back(low);
	        bad_channels.push_back(high);
            }
        }

        inf.close();
    
        cout << "Frequency pairs extracted: " << endl;
        for (int i = 0; i < length; i += 2) cout << bad_channels[i] << " " << bad_channels[i + 1] << endl;
    }

    // As explaned in rf_pipelines.hpp, the following four virtual functions in the base class
    // must be overridden, in order to define the badchannel_mask subclass.

    virtual void set_stream(const wi_stream &stream) override
    {
        this->nfreq = stream.nfreq;
    }

    virtual void start_substream(int isubstream, double t0) override
    {

    }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
	throw runtime_error("oops, C++ badchannel_mask class was requested, but not implemented yet");
    }

    virtual void end_substream() override
    {
	throw runtime_error("oops, C++ badchannel_mask class was requested, but not implemented yet");
    }
};


shared_ptr<wi_transform> make_badchannel_mask(const string &maskpath, int nt_chunk)
{
    return make_shared<badchannel_mask> (maskpath, nt_chunk);
}


}  // namespace rf_pipelines
