#include "rf_pipelines_internals.hpp"
#include <fstream>

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


struct badchannel_mask : public wi_transform {
    // Note: inherits { nfreq, nt_chunk, nt_prepad, nt_postpad } from base class wi_transform

    vector<double> m_bad_channels; // holds the bad frequencies, as specified from the input file
    vector<int> m_bad_indices;     // holds the final bad indices to be used to mask the weights array
    int m_len_indices;             // the size of bad_indices vector
   
    badchannel_mask(const string &maskpath, int nt_chunk)
    {
        stringstream ss;
	ss << "badchannel_mask_cpp(maskpath=" << maskpath << ", nt_chunk=" << nt_chunk << ")";
        this->name = ss.str();
	this->nt_prepad = 0;
	this->nt_postpad = 0;
	this->nt_chunk = nt_chunk;
	
	rf_assert(nt_chunk > 0);

        // Extract the channels to be removed into m_bad_channels
	get_bad_channels(maskpath, m_bad_channels);
    }

    void get_bad_channels(const string &maskpath, vector<double> &bad_channels)
    {
        ifstream inf(maskpath);

	if (inf.is_open())
	{
	    string line, freq;
	    double low, high;
	    while (getline(inf, line))
	    {
	        stringstream linestream(line);

		if (line[0] != '#')
		{
		    getline(linestream, freq, ',');
		    low = atof(freq.c_str());
		    getline(linestream, freq, ',');
		    high = atof(freq.c_str());
		    rf_assert(low < high && "The lower frequency must be less than the high frequency!");

		    bad_channels.push_back(low);
		    bad_channels.push_back(high);
		}
	    }
	    inf.close();
	}
	else
	    throw runtime_error("badchannel_mask: couldn't open file at the maskpath given!");
    }

    // As explaned in rf_pipelines.hpp, the following four virtual functions in the base class
    // must be overridden, in order to define the badchannel_mask subclass.
    virtual void set_stream(const wi_stream &stream) override
    {
        this->nfreq = stream.nfreq;
	double freq_lo_MHz = stream.freq_lo_MHz;
	double freq_hi_MHz = stream.freq_hi_MHz;

	// First, make sure all intervals are within the freq
	int vec_size = m_bad_channels.size();
	double low, high;
	vector<double> temp;
	for (int ilo=0; ilo < vec_size-1; ilo+=2)
        {
	    low = m_bad_channels[ilo];
	    high = m_bad_channels[ilo + 1];
	    // Both ends of the range are within the observation range
	    if (low > freq_lo_MHz && high < freq_hi_MHz)
	    {
	        temp.push_back(low);
	        temp.push_back(high);
	    }
	    // The low end is below the observation range
	    else if (low < freq_lo_MHz && high > freq_lo_MHz && high < freq_hi_MHz)
	    {
	        temp.push_back(freq_lo_MHz);
	        temp.push_back(high);
	    }
	    // The high end is above the observation range
	    else if (high > freq_hi_MHz && low < freq_hi_MHz && low > freq_lo_MHz)
	    {
	        temp.push_back(low);
	        temp.push_back(freq_hi_MHz);
	    }
	}

	// Here, we rescale the bad frequencies into bad indices. Adapted 
	// from Masoud's python transform comments...
	// First we need to scale our frequency mask so that it matches with
        // the number of channels in the chunk. Subtracting the max value 
	// from the mask (which runs from low to high values) leaves us
	// with an array that runs in reverse. We make sure that we include
	// any non-zero overlap with the mask. To this end, we take the
        // ceiling of the first element -- remember this is in a reversed
	// order, which means, e.g., a[0,1] has become a[1,0] -- and the
	// floor of the second element. Next, we convert them to integers
	// so they can be fed as indexes into the weights array.
	// We take the max at the end to ensure no values are below 0.
	double scale = stream.nfreq / (freq_hi_MHz - freq_lo_MHz);
	double factor = scale * freq_hi_MHz;
	m_len_indices = temp.size();
	for (int i = 0; i < m_len_indices; ++i)
	{
	    if (i % 2 == 0)
	      m_bad_indices.push_back(min(max(int(ceil(factor - temp[i]*scale)), 0), (int) stream.nfreq-1));
	    else
	      m_bad_indices.push_back(min(max(int(floor(factor - temp[i]*scale)), 0), (int) stream.nfreq-1));
	}
    }

    virtual  void start_substream(int isubstream, double t0) override
    {
        // Nothing necessary here
    }

    virtual void process_chunk(double t0, double t1, float *intensity, float *weights, ssize_t stride, float *pp_intensity, float *pp_weights, ssize_t pp_stride) override
    {
         // Loop over bad indices list
        for (int ibad_index=0; ibad_index < m_len_indices; ibad_index+=2)
        {
	    // Iterate over each frequency in the bad index range
	    for (int ifreq=m_bad_indices[ibad_index+1]; ifreq <= m_bad_indices[ibad_index]; ++ifreq)
	    {
	        // Set all weights in the channel to 0
	        for (int it=0; it < nt_chunk; ++it)
		    weights[ifreq*stride+it] = 0;
	    }
        }
    }

    virtual void end_substream() override
    {
        // Nothing necessary here
    }
};


shared_ptr<wi_transform> make_badchannel_mask(const string &maskpath, int nt_chunk)
{
    return make_shared<badchannel_mask> (maskpath, nt_chunk);
}

}  // namespace rf_pipelines
