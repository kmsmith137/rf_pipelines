#include "rf_pipelines_internals.hpp"
#include <fstream>

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


struct badchannel_mask : public wi_transform {
    // Note: inherits { nfreq, nt_chunk, nds } from base class wi_transform.

    // Constructor arguments
    const string mask_path;
    const vector<pair<double,double>> mask_ranges;

    vector<double> m_bad_channels; // holds the bad frequencies, as specified from the input file
    vector<int> m_bad_indices;     // holds the final bad indices to be used to mask the weights array
    int m_len_indices;             // the size of bad_indices vector
   
    badchannel_mask(const string &mask_path_, const vector<pair<double,double>> &mask_ranges_) :
	wi_transform("badchannel_mask"),
	mask_path(mask_path_),
	mask_ranges(mask_ranges_)
    {
        stringstream ss;
	ss << "badchannel_mask(mask_path=\"" << mask_path << "\"";
	if (mask_ranges.size() > 0)
	    ss << ", mask_ranges=<" << mask_ranges.size() << ">";
	ss << ")";

        this->name = ss.str();

        // Extract the channels to be removed into m_bad_channels
	if (mask_path.size() > 0)
	    get_bad_channels(mask_path, m_bad_channels);

	// Add mask_ranges.
	for (const auto &p: mask_ranges) {
	    double freq_lo = p.first;
	    double freq_hi = p.second;

	    if (freq_lo >= freq_hi)
		throw runtime_error("rf_pipelines::badchannel_mask constructor: expected freq_lo < freq_hi in 'mask_ranges'");

	    m_bad_channels.push_back(freq_lo);
	    m_bad_channels.push_back(freq_hi);
	}
    }

    // Helper function called by constructor.
    void get_bad_channels(const string &mask_path, vector<double> &bad_channels)
    {
        ifstream inf(mask_path);

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
	    throw runtime_error("badchannel_mask: couldn't open file at the mask_path given!");
    }

    // Called after 'nfreq' is initialized.
    virtual void _bind_transform(Json::Value &json_attrs) override
    {
	if (!json_attrs.isMember("freq_lo_MHz") || !json_attrs.isMember("freq_hi_MHz"))
	    throw runtime_error("badchannel_mask: expected json_attrs to contain members 'freq_lo_MHz' and 'freq_hi_MHz'");

	double freq_lo_MHz = json_attrs["freq_lo_MHz"].asDouble();
	double freq_hi_MHz = json_attrs["freq_hi_MHz"].asDouble();

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
	double scale = nfreq / (freq_hi_MHz - freq_lo_MHz);
	double factor = scale * freq_hi_MHz;
	m_len_indices = temp.size();
	for (int i = 0; i < m_len_indices; ++i)
	{
	    if (i % 2 == 0)
	      m_bad_indices.push_back(min(max(int(ceil(factor - temp[i]*scale)), 0), (int) nfreq-1));
	    else
	      m_bad_indices.push_back(min(max(int(floor(factor - temp[i]*scale)), 0), (int) nfreq-1));
	}
    }

    virtual void _process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) override
    {
         // Loop over bad indices list
        for (int ibad_index=0; ibad_index < m_len_indices; ibad_index+=2)
        {
	    // Iterate over each frequency in the bad index range
	    for (int ifreq=m_bad_indices[ibad_index+1]; ifreq <= m_bad_indices[ibad_index]; ++ifreq)
	    {
	        // Set all weights in the channel to 0
	        for (int it=0; it < nt_chunk; ++it)
		    weights[ifreq*wstride+it] = 0;
	    }
        }
    }

    virtual Json::Value jsonize() const override
    {
	Json::Value ret;
	ret["class_name"] = "badchannel_mask";
	ret["mask_path"] = mask_path;
	ret["mask_ranges"] = Json::Value(Json::arrayValue);

	for (const auto &p: mask_ranges) {
	    Json::Value v;
	    v.append(p.first);
	    v.append(p.second);
	    ret["mask_ranges"].append(v);
	}

	return ret;
    }

    static shared_ptr<badchannel_mask> from_json(const Json::Value &j)
    {
	string mask_path = string_from_json(j, "mask_path");
	vector<pair<double,double>> mask_ranges;

	Json::Value a = array_from_json(j, "mask_ranges");

	for (const Json::Value &v: a) {
	    if (!v.isArray() || (v.size() != 2) || !v[0].isDouble() || !v[1].isDouble())
		throw runtime_error("rf_pipelines::badchannel_mask::from_json: expected each element of 'mask_ranges' to be a length-2 floating-point array");
	    mask_ranges.push_back(make_pair(v[0].asDouble(), v[1].asDouble()));
	}
	
	return make_shared<badchannel_mask> (mask_path, mask_ranges);
    }
};


shared_ptr<wi_transform> make_badchannel_mask(const string &mask_path, const vector<pair<double,double>> &mask_ranges)
{
    return make_shared<badchannel_mask> (mask_path, mask_ranges);
}


namespace {
    struct _init {
	_init() {
	    pipeline_object::register_json_deserializer("badchannel_mask", badchannel_mask::from_json);
	}
    } init;
}


}  // namespace rf_pipelines
