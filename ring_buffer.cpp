#include "rf_pipelines_internals.hpp"

// Set to 1 to enable some very verbose debugging output
#define RF_RB_DEBUG 0

using namespace std;

namespace rf_pipelines {
#if 0
}  // emacs pacifier
#endif


// Only used for debugging.
static std::random_device g_rd;
static std::mt19937 g_rng(g_rd());


ring_buffer::ring_buffer(const vector<ssize_t> &cdims_, ssize_t nds_, bool debug_, const string &name_) :
    cdims(cdims_), 
    csize(prod(cdims_)),
    nds(nds_),
    debug(debug_),
    name(name_)
{
    check_cdims(cdims);

    if (nds <= 0)
	throw runtime_error("rf_pipelines::ring_buffer: expected nds > 0");
}


ring_buffer::~ring_buffer()
{
    free(buf);
    buf = nullptr;
}


void ring_buffer::update_params(ssize_t nt_contig_, ssize_t nt_maxlag_)
{
    rf_assert(buf == nullptr);
    rf_assert(nt_contig_ > 0);
    rf_assert(nt_maxlag_ >= nt_contig_);
    
    this->nt_contig = max(nt_contig, nt_contig_);
    this->nt_maxlag = max(nt_maxlag, nt_maxlag_);
}


void ring_buffer::allocate()
{
    rf_assert(nt_contig > 0);
    rf_assert(nt_maxlag >= nt_contig);
    rf_assert(ap == nullptr);

    // Double call to allocate() is not an error.
    if (buf != nullptr)
	return;

    this->_preallocate();

    // This strengthens unit tests a bit.
    if (debug)
	stride += 32 * randint(g_rng, 0, 8);
    
    this->buf = aligned_alloc<float> (csize * stride);

#if RF_RB_DEBUG
    cout << "ring_buffer::allocate(): nds=" << nds << ", nt_contig=" << nt_contig << ", nt_maxlag=" << nt_maxlag
	 << ", period=" << period << ", stride=" << stride << endl;
#endif
}


// _preallocate(): helper function which assigns 'period' and 'stride'.
// This used to be part of allocate(), but factored out so that it can be called from get_info().
// This is useful because get_info() can be called on an unallocated ring_buffer.

void ring_buffer::_preallocate()
{
    rf_assert(nt_contig > 0);
    rf_assert(nt_maxlag >= nt_contig);

    // The memory alignment heuristics below are intended to improve L1 cache associativity.
    // FIXME: some day, define a boolean flag for toggling this, to see how much it actually helps!

    this->period = (nt_maxlag + nds - 1) / nds;
    this->period = round_up(period, 32);

    this->stride = period + (nt_contig + nds - 2) / nds;
    this->stride = round_up(stride, 16);

    if (stride % 32 == 0)
	stride += 16;
}


void ring_buffer::deallocate()
{
    rf_assert(ap == nullptr);

    free(buf);
    buf = nullptr;
}


void ring_buffer::reset()
{
    this->curr_pos = 0;
    this->first_valid_sample = 0;
    this->last_valid_sample = 0;
    this->high_water_mark = 0;
    this->optimal_period = 0;
    this->nget_tot = 0;
    this->nget_mirror = 0;
}


float *ring_buffer::get(ssize_t pos0, ssize_t pos1, int mode)
{
#if RF_RB_DEBUG
    cout << "ring_buffer::get(" << access_mode_to_string(mode) << "): pos=(" << pos0 << "," << pos1
	 << "), valid=(" << first_valid_sample << "," << last_valid_sample << ")" << endl;
#endif

    // Argument checking
    rf_assert(pos0 >= 0);
    rf_assert(pos0 <= pos1);
    rf_assert(pos1 - pos0 <= nt_contig);
    rf_assert(pos0 % nds == 0);
    rf_assert(pos1 % nds == 0);
    rf_assert(mode != ACCESS_NONE);
    rf_assert(buf != nullptr);
    rf_assert(ap == nullptr);

    // Set ap_pos* before applying downsampling factor.
    // (Remaining fields 'ap' and 'ap_mode' will be set later.)
    this->ap_pos0 = pos0;
    this->ap_pos1 = pos1;

    // Apply downsampling factor
    pos0 /= nds;
    pos1 /= nds;

    if (mode == ACCESS_APPEND) {
	// Range check and advance buffer
	rf_assert(pos0 == curr_pos);
	curr_pos = pos1;
    }
    else {
	// Range check
	rf_assert(pos0 >= curr_pos - period);
	rf_assert(pos1 <= curr_pos);
    }

    // Sample range in memory
    ssize_t it0 = pos0 % period;
    ssize_t it1 = it0 + (pos1 - pos0);
    
    // Mirror data if necessary
    if (mode & ACCESS_READ) {
	_mirror_initial(it0);
	_mirror_final(it1);
    }
    else
	_mirror_initial(it1);

    this->ap = this->buf + it0;
    this->ap_mode = mode;
    this->high_water_mark = max(high_water_mark, it1);
    this->optimal_period = max(optimal_period, curr_pos - pos0);
    this->nget_tot += (it1 - it0);   // note: nget_mirror is updated in ring_buffer::_copy()

    return ap;
}


void ring_buffer::put(float *p, ssize_t pos0, ssize_t pos1, int mode)
{
#if RF_RB_DEBUG
    cout << "ring_buffer::put(" << access_mode_to_string(mode) << "): pos=(" << pos0 << "," << pos1
	 << "), valid=(" << first_valid_sample << "," << last_valid_sample << ")" << endl;
#endif

    rf_assert(ap == p);
    rf_assert(ap_pos0 == pos0);
    rf_assert(ap_pos1 == pos1);
    rf_assert(ap_mode == mode);

    this->ap = nullptr;

    if (!(mode & ACCESS_WRITE))
	return;

    // Cut-and-paste logic from ring_buffer::get(), for determining (it0,it1).
    pos0 /= nds;
    pos1 /= nds;
    ssize_t it0 = pos0 % period;
    ssize_t it1 = it0 + (pos1 - pos0);

#if RF_RB_DEBUG
    ssize_t save_first_valid = first_valid_sample;
    ssize_t save_last_valid = last_valid_sample;
#endif
    
    if (it0 < first_valid_sample) {
	rf_assert(first_valid_sample <= it1);
	first_valid_sample = it0;
    }
	
    if (it1 > last_valid_sample) {
	rf_assert(last_valid_sample >= it0);
	last_valid_sample = it1;
    }

    last_valid_sample = min(last_valid_sample, it0 + period);
    first_valid_sample = max(first_valid_sample, it1 - period);

#if RF_RB_DEBUG
    if ((first_valid_sample != save_first_valid) || (last_valid_sample != save_last_valid)) {
	cout << "    update valid: (" << save_first_valid << "," << save_last_valid
	     << ") -> (" << first_valid_sample << "," << last_valid_sample << ")" << endl;
    }
#endif
}


ssize_t ring_buffer::get_stride() const
{
    rf_assert(buf != nullptr);
    return stride;
}


Json::Value ring_buffer::get_info()
{
    this->_preallocate();

    Json::Value j;
    j["name"] = name;
    j["cdims"] = Json::Value(Json::arrayValue);
    j["csize"] = Json::Int64(csize);
    j["nds"] = Json::Int64(nds);
    j["nt_contig"] = Json::Int64(nt_contig);
    j["nt_maxlag"] = Json::Int64(nt_maxlag);
    j["period"] = Json::Int64(period);
    j["stride"] = Json::Int64(stride);
    j["high_water_mark"] = Json::Int64(high_water_mark);
    j["optimal_period"] = Json::Int64(optimal_period);
    j["nget_tot"] = Json::Int64(nget_tot);
    j["nget_mirror"] = Json::Int64(nget_mirror);
    j["mb"] = 4.0e-6 * double(stride) * double(csize);

    for (ssize_t d: cdims)
	j["cdims"].append(Json::Int64(d));

    return j;
}


void ring_buffer::_mirror_initial(ssize_t it0)
{
    if (it0 < first_valid_sample) 
    {
#if RF_RB_DEBUG
    cout << "    _mirror_initial: valid=(" << first_valid_sample << "," << last_valid_sample
	 << ") -> (" << it0 << "," << last_valid_sample << ")" << endl;
#endif

	rf_assert(last_valid_sample >= first_valid_sample + period);
	_copy(it0, it0 + period, first_valid_sample - it0);
	first_valid_sample = it0;
    }
}


void ring_buffer::_mirror_final(ssize_t it1)
{
    if (it1 > last_valid_sample) 
    {
#if RF_RB_DEBUG
    cout << "    _mirror_final: valid=(" << first_valid_sample << "," << last_valid_sample
	 << ") -> (" << first_valid_sample << "," << it1 << ")" << endl;
#endif

	rf_assert(first_valid_sample <= last_valid_sample - period);
	_copy(last_valid_sample, last_valid_sample - period, it1 - last_valid_sample);
	last_valid_sample = it1;
    }
}


void ring_buffer::_copy(ssize_t it_dst, ssize_t it_src, ssize_t n)
{
#if RF_RB_DEBUG
    cout << "    _copy: dst=" << it_dst << ", src=" << it_src << ", n=" << n << endl;
#endif

    for (ssize_t i = 0; i < csize; i++)
	memcpy(buf + i*stride + it_dst, buf + i*stride + it_src, n * sizeof(float));
    
    this->nget_mirror += n;
}


// static member function
string ring_buffer::access_mode_to_string(int access_mode)
{
    if (access_mode == ACCESS_NONE)
	return "ACCESS_NONE";
    if (access_mode == ACCESS_READ)
	return "ACCESS_READ";
    if (access_mode == ACCESS_WRITE)
	return "ACCESS_WRITE";
    if (access_mode == ACCESS_RW)
	return "ACCESS_RW";
    if (access_mode == ACCESS_APPEND)
	return "ACCESS_APPEND";

    throw runtime_error("rf_pipelines: internal error: bad argument to ring_buffer::access_mode_to_string()");
}


// static member functon
void ring_buffer::check_cdims(const vector<ssize_t> &cdims)
{
    if (cdims.size() >= 6)
	throw runtime_error("rf_pipelines: attempt to construct high-dimensional ring buffer is probably unintentional");

    for (size_t i = 0; i < cdims.size(); i++) {
	if (cdims[i] <= 0)
	    throw runtime_error("rf_pipelines::ring_buffer: expected all dimensions > 0");
    }
}


}  // namespace rf_pipelines
