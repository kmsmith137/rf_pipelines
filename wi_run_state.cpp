#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


wi_run_state::wi_run_state(const wi_stream &stream, const std::vector<std::shared_ptr<wi_transform> > &transforms_) :
    nfreq(stream.nfreq),
    freq_lo_MHz(stream.freq_lo_MHz),
    freq_hi_MHz(stream.freq_hi_MHz),
    dt_sample(stream.dt_sample),
    nt_stream_maxwrite(stream.nt_maxwrite),
    ntransforms(transforms_.size()),
    transforms(transforms_),
    transform_ipos(transforms_.size(), 0),
    stream_ipos(0),
    prepad_buffers(transforms_.size()),
    is_running(false)
{
    if (!nfreq)
	throw runtime_error("rf_pipelines: run() called on uninitialized stream (did you forget to call wi_stream::construct()?)");
    if (ntransforms < 1)
	throw runtime_error("rf_pipelines: run() called on empty transform list");
}


void wi_run_state::start_stream()
{
    if (this->is_running)
	throw runtime_error("wi_run_state::start_stream() called on already-running stream (maybe a call to end_stream() is missing somewhere?)");
   
    for (int it = 0; it < ntransforms; it++) {
	transforms[it]->start_stream(nfreq, freq_lo_MHz, freq_hi_MHz, dt_sample);
	transforms[it]->check_invariants();
    }

    // Allocate main buffer

    int nt_contig = nt_stream_maxwrite;
    for (int it = 0; it < ntransforms; it++)
	nt_contig = max(nt_contig, transforms[it]->nt_chunk + transforms[it]->nt_postpad);

    int nt_ring = transforms[0]->nt_chunk + transforms[0]->nt_postpad + nt_stream_maxwrite;

    for (int it = 1; it < ntransforms; it++) {
	int g = transforms[it-1]->nt_chunk;
	g = gcd(g, transforms[it]->nt_chunk);
	g = gcd(g, transforms[it]->nt_postpad);
	nt_ring += transforms[it]->nt_chunk + transforms[it]->nt_postpad - g;
    }
    
    this->main_buffer.construct(nfreq, nt_contig, nt_ring);

    //
    // Allocate prepad buffers
    //
    // Note that we shift the prepad buffer by nt_prepad, i.e. "time" i in the prepad buffer 
    // is actually a saved sample from time (i-nt_prepad).
    //
    // FIXME wraparaound_buf as currently implemented is a little suboptimal here.  I'm leaving this
    // as something to fix in the future since I doubt the suboptimality is very important.
    //
    for (int it = 0; it < ntransforms; it++) {
	int n0 = transforms[it]->nt_prepad;
	int n1 = transforms[it]->nt_chunk;

	if (!n0)
	    continue;

	int nt_contig = max(n0, n1);
	int nt_ring = n0 + n1;
	this->prepad_buffers[it].construct(nfreq, nt_contig, nt_ring);

	// shift prepad buffer by n0 as noted above
	this->prepad_buffers[it].append_zeros(n0);
    }

    this->is_running = true;
}


void wi_run_state::setup_write(int nt, float* &intensityp, float* &weightp, int &stride, bool zero_flag)
{
    if (!this->is_running)
	throw runtime_error("wi_run_state::setup_write() called on non-running stream (did you forget a call to wi_run_state::start_stream()?)");
    if (nt <= 0)
	throw runtime_error("wi_run_state::setup_write(): expected nt > 0");
    if (nt > this->nt_stream_maxwrite)
	throw runtime_error("wi_run_state::setup_write(): expected nt <= nt_stream_maxwrite");

    this->main_buffer.setup_append(nt, intensityp, weightp, stride, zero_flag);
}


void wi_run_state::finalize_write(int nt)
{
    if (!this->is_running)
	throw runtime_error("wi_run_state::finalize_write() called on non-running stream");

    this->main_buffer.finalize_append(nt);
    this->stream_ipos += nt;

    int curr_ipos = stream_ipos;

    for (int it = 0; it < ntransforms; it++) {
	int n0 = transforms[it]->nt_prepad;
	int n1 = transforms[it]->nt_chunk;
	int n2 = transforms[it]->nt_postpad;

	while (transform_ipos[it] + n1 + n2 <= curr_ipos) {

	    // Logic for running transform is here

	    float *intensity = nullptr;
	    float *weights = nullptr;
	    int stride = 0;

	    float *pp_intensity = nullptr;
	    float *pp_weights = nullptr;
	    int pp_stride = 0;

	    // Note (n1+n2) here, versus (n1) in call to finalize_write() below.
	    main_buffer.setup_write(transform_ipos[it], n1+n2, intensity, weights, stride);
	    
	    if (n0 > 0) {
		// Transform is prepadded.  We need to get pointers to the prepadded data, and
		// we also need to save data which will be overwritten by running the transform,
		// since it will be prepadded data in a future iteration.

		bool zero_flag = false;
		prepad_buffers[it].setup_append(n1, pp_intensity, pp_weights, pp_stride, zero_flag);

		for (int ifreq = 0; ifreq < nfreq; ifreq++) {
		    memcpy(pp_intensity + ifreq*pp_stride, intensity + ifreq*stride, n1 * sizeof(float));
		    memcpy(pp_weights + ifreq*pp_stride, weights + ifreq*stride, n1 * sizeof(float));
		}

		prepad_buffers[it].finalize_append(n1);

		// Now get pointers to the prepadded data which will be needed for the transform.
		prepad_buffers[it].setup_write(transform_ipos[it], n0, pp_intensity, pp_weights, pp_stride);
	    }

	    // Transform is called here.
	    transforms[it]->process_chunk(intensity, weights, stride, pp_intensity, pp_weights, pp_stride);

	    // Note (n1) here, versus (n1+n2) in call to finalize_write() below.
	    main_buffer.finalize_write(transform_ipos[it], n1);
	    transform_ipos[it] += n1;
	}
	
	curr_ipos = transform_ipos[it];
    }
}


void wi_run_state::end_stream()
{
    if (!this->is_running)
	throw runtime_error("wi_run_state::end_stream() called on non-running stream");

    // We pad the stream with fake weight-zero data, until every transform has "seen"
    // every sample of real data.  First we need to compute the needed amount of padding.
    
    int save_ipos = stream_ipos;
    int target_ipos = stream_ipos;

    for (int it = ntransforms-1; it >= 0; it--) {
	int n1 = round_up(target_ipos - transform_ipos[it], transforms[it]->nt_chunk);
	int n2 = transforms[it]->nt_postpad;

	// target for (it-1)-th transform
	target_ipos = transform_ipos[it] + n1 + n2;
    }

    while (stream_ipos < target_ipos) {
	float *dummy_intensity;
	float *dummy_weights;
	int dummy_stride;
	
	// To add fake weight-zero data to the stream, we call setup_write() with
	// zero_flag=true, followed by finalize_write().
	int nt = min(target_ipos - stream_ipos, nt_stream_maxwrite);
	this->setup_write(nt, dummy_intensity, dummy_weights, dummy_stride, true);
	this->finalize_write(nt);
    }

    // Check on padding calculation
    rf_assert(transform_ipos[ntransforms-1] >= save_ipos);

    for (int it = 0; it < ntransforms; it++) {
	transforms[it]->end_stream();
	prepad_buffers[it].reset();
    }

    this->main_buffer.reset();
    this->is_running = false;
}


}  // namespace rf_pipelines
