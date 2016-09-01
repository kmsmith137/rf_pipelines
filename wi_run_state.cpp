// Note: I haven't systematically documented the C++ interface to rf_pipelines,
// so the level of documentation will be hit-or-miss.  Also please note that the
// python-wrapping in rf_pipelines_c.cpp is kind of a mess which I hope to improve
// soon.  In the meantime if you want to python-wrap a C++ class, just email me
// and I'll help navigate the mess!

#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif


wi_run_state::wi_run_state(const wi_stream &stream, const vector<shared_ptr<wi_transform> > &transforms_, const shared_ptr<outdir_manager> &manager_, bool noisy_) :
    nfreq(stream.nfreq),
    nt_stream_maxwrite(stream.nt_maxwrite),
    manager(manager_),
    ntransforms(transforms_.size()),
    transforms(transforms_),
    dt_sample(stream.dt_sample),
    substream_start_time(0.0),
    stream_curr_time(0.0),
    transform_ipos(transforms_.size(), 0),
    stream_ipos(0),
    state(0),
    isubstream(0),
    nt_pending(0),
    noisy(noisy_),
    prepad_buffers(transforms_.size())
{
    if (!nfreq)
	throw runtime_error("wi_run_state constructor called on uninitialized stream");
    if (ntransforms < 1)
	throw runtime_error("wi_run_state constructor called on empty transform list");
    if (!manager)
	throw runtime_error("wi_run_state constructor called with empty manager pointer");
}


void wi_run_state::start_substream(double t0)
{
    if ((this->state > 0) && (this->state < 4))
	throw runtime_error("rf_transforms: logic error in stream: double call to start_substream() (maybe a call to end_substream() is missing somewhere?)");

    this->substream_start_time = t0;
    this->stream_curr_time = t0;

    for (int it = 0; it < ntransforms; it++)
	transforms[it]->start_substream(this->isubstream, t0);

    // Allocate main buffer

    ssize_t nt_contig = nt_stream_maxwrite;
    for (ssize_t it = 0; it < ntransforms; it++)
	nt_contig = max(nt_contig, transforms[it]->nt_chunk + transforms[it]->nt_postpad);

    ssize_t nt_ring = transforms[0]->nt_chunk + transforms[0]->nt_postpad + nt_stream_maxwrite;

    for (int it = 1; it < ntransforms; it++) {
	ssize_t g = transforms[it-1]->nt_chunk;
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
	ssize_t n0 = transforms[it]->nt_prepad;
	ssize_t n1 = transforms[it]->nt_chunk;

	if (!n0)
	    continue;

	ssize_t nt_contig = max(n0, n1);
	ssize_t nt_ring = n0 + n1;
	this->prepad_buffers[it].construct(nfreq, nt_contig, nt_ring);

	// shift prepad buffer by n0 as noted above
	this->prepad_buffers[it].append_zeros(n0);
    }

    this->clear_per_substream_data();
    this->state = 1;
}


void wi_run_state::setup_write(ssize_t nt, float* &intensityp, float* &weightp, ssize_t &stride, bool zero_flag, double t0)
{
    // states 1,3 are OK
    if ((this->state == 0) || (this->state == 4))
	throw runtime_error("rf_transforms: logic error in stream: setup_write() was called without prior call to start_substream()");
    if (this->state == 2)
	throw runtime_error("rf_transforms: logic error in stream: double call to setup_write(), without call to finalize_write() in between");

    if (nt <= 0)
	throw runtime_error("rf_transforms: logic error in stream: setup_write() was called with nt <= 0");
    if (nt > this->nt_stream_maxwrite)
	throw runtime_error("rf_transforms: logic error in stream: setup_write() was called with nt > nt_maxwrite");

    if (fabs(t0 - stream_curr_time) >= 1.0e-2 * dt_sample)
	throw runtime_error("rf_transforms: timestamp jitter is not allowed to exceed 1% of the sample length");

    this->main_buffer.setup_append(nt, intensityp, weightp, stride, zero_flag);
    this->stream_curr_time = t0;
    this->state = 2;
    this->nt_pending = nt;
}


void wi_run_state::setup_write(ssize_t nt, float* &intensityp, float* &weightp, ssize_t &stride, bool zero_flag)
{
    double t0 = stream_curr_time;
    this->setup_write(nt, intensityp, weightp, stride, zero_flag, t0);
}


void wi_run_state::finalize_write(ssize_t nt)
{
    if (this->state == 3)
	throw runtime_error("rf_transforms: logic error in stream: double call to finalize_write()");
    if (this->state != 2)
	throw runtime_error("rf_transforms: logic error in stream: finalize_write() was called without prior call to setup_write()");
    if (nt != this->nt_pending)
	throw runtime_error("rf_transforms: logic error in stream: values of 'nt' in setup_write() and finalize_write() don't match");

    // stream_ipos and stream_curr_time get updated at the end
    this->main_buffer.finalize_append(nt);

    ssize_t curr_ipos = this->stream_ipos + nt;

    for (int it = 0; it < ntransforms; it++) {
	ssize_t n0 = transforms[it]->nt_prepad;
	ssize_t n1 = transforms[it]->nt_chunk;
	ssize_t n2 = transforms[it]->nt_postpad;

	while (transform_ipos[it] + n1 + n2 <= curr_ipos) {

	    // Logic for running transform is here

	    float *intensity = nullptr;
	    float *weights = nullptr;
	    ssize_t stride = 0;

	    float *pp_intensity = nullptr;
	    float *pp_weights = nullptr;
	    ssize_t pp_stride = 0;

	    // Note (n1+n2) here, versus (n1) in call to finalize_write() below.
	    main_buffer.setup_write(transform_ipos[it], n1+n2, intensity, weights, stride);
	    
	    if (n0 > 0) {
		// Transform is prepadded.  We need to get pointers to the prepadded data, and
		// we also need to save data which will be overwritten by running the transform,
		// since it will be prepadded data in a future iteration.

		bool zero_flag = false;
		prepad_buffers[it].setup_append(n1, pp_intensity, pp_weights, pp_stride, zero_flag);

		for (ssize_t ifreq = 0; ifreq < nfreq; ifreq++) {
		    memcpy(pp_intensity + ifreq*pp_stride, intensity + ifreq*stride, n1 * sizeof(float));
		    memcpy(pp_weights + ifreq*pp_stride, weights + ifreq*stride, n1 * sizeof(float));
		}

		prepad_buffers[it].finalize_append(n1);

		// Now get pointers to the prepadded data which will be needed for the transform.
		prepad_buffers[it].setup_write(transform_ipos[it], n0, pp_intensity, pp_weights, pp_stride);
	    }

	    // Transform is called here.
	    double t0 = this->stream_curr_time + dt_sample * (transform_ipos[it] - stream_ipos);
	    double t1 = this->stream_curr_time + dt_sample * (transform_ipos[it] - stream_ipos + n1);

	    struct timeval tv0 = get_time();
	    transforms[it]->process_chunk(t0, t1, intensity, weights, stride, pp_intensity, pp_weights, pp_stride);
	    transforms[it]->time_spent_in_transform += time_diff(tv0, get_time());

	    // Note (n1) here, versus (n1+n2) in call to finalize_write() below.
	    main_buffer.finalize_write(transform_ipos[it], n1);
	    transform_ipos[it] += n1;
	}
	
	curr_ipos = transform_ipos[it];
    }

    this->stream_curr_time += dt_sample * nt;
    this->stream_ipos += nt;
    this->state = 3;
    this->nt_pending = 0;
}


void wi_run_state::end_substream()
{
    if (this->state == 4)
	throw runtime_error("rf_transforms: logic error in stream: double call to end_substream()");
    if (this->state != 3)
	throw runtime_error("rf_transforms: logic error in stream: call to end_substream() without prior call to start_substream()");

    // We pad the stream with fake weight-zero data, until every transform has "seen"
    // every sample of real data.  First we need to compute the needed amount of padding.
    
    ssize_t save_ipos = stream_ipos;
    ssize_t target_ipos = stream_ipos;

    for (int it = ntransforms-1; it >= 0; it--) {
	ssize_t n1 = round_up(target_ipos - transform_ipos[it], transforms[it]->nt_chunk);
	ssize_t n2 = transforms[it]->nt_postpad;

	// target for (it-1)-th transform
	target_ipos = transform_ipos[it] + n1 + n2;
    }

    while (stream_ipos < target_ipos) {
	float *dummy_intensity;
	float *dummy_weights;
	ssize_t dummy_stride;
	
	// To add fake weight-zero data to the stream, we call setup_write() with
	// zero_flag=true, followed by finalize_write().
	ssize_t nt = min(target_ipos - stream_ipos, nt_stream_maxwrite);
	this->setup_write(nt, dummy_intensity, dummy_weights, dummy_stride, true);
	this->finalize_write(nt);
    }

    // Check on padding calculation
    rf_assert(transform_ipos[ntransforms-1] >= save_ipos);

    for (int it = 0; it < ntransforms; it++)
	transforms[it]->end_substream();

    // Write outputs
    if (noisy) {
	cerr << ("rf_pipelines: processed " + to_string(save_ipos) + " samples\n");
	for (int it = 0; it < ntransforms; it++)
	    cerr << "    Transform " << it << ": " << transforms[it]->time_spent_in_transform << " sec  [" << transforms[it]->get_name() << "]\n";
    }

    this->write_per_substream_json_file();
    this->clear_per_substream_data();

    // Deallocate buffers and advance state
    this->main_buffer.reset();
    for (int it = 0; it < ntransforms; it++)
	this->prepad_buffers[it].reset();

    this->state = 4;
    this->isubstream++;
}


void wi_run_state::write_per_substream_json_file()
{
    Json::Value json_all;
    json_all["nsamples"] = int64_t(stream_ipos);
    // more things will go here!

    for (const shared_ptr<wi_transform> &t: transforms) {
	string transform_name = "<error in wi_transform::get_name()>";

	try {
	    transform_name = t->get_name();
	} catch (...) {
	    cerr << "warning: wi_transform::get_name() threw exception";
	}

	Json::Value &json_t = t->json_misc;   // includes everything except "time" and "plots"
	json_t["name"] = transform_name;
	json_t["time"] = t->time_spent_in_transform;

	for (const shared_ptr<plot_group> &g: t->plot_groups) {
	    if (g->is_empty)
		continue;

	    Json::Value json_g;
	    json_g["name"] = g->name;
	    json_g["nt_per_pix"] = g->nt_per_pix;
	    json_g["ny"] = g->ny;
	    json_g["it0"] = g->curr_it0;
	    json_g["it1"] = g->curr_it1;
	    json_g["files"].append(g->files);

	    json_t["plots"].append(json_g);
	}

	json_all["transforms"].append(json_t);
    }

    manager->write_per_substream_json_file(isubstream, json_all, noisy);
}


void wi_run_state::clear_per_substream_data()
{
    for (const shared_ptr<wi_transform> &t: transforms) {
	t->time_spent_in_transform = 0.0;
	t->json_misc.clear();

	for (const shared_ptr<plot_group> &g: t->plot_groups) {
	    g->is_empty = true;
	    g->files.clear();
	}
    }
}


}  // namespace rf_pipelines
