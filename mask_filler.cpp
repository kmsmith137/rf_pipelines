#include "rf_pipelines_internals.hpp"
#include "bonsai_inlines.hpp"
#include "chlog.hpp"

#include <rf_kernels/online_mask_filler.hpp>

using namespace std;

#define square bonsai::square
#define min3   bonsai::min3
#define max3   bonsai::max3

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode
#endif

mask_filler::mask_filler(const bonsai::config_params &config_) :
    wi_transform("mask_filler"),
    kernel(config_.nfreq)
{
    // This copy is silly, but helps navigate the 'const' minefield!
    bonsai::config_params config(config_);
    config.finalize(0);  // verbosity = 0

    if (config.variance_timescale <= 0.0)
	throw runtime_error("bonsai: in order to use the 'rfi_mask_fill' functionality, you need to initialize 'variance_timescale' in the bonsai config file");

    // Just checking...
    bonsai_assert(config.variance_estimation_timescale > 0.0);
    bonsai_assert(config.reweighting_timescale > 0.0);

    // FIXME -- some of these values should be set during _bind of the pipeline!

    this->nfreq_f = config.nfreq;
    this->nfreq_c = config.variance_nfreq_bins;
    chlog("mask_filler: Setting nt_chunk from config file: " << config.nt_chunk);
    this->nt_chunk = config.nt_chunk;
    this->freq_lo = config.freq_lo_MHz;
    this->freq_hi = config.freq_hi_MHz;
    this->freq_bin_delim = bonsai::make_delim(0, nfreq_f, nfreq_c);

    kernel.v1_chunk = 32;
    kernel.var_weight = (config.dt_sample * kernel.v1_chunk) / config.variance_estimation_timescale;
    kernel.w_clamp = (config.dt_sample * kernel.v1_chunk) / config.reweighting_timescale;
    kernel.w_cutoff = 0.5;   // FIXME add this parameter to config_params or dedisperser::initializer?

    if (kernel.var_weight >= 0.2)
	throw runtime_error("bonsai: 'variance_timescale' parameter is too small");

    // Determine ring buffer capacity, based on max dedispersion delay.
    rb_capacity = 0;
    for (int itree = 0; itree < config.ntrees; itree++) {
	int n0 = config.nt_ftree_pad[itree] + config.tree_size[itree];        // max delay in "tree samples"
	int n1 = (n0 * config.nds[itree]) / (config.nups[itree] * nt_chunk);  // max delay in "input chunks"
	rb_capacity = max(rb_capacity, n1+4);               // the +4 accounts for various boundary effects
    }

    this->frac_delay.resize(nfreq_c + 1, 0.0);

    // Initialize frac_delay: a length-(nfreq_c+1) array containing dispersion delays at frequency channel
    // boundaries, as fraction of total dispersion delay across the band.  Thus frac_delay[0] = 1 and
    // frac_delay[nfreq_c] = 0.  (Note that "delay" means t(nu_lo)-t(nu), and channels are ordered from
    // highest frequency to lowest.)
    
    for (int ifreq_c = 0; ifreq_c <= nfreq_c; ifreq_c++) {
	double freq = (ifreq_c * freq_lo + (nfreq_c-ifreq_c) * freq_hi) / nfreq_c;
	double t = 1.0/square(freq_lo) - 1.0/square(freq_hi);
    
	frac_delay[ifreq_c] = (1.0/square(freq_lo) - 1.0/square(freq)) / t;
    }
}

void mask_filler::_bind_transform_rb(ring_buffer_dict &rb_dict) {
    int nds = nt_chunk;
    vector<ssize_t> cdims;
    cdims.push_back(nfreq_f);
    this->rb_variance = create_buffer(rb_dict, "RUNNING_VARIANCE", cdims, nds);
    this->rb_weight   = create_buffer(rb_dict, "RUNNING_WEIGHT", cdims, nds);
    cdims.clear();
    cdims.push_back(nfreq_c);
    this->rb_wvar     = create_buffer(rb_dict, "RUNNING_WVAR", cdims, nds);
    // Set ring buffer nt_maxlag
    int nhist = this->rb_capacity * nt_chunk;
    this->rb_variance->update_params(1, nhist);
    this->rb_weight->update_params(1, nhist);
    this->rb_wvar->update_params(1, nhist);
    this->rb_variance->dense = true;
    this->rb_weight->dense = true;
    this->rb_wvar->dense = true;
}

void mask_filler::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos) {
    // The kernel does the following:
    //   - Updates its internal estimates of the running variance.
    //   - Updates the running_weights which it applies.

    //chlog("mask_filler::process_chunk mask_fill_in_place.  nt_chunk " << nt_chunk << ", nt_chunk_out " << nt_chunk_out);
    //chlog("rb_weight: " << rb_weight->get_info());
    //chlog("rb_variance: " << rb_variance->get_info());
    //chlog("rb_wvar: " << rb_wvar->get_info());
    
    kernel.mask_fill_in_place(nt_chunk, intensity, istride, weights, wstride);

    //   - Fills the output 
    //   - Multiplies the intensity by the running_weights.
    //kernel.mask_fill_and_multiply(nt_chunk, out, ostride, intensity, istride, weights, wstride);
    
    // These point to the kernel's running variance estimate and applied weights.
    const float *running_weights = kernel.running_weights.get();
    const float *running_var = kernel.running_var.get();

    ring_buffer_subarray subw(rb_weight,   pos, pos+nt_chunk, ring_buffer::ACCESS_APPEND);
    ring_buffer_subarray subv(rb_variance, pos, pos+nt_chunk, ring_buffer::ACCESS_APPEND);
    ring_buffer_subarray subwv(rb_wvar,    pos, pos+nt_chunk, ring_buffer::ACCESS_APPEND);

    for (ssize_t i = 0; i < rb_weight->csize; i++) {
        subw.data[i*subw.stride] = running_weights[i];
        subv.data[i*subv.stride] = running_var[i];
     }

    // We ring-buffer the (running_weights)^2 * (running_variance).
    // This is the variance of the data after the running_weights have been applied.
    
    // Outer loop over coarse frequencies.
    for (int ifreq_c = 0; ifreq_c < nfreq_c; ifreq_c++) {
	// 1d array of length rb_capacity (with periodic ring-buffer index)

	// Inner loop over fine frequencies, accumulating (weight^2 * variance).
	float var = 0.0;
	for (int ifreq_f = freq_bin_delim[ifreq_c]; ifreq_f < freq_bin_delim[ifreq_c+1]; ifreq_f++)
	    var += square(running_weights[ifreq_f]) * running_var[ifreq_f];

	var /= (freq_bin_delim[ifreq_c+1] - freq_bin_delim[ifreq_c]);

        subwv.data[ifreq_c * subwv.stride] = var;
    }

    this->nchunks_processed++;
    chlog("mask_filler: process_chunk done!");
}


// Note: returns 0 if no variance estimate is available.
//
// FIXME it would be slightly better to have eval_weighted_variance() return a min/max
// variance, based on rf_kernels::mask_filler::chunk_{min,max}_weight.  Is this
// worth implementing?

float mask_filler::eval_weighted_variance(int ifreq_c, double ns) const
{
    bonsai_assert(ifreq_c >= 0);
    bonsai_assert(ifreq_c < nfreq_c);
    
    //const float *var_1d = &rb_variance[ifreq_c * rb_capacity];

    // Convert timestamp from samples to chunks, and range-check timestamp.
    double nc = ns / nt_chunk;
    bonsai_assert(nc > nchunks_processed - rb_capacity + 2.001);
    bonsai_assert(nc < nchunks_processed + 0.001);

    if (ns <= 0.0)
	return 0.0;
    
    // Decided to linearly interpolate, but I doubt this detail affects anything!
    // Interpolation will involve indices ic-1,ic (not ic,ic+1).
    // This is because ring buffer index 'ic' stores the variance estimate
    // at the end (not beginning) of the ic-th chunk.
    
    int ic = min(int(nc), nchunks_processed-1);

    double v0,v1;
    if (ic > 0) {
        ring_buffer_subarray subwv(rb_wvar, (ic-1)*nt_chunk, (ic+1)*nt_chunk, ring_buffer::ACCESS_READ);
        v0 = subwv.data[ifreq_c * subwv.stride];
        v1 = subwv.data[(ifreq_c+1) * subwv.stride];
    } else {
        ring_buffer_subarray subwv(rb_wvar, ic*nt_chunk, (ic+1)*nt_chunk, ring_buffer::ACCESS_READ);
        v0 = 0.0;
        v1 = subwv.data[ifreq_c * subwv.stride];
    }
    //double v0 = (ic > 0) ? var_1d[(ic-1) % rb_capacity] : 0.0;
    //double v1 = var_1d[(ic) % rb_capacity];
    double xc = nc - ic;
    
    return (1.0-xc)*v0 + (xc)*v1;
}


// The output array 'vout' has shape (nfreq_c, 2).
void mask_filler::eval_sweep_variance(double ns_d, double ns_f, float *vout, bool update) const
{
    for (int ifreq_c = 0; ifreq_c < nfreq_c; ifreq_c++) {
	// (ns0, ns1) = (lower, upper) limit of pulse sweep in channel.
	double ns0 = ns_f - ns_d * frac_delay[ifreq_c];
	double ns1 = ns_f - ns_d * frac_delay[ifreq_c+1];

	// Reminder: Either of these calls can return 0.0 if variance estimate is unavailable.	
	float var0 = this->eval_weighted_variance(ifreq_c, ns0);
	float var1 = this->eval_weighted_variance(ifreq_c, ns1);

	if (update) {
	    vout[2*ifreq_c] = min3(vout[2*ifreq_c], var0, var1);
	    vout[2*ifreq_c+1] = max3(vout[2*ifreq_c+1], var0, var1);
	}
	else {
	    vout[2*ifreq_c] = min(var0, var1);
	    vout[2*ifreq_c+1] = max(var0, var1);
	}
    }
}

void mask_filler::get_weights_and_variances(vector<float>* weights,
                                                   vector<float>* variances) const
{
    if (nchunks_processed == 0)
      return;
    
    if (weights) {
        const float *running_weights = kernel.running_weights.get();
        weights->insert(weights->end(), running_weights, running_weights + nfreq_f);
    }
    if (variances) {
        const float *running_var = kernel.running_var.get();
        variances->insert(variances->end(), running_var, running_var + nfreq_f);
    }
}

//shared_ptr<wi_transform> make_bonsai_dedisperser(const std::string &config_filename, const bonsai_initializer &ini_params) {}


} // namespace
