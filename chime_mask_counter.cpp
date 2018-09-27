#include "rf_pipelines_internals.hpp"
#include "rf_pipelines_inventory.hpp"
#include "ch_frb_io.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}; // pacify emacs c-mode
#endif

chime_mask_counter::chime_mask_counter(string where_) :
    // First argument is 'nt_chunk'.
    mask_counter_transform(ch_frb_io::constants::nt_per_assembled_chunk, where_, "chime_mask_counter")
{}

void chime_mask_counter::set_stream(std::shared_ptr<ch_frb_io::intensity_network_stream> _stream, int _beam)
{
    if (stream)
	throw runtime_error("double call to chime_mask_counter::set_stream()");
    if (this->state != UNBOUND)
	throw runtime_error("chime_mask_counter::set_stream() must be called before bind()");
    
    // Note: we don't assert(stream->nrfifreq == this->nfreq) here, because set_stream()
    // is called before chime_mask_counter "knows" the value of nfreq.
    stream = _stream;
    beam = _beam;
}

void chime_mask_counter::_bind_transform(Json::Value &json_attrs)
{
    // _bind_transform() is the earliest place we can put this assert.
    if (stream && (stream->ini_params.nrfifreq != this->nfreq))
	throw runtime_error("chime_mask_counter: value of 'nrfifreq' does not match the frequency resolution in the RFI transform chain");

    // Should be redundant with asserts elsewhere in rf_pipelines, but just being paranoid!
    rf_assert(this->nds == 1);
}

void chime_mask_counter::_start_pipeline(Json::Value &j)
{
    this->initial_fpga_count = uint64_t_from_json(j, "initial_fpga_count");
    this->fpga_counts_per_sample = int_from_json(j, "fpga_counts_per_sample");
    this->fpga_counts_initialized = true;
}

void chime_mask_counter::_process_chunk(float *intensity, ssize_t istride, float *weights, ssize_t wstride, ssize_t pos)
{
    if (!stream) {
        cout << "chime_mask_counter: processing chunk, but stream not set" << endl;
        mask_counter_transform::_process_chunk(intensity, istride, weights, wstride, pos);
        return;
    }

    // Reminder: previous asserts have already checked that
    //   this->nfreq == stream->nrfifreq
    //   this->nt_chunk == ch_frb_io::nt_per_assembled_chunk
    //   this->nds == 1
    
    if (!fpga_counts_initialized)
	throw runtime_error("rf_pipelines::chime_mask_counter internal error: fpga count fields were not initialized as expected");
    
    cout << "chime_mask_counter: finding chunk for pos " << pos << endl;

    // The 'pos' argument is the current pipeline position in units of time samples (not FPGA counts)
    uint64_t fpga_counts = pos * this->fpga_counts_per_sample + this->initial_fpga_count;
    //cout << "FPGA counts: " << fpga_counts << endl;

    // The last argument in find_assembled_chunk() is 'toplevel'.
    shared_ptr<ch_frb_io::assembled_chunk> chunk = stream->find_assembled_chunk(beam, fpga_counts, true);

    // These should all be redundant with asserts in ch_frb_io, but a little paranoia never hurts.
    if (!chunk)
	throw runtime_error("chime_mask_counter: find_assembled_chunk() returned empty pointer");
    if (!chunk->rfi_mask)
	throw runtime_error("chime_mask_counter: find_assembled_chunk() returned chunk with no RFI mask");
    if (chunk->nrfifreq != nfreq)
	throw runtime_error("chime_mask_counter: find_assembled_chunk() returned chunk with mismatched 'nrfifreq'");
    if (chunk->has_rfi_mask)
	throw runtime_error("chime_mask_counter: find_assembled_chunk() returned chunk with has_rfi_mask=true");
    if (chunk->binning != 1)
	throw runtime_error("chime_mask_counter: find_assembled_chunk() returned chunk with binning != 1");

    mask_counter_measurements meas;
    meas.pos = pos;
    meas.nsamples = nfreq*nt_chunk;
    meas.nsamples_masked = 0;
    meas.nt = nt_chunk;
    meas.nt_masked = 0;
    meas.nf = nfreq;
    meas.nf_masked = 0;
    meas.freqs_masked = shared_ptr<uint16_t>((uint16_t*)calloc(nfreq, sizeof(uint16_t)), free);
    meas.times_masked = shared_ptr<uint16_t>((uint16_t*)calloc(nt_chunk, sizeof(uint16_t)), free);

    uint16_t* fm = meas.freqs_masked.get();
    uint16_t* tm = meas.times_masked.get();

    uint8_t* rfimask = chunk->rfi_mask;

    for (int i_f=0; i_f<nfreq; i_f++) {
        for (int i_t=0; i_t<nt_chunk/8; i_t++) {
            uint8_t m_out = 0;
            for (int j=0; j<8; j++) {
                if (weights[i_f*wstride + 8*i_t + j] == 0) {
                    meas.nsamples_masked++;
                    fm[i_f]++;
                    tm[8*i_t+j]++;
                } else
                    m_out |= (1 << j);
            }
            rfimask[i_f*nt_chunk/8 + i_t] = m_out;
        }
    }

    chunk->has_rfi_mask = true;

    for (int i_f=0; i_f<nfreq; i_f++)
        if (fm[i_f] == nt_chunk)
            meas.nf_masked++;
    for (int i_t=0; i_t<nt_chunk; i_t++)
        if (tm[i_t] == nfreq)
            meas.nt_masked++;
            
    cout << "chime_mask_counter " << where << ", pos " << pos << ": N samples masked: " << meas.nsamples_masked << "/" << (meas.nsamples) << "; n times " << meas.nt_masked << "/" << meas.nt << "; n freqs " << meas.nf_masked << "/" << meas.nf << endl;
    for (const auto &cb : callbacks)
        cb->mask_count(meas);

    // Notify stream's output_devices that a chunk has had its
    // rfi_mask filled in.
    for (auto od : stream->ini_params.output_devices)
        od->filled_rfi_mask(chunk);
}


Json::Value chime_mask_counter::jsonize() const
{
    Json::Value ret;
    
    ret["class_name"] = "chime_mask_counter";
    ret["where"] = where;
    return ret;
}


shared_ptr<chime_mask_counter>
chime_mask_counter::from_json(const Json::Value &j)
{
    string where = string_from_json(j, "where");
    return make_shared<chime_mask_counter> (where);
}

namespace {
    struct _init {
        _init() {
            pipeline_object::register_json_deserializer("chime_mask_counter", chime_mask_counter::from_json);
        }
    } init;
}

// Externally callable
shared_ptr<wi_transform> make_chime_mask_counter(string where)
{
    return make_shared<chime_mask_counter> (where);
}


}  // namespace rf_pipelines



