
//
// Canonical setup:
//  Construct stream
//  get params

struct stream_params {
    const int nfreq;
    const double freq_lo_MHz;
    const double freq_hi_MHz;
    const double dt_sample;

    stream_params(int nfreq, double freq_lo_MHz, double freq_hi_Mhz, double dt_sample);
};


// Virtual base class
struct wi_stream {
    virtual stream_params get_params() = 0;
    virtual void run(wi_transform_chain &chain) = 0;

    virtual ~wi_stream() { }
};


// Always used through shared_ptr<>
struct wi_transform {
    stream_params params;
    const int nt_chunk;
    const int nt_prepad;
    const int nt_postpad;

    wi_transform(const stream_params &p, int nt_chunk, int nt_prepad, int nt_postpad);

    // Warning: pp_data, pp_wt have stride nt_prepad!
    // If nt_prepad=0, then they will be NULL pointers.
    virtual void process_chunk(float *data, float *weight, int data_stride, float *pp_data, float *pp_weight) = 0;
};


// Buffers live here
struct wi_transform_chain : noncopyable {
    int ntransforms;
    int nt_headroom;

    std::vector<std::shared_ptr<wi_transform> > transforms;
    std::vector<int> buf_it;
    int stream_it;

    int nt_mbuf;
    float *intensity_mbuf;
    float *weight_mbuf;

    // Called by the stream class
    // Pointers are only valid until next call to stream_advance().
    float *get_stream_intensity_ptr(int nt);
    float *get_stream_weight_ptr(int nt);
    void stream_advance(int nt);

    // Internal helpers
    void _run_transforms();

    // Hmm a bit tricky
    void finalize();
};


void wi_transform_chain::stream_advance(int nt)
{
    if (nt <= 0)
	throw runtime_error("wi_transform_chain::stream_advance(): expected nt > 0");
    if (stream_it + nt > nt_mbuf)
	throw runtime_error("wi_transform_chain::stream_advance(): stream attempted to overwrite buffer");

    this->stream_it += nt;
    this->_run_transforms();
}


void wi_transform_chain::_run_transforms()
{
    int prev_it = stream_it;
    
    for (int i = 0; i < ntransforms; i++) {
	wi_transform *tp = transforms[i].get();

	while (buf_it[i] + tp->nt_chunk + tp->nt_postpad <= prev_it) {
	    tp->process_chunk(xxx);
	    buf_it[i] += tp->nt_chunk;
	}

	prev_it = buf_it[i];
    }
}

    
