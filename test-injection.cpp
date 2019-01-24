#include "rf_pipelines_internals.hpp"

using namespace std;
using namespace rf_pipelines;

int main(int argc, char **argv) {

    const int niter = 50;

    std::random_device rd;
    std::mt19937 rng(rd());
    
    for (int iter = 0; iter < niter; iter++) {
	//if (iter % 50 == 0)
        cout << "test-injection: iteration " << iter << "/" << niter << endl;

	ssize_t nfreq = randint(rng, 1, 1024);
        // To ensure we can use fast_assembled_chunk
        ssize_t nt_per_packet = 16;
        ssize_t nt_chunk = nt_per_packet * randint(rng, 1, 1024);
	ssize_t nt_total = nt_chunk * randint(rng, 1, 10);

        double freq_lo_mhz = 400.;
        double freq_hi_mhz = 800.;
        double dt_sample = 384 * 2.56e-6;

        double sample_rms = 1.0;
        
        auto noise = make_gaussian_noise_stream(nfreq, nt_total,
                                                freq_lo_mhz, freq_hi_mhz,
                                                dt_sample, sample_rms,
                                                nt_chunk);

        vector<string> spool_bufnames;
        spool_bufnames.push_back("INTENSITY");
        auto spool1 = make_shared<pipeline_spool>(spool_bufnames);

        auto inj = make_intensity_injector(nt_chunk);
        
        auto spool2 = make_shared<pipeline_spool>(spool_bufnames);

        auto pipeline = make_shared<rf_pipelines::pipeline> ();
        pipeline->add(noise);
        pipeline->add(spool1);
        pipeline->add(inj);
        pipeline->add(spool2);

        rf_pipelines::run_params rparams;
        rparams.outdir = "";  // disables
        rparams.verbosity = 0;

        // MUST BIND so that the "inj" pipeline stage knows what *nfreq* is.
        pipeline->bind(rparams);

        vector<int> all_offsets;
        vector<int> all_ndata;
        vector<double> all_data;
        
        int ninject = randint(rng, 0, 10);

        for (int j=0; j<ninject; j++) {
            auto injdata = make_shared<inject_data>();
            injdata->mode = 0;
            injdata->sample0 = randint(rng, 0, nt_chunk);

            injdata->sample_offset.resize(nfreq);
            injdata->ndata.resize(nfreq);

            auto uni = std::uniform_real_distribution<>();
            
            int totinj = 0;
            for (int k=0; k<nfreq; k++) {
                int nchunk = randint(rng, 1, 100);
                injdata->sample_offset[k] = randint(rng, 0, nt_total);
                injdata->ndata[k] = nchunk;
                injdata->data.resize(totinj + nchunk);
                for (int m=0; m<nchunk; m++) {
                    injdata->data[totinj + m] = uni(rng);
                }
                totinj += nchunk;
            }
            inj->inject(injdata);

            for (int k=0; k<nfreq; k++)
                all_offsets.push_back(injdata->sample_offset[k] +
                                      injdata->sample0);
            all_ndata.insert(all_ndata.end(),
                             injdata->ndata.begin(),
                             injdata->ndata.end());
            all_data.insert(all_data.end(),
                            injdata->data.begin(),
                            injdata->data.end());

            cout << "Total of " << all_offsets.size() << " chunks, " << all_data.size() << " samples to inject" << endl;
        }

        pipeline->run(rparams);

        string bufname = "INTENSITY";
        auto buf1 = spool1->get_spooled_buffer(bufname);
        auto buf2 = spool2->get_spooled_buffer(bufname);
        
        assert(buf1->cdims.size() == 1);
        assert(buf2->cdims.size() == 1);
        assert(buf1->cdims[0] == nfreq);
        assert(buf2->cdims[0] == nfreq);

        int nt_spool = buf1->nt;
        assert(buf1->nt == buf2->nt);

        int ndiff = 0;
        for (int i=0; i<nt_spool * buf1->csize; i++) {
            if (buf1->data[i] != buf2->data[i])
                ndiff++;
        }
        cout << "Spooled: " << ndiff << " samples differ" << endl;

        ///// add injected data to buf1...
        int data_index = 0;
        for (int i=0; i<all_offsets.size(); i++) {
            // each "inject_data" struct has "nfreq" strings of samples
            int ifreq = i % nfreq;
            int offset = all_offsets[i];
            int nsamples = all_ndata[i];
            for (int j=0; j<nsamples; j++) {
                int itime = offset + j;
                if (itime < 0)
                    continue;
                if (itime >= nt_spool)
                    continue;
                double sample = all_data[data_index + j];
                // FIXME -- downsampling not handled here!!
                buf1->data[ifreq * nt_spool + itime] += sample;
            }
            data_index += nsamples;
        }

        ndiff = 0;
        for (int i=0; i<nt_spool * buf1->csize; i++) {
            if (buf1->data[i] != buf2->data[i])
                ndiff++;
        }
        cout << "Predicted: " << ndiff << " samples differ" << endl;
        assert(ndiff == 0);
        
        // assert equality between buf1 and buf2 (within machine eps)
    }

    cout << "test-injection: pass" << endl;
    return 0;
}

