#include <string>
#include <iostream>

#include "ch_frb_io.hpp"
#include "assembled_chunk_msgpack.hpp"

#include "rf_pipelines.hpp"

using namespace std;
using namespace ch_frb_io;
using namespace rf_pipelines;

int main() {
    string fn = "00000042.h5";
    intensity_hdf5_file fin(fn, true);

    float* intensity = reinterpret_cast<float*>(malloc(sizeof(float) * fin.nfreq * fin.nt_logical));
    float* weight    = reinterpret_cast<float*>(malloc(sizeof(float) * fin.nfreq * fin.nt_logical));
    int t0, nt;
    t0 = 0;
    nt = fin.nt_logical;

    fin.get_unpolarized_intensity(intensity, weight, t0, nt);

    // There's a slight mismatch here: 1017 time samples in the h5 files,
    // vs 1024 demanded by the assembled_chunk.

    cout << "dt_sample: " << fin.dt_sample << endl;
    
    int beam = 1;
    int nupfreq = (fin.nfreq / ch_frb_io::constants::nfreq_coarse_tot);
    cout << "Nupfreq: " << nupfreq << endl;

    int nt_per_packet = 16; // ??

    int fpga_counts_per_sample = fin.dt_sample / 2.56e-6;
    cout << "FPGA counts per sample: " << fpga_counts_per_sample << endl;

    int ichunk = fin.time_lo * fin.dt_sample / (fpga_counts_per_sample * ch_frb_io::constants::nt_per_assembled_chunk);
    cout << "ichunk " << ichunk << endl;

    auto chunk = assembled_chunk::make(beam, nupfreq, nt_per_packet,
                                       fpga_counts_per_sample, ichunk);

    vector<string> fns;
    fns.push_back("00000040.h5");
    fns.push_back("00000041.h5");
    fns.push_back("00000042.h5");
    fns.push_back("00000043.h5");
    fns.push_back("00000044.h5");
    fns.push_back("00000045.h5");
    fns.push_back("00000046.h5");
    fns.push_back("00000047.h5");
    fns.push_back("00000048.h5");
    fns.push_back("00000049.h5");
    auto stream = make_chime_stream_from_filename_list(fns);

    string dest = "127.0.0.1:10252";

    // For full CHIME, we expect to use (nbeam, nfreq_coarse, nupfreq, ntsamp) = (8, 4, 16, 16)

    int nfreq_coarse_per_packet = 4;
    int nt_per_chunk = ch_frb_io::constants::nt_per_assembled_chunk;
    //nt_per_packet = 16;
    //nt_per_packet = 4;
    nt_per_packet = 2;
    float wt_cutoff = 1e6;
    float target_gbps = 0.;
    
    auto packetizer = make_chime_packetizer(dest, nfreq_coarse_per_packet,
                                            nt_per_chunk, nt_per_packet,
                                            wt_cutoff, target_gbps);


    vector<shared_ptr<wi_transform> > transforms;
    transforms.push_back(packetizer);
    
    stream->run(transforms);
    
}
