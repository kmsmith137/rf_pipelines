import rf_pipelines
import pandas as pd
import random
import matplotlib
matplotlib.use('Agg')

injections = 1000
snr_list = []; arrival_time_list = []; dm_list = []

snr_list_bonsai = []; dm_list_bonsai = []; time_list_bonsai = []

snr_list_fdmt = []; dm_list_fdmt = []; time_list_fdmt = []

const = 4.14e9
freq_fac = (1.0/400)**2 - (1.0/800)**2
dt_sample = 1e-3

for i in range(injections):
    snr = float(random.random()*15)
    dm = 15
    s = rf_pipelines.gaussian_noise_stream(nfreq = 1024, nt_tot = 65536, nt_chunk = 65536,freq_lo_MHz = 400.0, freq_hi_MHz = 800.0,     dt_sample = 1.0e-3)

    delay_across_band = 4.15e3 * dm * freq_fac
    arrival_delay = 4.15e3 * dm * ((1./800)**2)
    margin = 1.0
    arrival_time = 15 #float(random.randrange(int(0 + arrival_delay + margin), int( 64 - delay_across_band - margin)))
    snr_list.append(snr)
    arrival_time_list.append(arrival_time)
    dm_list.append(dm)


    t = []
    t.append( rf_pipelines.frb_injector_transform(snr = snr, undispersed_arrival_time = arrival_time, dm = dm, intrinsic_width = 1.0e-3))
    t.append(rf_pipelines.bonsai_dedisperser('bonsai_config.txt', track_global_max=True ))
    t.append(rf_pipelines.variance_estimator())
    s.run(t, "/scratch/k/kmsmith/ugiri")
    df = pd.read_json("/scratch/k/kmsmith/ugiri/rf_pipeline_0.json")

    snr_list_bonsai.append(df['transforms'][2]["frb_global_max_trigger"])
    dm_list_bonsai.append(df['transforms'][2]["frb_global_max_trigger_dm"])
    time_list_bonsai.append(df['transforms'][2]["frb_global_max_trigger_tfinal"])

    snr_list_fdmt.append(snr)
    with open("/scratch/k/kmsmith/ugiri/fdmt_deltat.out") as f:
        deltat = float( f.readline() )

        dm_list_fdmt.append(deltat / ( const * freq_fac * dt_sample))
    with open("/scratch/k/kmsmith/ugiri/fdmt_time.out") as f:
        time = float(f.readline())
        print time
        time_list_fdmt.append(time/1024.0)

    df = pd.DataFrame({"snr" : snr_list_bonsai, "time" : time_list_bonsai, "dm" : dm_list_bonsai})
    df.to_csv("/scratch/k/kmsmith/ugiri/bonsai.csv")

    df = pd.DataFrame({"snr" : snr_list_fdmt, "time" : time_list_fdmt, "dm" : dm_list_fdmt})
    df.to_csv("/scratch/k/kmsmith/ugiri/fdmt.csv")

    df = pd.DataFrame({"snr" : snr_list, "time" : arrival_time_list, "dm" : dm_list})
    df.to_csv("/scratch/k/kmsmith/ugiri/simulated.csv")


