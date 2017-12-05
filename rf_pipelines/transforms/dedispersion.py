#!/usr/bin/env python
import rf_pipelines
import pandas as pd
import random
import matplotlib
matplotlib.use('Agg')

simulations = 1000

fdmtConst = 4.14e9
fdmtFreqFactor = (1.0/400)**2 - (1.0/800)**2
dtSample = 1e-3

snrList = []; arrivalTimeList = []; dmList = []
snrListBonsai = []; dmListBonsai = []; timeListBonsai = []
snrListFdmt = []; dmListFdmt = []; timeListFdmt = []
snrListBb = []; dmListBb = []; timeListBb = []


def writeOutputs():

    df = pd.DataFrame({"snr" : snrList, "time" : arrivalTimeList, "dm" : dmList})
    df.to_csv("/scratch/k/kmsmith/ugiri/simulated.csv")

    df = pd.DataFrame({"snr" : snrListBonsai, "time" : timeListBonsai, "dm" : dmListBonsai})
    df.to_csv("/scratch/k/kmsmith/ugiri/bonsai.csv")

    df = pd.DataFrame({"snr" : snrListBb, "time" : timeListBb, "dm" : dmListBb})
    df.to_csv("/scratch/k/kmsmith/ugiri/bb.csv")

    #df = pd.DataFrame({"snr" : snrListFdmt, "time" : timeListFdmt, "dm" : dmListFdmt})
    #df.to_csv("/scratch/k/kmsmith/ugiri/fdmt.csv")


for i in xrange( simulations ):

    transform = []

    SNR =  100#float(random.random()*30)
    DM  = float(random.randrange(0, 1000))

    delayAcrossBand = 4.15e3 * DM * ( ( 1./400.)**2 - (1./800.)**2 )
    arrivalDelay = 15#4.15e3 * DM * ( ( 1./800.)**2 )
    margin = 1.0

    arrivalTime = float( random.randrange(int(0 + arrivalDelay + margin),int( 64 - delayAcrossBand - margin)) )

    snrList.append(SNR); arrivalTimeList.append(arrivalTime); dmList.append(DM)

    stream = rf_pipelines.gaussian_noise_stream(nfreq = 1024, nt_tot = 65536,
                                           nt_chunk = 65536,freq_lo_MHz = 400.0,
                                           freq_hi_MHz = 800.0,     dt_sample = 1.0e-3)

    transform.append( rf_pipelines.frb_injector_transform(snr = SNR,
                                         		  undispersed_arrival_time = arrivalTime,
                                         		  dm = DM,
                                         		  intrinsic_width = 1.0e-3) )
    transform.append( rf_pipelines.bonsai_dedisperser('bonsai_config.txt', track_global_max=True ))

    transform.append(rf_pipelines.bb_dedisperser(dm_start = 0,
                                		 dm_end = 1500,
                                 		 dm_tol = 1.25,    # value recommended in dedisp example code
                                 		 pulse_width_ms = 2.0))

    #transform.append(rf_pipelines.variance_estimator())

    stream.run(transform, "/scratch/k/kmsmith/ugiri/")


    df = pd.read_json("/scratch/k/kmsmith/ugiri/rf_pipeline_0.json")

    snrListBonsai.append(df['transforms'][2]["frb_global_max_trigger"])
    dmListBonsai.append(df['transforms'][2]["frb_global_max_trigger_dm"])
    timeListBonsai.append(df['transforms'][2]["frb_global_max_trigger_tfinal"])

    snrListBb.append(df['transforms'][3]["frb_global_max_trigger"])
    dmListBb.append(df['transforms'][3]["frb_global_max_trigger_dm"])
    timeListBb.append(df['transforms'][3]["frb_global_max_trigger_tfinal"])

    """snrListFdmt.append(SNR)
    with open("/scratch/k/kmsmith/ugiri/fdmt_deltat.out") as f:
        deltat = float( f.readline() )
        dmListFdmt.append(deltat / ( fdmtConst * fdmtFreqFactor * dtSample))
    with open("/scratch/k/kmsmith/ugiri/fdmt_time.out") as f:
        time = float(f.readline())
        timeListFdmt.append(time/1024.0)

    """
    writeOutputs()


