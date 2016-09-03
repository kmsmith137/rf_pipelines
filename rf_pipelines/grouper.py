#!/usr/bin/env python


import sys
import numpy as np
import h5py

event_dtype = [('dm',float),('time',float),('snr',float),
               ('beta',float),('sm',float),('itree',int),
               ('curve','11float')]

def group_bonsai_output(thr=8., fn='triggers.hdf5'):
    """
    This function is used to group triggers from bonsai to form L1 events.
    It is designed to mimic the eventual chunked form of the bonsai output.


    Usage:
        
        events = rf_pipelines.group_bonsai_output(thr=8.,fn='triggers.hdf5')

        "thr=8." is the threshold above which an event may be created

        "fn=triggers.hdf5'" is the filename of the bonsai output


    Output:

        recarray with dtype = [('dm',float),('time',float),('snr',float),
                               ('beta',float),('sm',float),('itree',int),
                               ('curve','11float')]


    Warnings:

        - To emulate chunking, a "bonsai_config.txt" file is expected to exist

        - Right now, only a single threshold is applied and no sifting is done.
        
        - Event datatype is subject to change.  (e.g. 'time' -> np.datetime64)
        
        - dm, time, beta, sm are not actually calculated.  Instead, indices from
          the trigger arrays are used.
        
        - 'curve' size is subject to change (would like to do some more tests)
        
        - Right now, the number of coarse-grained dm intervals for each tree is
          assumed to be the same. 
        
        - For time intervals per chunk for each tree, two options are available
          depending on how nt_per_trigger is set in bonsai_config.txt
            1. single value (e.g. 64)
            2. list that gives equal length triggers in time (e.g. [128, 64, 32])
    """
    try:
        h5 = h5py.File(fn,'r')
    except IOError:
        print "\nERROR! '"+fn+"' doesn't exist."
        sys.exit(1)
    T  = [h5[key]['TRIGGERS'][:].T for key in h5.keys() if 'TREE' in key]
    gulps  = get_gulps()
    stream = multitree_stream(T,gulps)
    events = []
    event_buffer = []
    for itree, t0, chunk in stream:
        above = np.argwhere(chunk[1:-1] > thr)
        for t,beta,sm,dm in above:
            t += 1
            nhood = chunk[t-1:t+2,beta,sm,dm-5:dm+6]
            if nhood.shape != (3,11) or np.argmax(nhood) != nhood.size/2: continue
            time = t0+t*2**itree if gulps[0] == gulps[1] else t0+t
            dm *= 2**itree
            e = np.array((dm, time, nhood.max(), beta, sm,
                          itree, nhood.max(0)),event_dtype).view(np.recarray)
            insert_event(e,event_buffer)
        if itree == len(T)-1:
            events += event_buffer
            event_buffer = []
    return np.array(events).view(np.recarray)


def insert_event(e1, event_buffer):
    """ 
    Adds event to buffer if it is new or superior in S/N to a
    previous event in the buffer (different downsampling)
    """
    for i, e2 in enumerate(event_buffer): 
        if abs(e1.time-e2.time) < 2 and abs(e1.dm-e2.dm) < 5:
            if e1.snr > e2.snr: event_buffer[i] = e1
            return
    event_buffer.append(e1)


def get_gulps():
    """
    Reads 'bonsai_config.txt' to determine the shape of the bonsai output
    What gulping is used will depend on nt_tree and nt_per_trigger.
    """
    try:
        with open('bonsai_config.txt') as f:
            for l in f: exec l
        ntrees = len(nds) if type(nds) is list else 1
        if type(nt_per_trigger) is list:
            return [nt_tree/x for x in nt_per_trigger]
        else:
            return [nt_tree/nt_per_trigger for _ in xrange(ntrees)]
    except:
        print "\nERROR! Couldn't determine the chunk sizes of bonsai output."
        print "Make sure 'bonsai_config.txt' is in the working directory,",
        print "and nt_tree, nt_per_trigger, and nds are all defined within."
        sys.exit(1)


def multitree_stream(trees,gulps):
    """ 
    Generator that uses tree_picker and tree_stream
    to yield buffered chunks from all trees in the right order
    """
    tree_ids = tree_picker(len(trees))
    streams  = [tree_stream(t,g) for t,g in zip(trees,gulps)]
    while True:
        itree = tree_ids.next()
        t0, chunk = streams[itree].next()
        yield itree, t0*2**itree, chunk


def tree_picker(ntrees):
    """ 
    Generator that tells us which tree just completed.
    Example sequence: 0, 0, 1, 0, 0, 1, 2...
    """
    i = 0
    while True:
        i += 1
        for itree in xrange(ntrees):
            if i % 2**itree == 0: yield itree 


def tree_stream(tree,gulp,buff=1):
    """ 
    Generator that yields buffered chunks of 
    coarse-grained triggers for a single tree
    """
    b = np.zeros((buff*2+gulp,)+tree.shape[1:])
    i = 0
    while i+gulp < tree.shape[0]:
        b[:2*buff] = b[-2*buff:]
        b[2*buff:] = tree[i:i+gulp]
        yield i, b
        i += gulp
