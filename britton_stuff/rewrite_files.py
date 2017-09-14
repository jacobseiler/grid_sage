#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')

import h5py
import numpy as np
import pylab as plt
from random import sample

import sys
sys.path.append('/home/jseiler/SAGE-stuff/output/')

import AllVars
import PlotScripts
import ReadScripts

snaplow = 8
snaphigh = 10
num_cores = 256
filepath = '/lustre/projects/p004_swin/bsmith/1.6Gpc/means/halo_1721228/dm_gadget/data/'
TypeMax = 6

def write_fof_header(snapshot_idx):

    numpart = np.zeros((TypeMax), dtype = np.int32)
    for core_idx in xrange(0, num_cores):    
        tmp = "groups_%03d/fof_tab_%03d.%d.hdf5" %(snapshot_idx, snapshot_idx, core_idx)
        fname = filepath + tmp

        with h5py.File(fname, 'r') as f:
            try:
                grouplen = f['Group']['GroupLenType']
            except KeyError:
                pass
            else:
                if (np.shape(grouplen)[1] != TypeMax): # Ensure that we're getting the correct number of particle entries.
                    raise ValueError("For file %s GroupLenType had %d entries for the number of particle types when we expect there to be %d." %(fname, np.shape(grouplen)[1], TypeMax))
                 
                local_numpart = np.sum(grouplen[:], axis = 0) # Sums each particle type for this core.
                local_numpart.astype(np.int32) # Recast as explicit integer.
              
                numpart = np.add(numpart,local_numpart) # Add to running total for this snapshot.
        
    ## At this point we have the number of particles (of each type) within this snapshot. ## 

    exit()

def write_snapshot_header(snapshot_idx):

    numpart = np.zeros((TypeMax), dtype = np.int64)
    for core_idx in xrange(0, num_cores):    
        tmp = "snapdir_%03d/snapshot_%03d.%d.hdf5" %(snapshot_idx, snapshot_idx, core_idx)
        fname = filepath + tmp
 
        with h5py.File(fname, 'r') as f:

            for type_idx in xrange(0, TypeMax):
               tmp = "PartType%d" %(type_idx)
               try: 
                   part = f[tmp]['Coordinates']
               except KeyError:
                   pass
               else:
                   numpart[type_idx] += np.shape(part)[0]
                   
    ## At this point we have the number of particles (of each type) within this snapeshot. ##

     
   
   
  
              

if __name__ == '__main__':

    for snapshot_idx in xrange(snaplow, snaphigh + 1):
        #write_fof_header(snapshot_idx) 
        write_snapshot_header(snapshot_idx)
