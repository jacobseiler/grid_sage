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

    '''
    numpart_total = np.zeros((TypeMax), dtype = np.int64)
    for core_idx in xrange(0, num_cores):    
        numpart_thisfile = np.zeros((TypeMax), dtype = np.int32)

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
                   numpart_total[type_idx] += np.shape(part)[0] # Running total for number of particles in this snapshot.
                   numpart_thisfile[type_idx] = np.shape(part)[0] # Number of particles for this file.
            try:
                 
            except KeyError:
                pass # If it hasn't been created, move on.
            else:
                del f['Header']['NumPartThisFile'] # If it has, delete all datasets within the header group so we can update them. 
            
            dset = f.create_dataset("Header/NumPartThisFile", (TypeMax,), dtype = np.int32, data = numpart_thisfile) # Write the new header to the file.
    '''     
    ## At this point we have the number of particles (of each type) within this snapeshot. ##
    ## We have also written the Header folder containing the number of particles contained within each file. ##

    # Now open up each file once again to write out the final header fields.

    for core_idx in xrange(0, num_cores):
        
        

if __name__ == '__main__':

    for snapshot_idx in xrange(snaplow, snaphigh + 1):
        #write_fof_header(snapshot_idx) 
        write_snapshot_header(snapshot_idx)
