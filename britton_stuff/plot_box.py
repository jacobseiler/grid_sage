import h5py
import numpy as np


import sys
sys.path.append('/home/jseiler/SAGE-stuff/output/')

import AllVars
import PlotScripts
import ReadScripts

if __name__ == '__main__':

    filepath = '/lustre/projects/p004_swin/bsmith/1.6Gpc/means/halo_1621228/dm_gadget/data/snapdir_000/'
    
    for snapshot_idx in np.arange(0, 11):
        tmp = "snapshot_%03d.129.hdf5" %(snapshot_idx)
        fname = filepath + tmp

        with h5py.File(fname, 'r') as f:

            part_one = f['PartType1']

            print part_one
            exit()
