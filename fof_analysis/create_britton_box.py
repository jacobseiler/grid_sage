#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import argparse
import sys
import os
import h5py
from tqdm import tqdm

def copy_group(file_in, file_out, key):
    """
    Copies HDF5 group into a new HDF5 file (with same data-structure).

    Parameters
    ----------

    file_in, file_out: Open HDF5 files.  Required.
        HDF5 files for the data being copied (file_in) and the file the
        data is being copied to (file_out).

    key: String.  Required.
        Name of the HDF5 group being copied.

    Returns
    ----------

    None.
    """

    group_path = file_in[key].parent.name  # Name of the group path.
    group_id = file_out.require_group(group_path)  # Create the group.
    name = "{0}".format(key)  # Name the group.
    file_in.copy(name, group_id, name=key)  # Copy over the data.

def scale_positions(snap, chunk, snapdir, outdir):

    if snap < 66: # Slightly different naming scheme depending upon the snapshot. 
        fname = "{0}/snapdir_{1:03d}/snapdir_{1:03d}/snapshot_{1:03d}.{2}.hdf5".format(snapdir,
                                                                                       snap,
                                                                                       chunk)
    else:
        fname = "{0}/snapdir_{1:03d}/snapshot_{1:03d}.{2}.hdf5".format(snapdir, 
                                                                       snap,
                                                                       chunk)
                                         

    fout_fname = "{0}/snapshot_{1:03d}.{2}.hdf5".format(outdir, snap, chunk) 

    with h5py.File(fname, "r") as f_in,\
         h5py.File(fout_fname, "w") as f_out:

        tmp = "PartType1"
        try:
            particles = f_in[tmp]
        except KeyError:
            numpart = 0
            min_value = 999
            max_value = -999 
            
            pass
        else:
            copy_group(f_in, f_out, tmp)
            copy_group(f_in, f_out, "Header")

            coordinates = particles['Coordinates'][:]
            min_value = np.min(coordinates)
            max_value = np.max(coordinates)

            f_out[tmp]['Coordinates'][:] = coordinates - min_value 

            numpart = len(coordinates)

    return numpart, min_value, max_value 

def create_smallbox(SnapLow, SnapHigh, num_chunks, snapdir, outdir):

    for snap in range(SnapLow, SnapHigh + 1):
        global_min = 1000 
        global_max = -999
        N = 0
        for chunk in tqdm(range(num_chunks)):
            numpart, min_value, max_value = scale_positions(snap, chunk, 
                                                            snapdir, outdir)
            N += numpart
            if (min_value < global_min):
                global_min = min_value
            if (max_value > global_max):
                global_max = max_value

        print("For snapshot {0} there were {1} Highres particles.  The "
              "minimum spatial position was {2} and the maximum was {3}."\
              .format(snap, N, global_min, global_max))

def adjust_header(SnapLow, SnapHigh, num_chunks, outdir):

    for snap in range(SnapLow, SnapHigh + 1):
        for chunk in tqdm(range(num_chunks)):

            fout_fname = "{0}/snapdir_{1:03d}/snapshot_{1:03d}.{2}.hdf5".format(outdir, snap, chunk) 
            with h5py.File(fout_fname, "r+") as f:

                tot_part = f['Header'].attrs['NumPart_Total']
                this_part = f['Header'].attrs['NumPart_ThisFile']
                
                for i in range(0,6):
                    if i == 1:
                        continue
                    tot_part[i] = 0 
                    this_part[i] = 0 

                f['Header'].attrs.modify('NumPart_Total', tot_part)
                f['Header'].attrs.modify('NumPart_ThisFile', this_part)

if __name__ == "__main__":

    SnapLow = 40
    SnapHigh = 40
    num_chunks = 256 
    snapdir = "/lustre/projects/p134_swin/jseiler/simulations/1.6Gpc/means/halo_1721673/dm_gadget/data"
    outdir = "/lustre/projects/p004_swin/jseiler/britton_smallbox"

    #create_smallbox(SnapLow, SnapHigh, num_chunks, snapdir, outdir)
    adjust_header(SnapLow, SnapHigh, num_chunks, outdir) 
