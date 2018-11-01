#!/usr/bin/env python
from __future__ import print_function

import sys
print(sys.path)

import matplotlib
matplotlib.use('Agg')

import h5py
import numpy as np
import pylab as plt
from random import sample

import yt

sys.path.append('/home/jseiler/SAGE-stuff/output/')

import AllVars
import PlotScripts
import ReadScripts

from mpi4py import MPI

from tqdm import tqdm, trange
comm= MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

num_cores = 256
groupdir = '/lustre/projects/p004_swin/bsmith/1.6Gpc/means/halo_1721673/dm_gadget/data/'
snapdir = '/lustre/projects/p134_swin/jseiler/simulations/1.6Gpc/means/halo_1721673/dm_gadget/data/'
linking_outdir = '/lustre/projects/p134_swin/jseiler/simulations/my_britton_fof_full/'
TypeMax = 6

subfind_dir = '/lustre/projects/p134_swin/jseiler/subfind_britton/'

def init(gridsize):
    
    local_grid = np.zeros((pow(gridsize,3)), dtype=np.float64)

    return local_grid


def map_particles(local_grid, fname):

    with h5py.File(fname, "r") as file_snap:
        


if __name__ == '__main__':

    # Initialize
    # Loop over core_idx
        # Grid snapshot chunk
    # Pass arrays back to root process
    # Normalize
    # Print

    ## Call arguments will be: Snapshot directory (snapshot base name), grid size, output directory (with base name) ##   
 
    if len(sys.argv) != 5:
        print("This program will construct a density grid for a GADGET HDF5 file.") 
        print("Usage: grid_density <snapshot directory (with base name)> <grid size> <output directory> <output base name>")
        exit()

    snapdir = sys.argv[1]
    gridsize = int(sys.argv[2])
    outputdir = sys.argv[3]
    out_name = sys.argv[4]

    if not os.path.exists(outputdir):
        print("{0} directory does not exist, creating it.".format(outputdir))
        os.makedirs(outputdir)

    local_grid = init(gridsize) # Initialize the local grid.
 
    fname = "{0}0.hdf5".format(snapdir)
    with h5py.File(fname, "r") as file_snap:
        num_files = file_snap['Header'].attrs['NumFilesPerSnapshot']
         
    print("There are {0} files per snapshots.".format(num_files))

    if(size > num_files):
        print("You are attempting to run with {0} processors. The number of files per snapshots is only {1}".format(size, num_files))
        exit()
    
    for core_idx in range(rank, num_files, size):
 
        print("I am rank {0} and I am reading chunk {1}.".format(rank, core_idx))
        fname = "{0}{1}.hdf5".format(snapdir)

        local_grid = map_particles(local_grid, fname)


