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

snaplow = 8
snaphigh = 10
num_cores = 256
filepath = '/lustre/projects/p004_swin/bsmith/1.6Gpc/means/halo_1721228/dm_gadget/data/'
TypeMax = 6

def write_fof_header(snapshot_idx):

    tmp = "snapdir_000/snapshot_000.0.hdf5"
    fname = filepath + tmp

#    mass_table = np.zeros((TypeMax), dtpye = np.float32)

    with h5py.File(fname, 'r') as f: 
        mass_table = f['Header'].attrs['MassTable']
       
    
    numpart_total = np.zeros((TypeMax), dtype = np.int64)
    for core_idx in range(0, num_cores):
        numpart_thisfile = np.zeros((TypeMax), dtype = np.int32)
        #tmp = "groups_%03d/fof_tab_%03d.%d.hdf5" %(snapshot_idx, snapshot_idx, core_idx) 
        #fname = filepath + tmp
 
        tmp = "/lustre/projects/p134_swin/jseiler/tinkering/groups_045/fof_tab_045.%d.hdf5" %(core_idx) # Directory for writable files.
        fname = tmp 

        with h5py.File(fname, 'r+') as f:
            try: 
                grouplen = f['Group']['GroupLenType']
            except KeyError:
                pass
            else:
                local_numpart = np.sum(grouplen[:], axis = 0) # Sums each particle type for this file.
                local_numpart.astype(np.int32) # Recast as explicit integer
                numpart_total = np.add(numpart_total, local_numpart)

            f['Header'].attrs.create("NumPart_ThisFile", local_numpart, dtype = np.int32)
            f['Header'].attrs.create("MassTable", mass_table, dtype = np.float32)
            #print "Done ThisFile stuff for Core %d" %(core_idx)

    ## At this point numpart_total contains the total number of particles (of each type) for the groups within the snapshot. ## 

    for core_idx in range(0, num_cores):
        numpart_highword = np.zeros((TypeMax), dtype = np.int32)
        #tmp = "groups_%03d/fof_tab_%03d.%d.hdf5" %(snapshot_idx, snapshot_idx, core_idx) 
        #fname = filepath + tmp
 
        tmp = "/lustre/projects/p134_swin/jseiler/tinkering/groups_045/fof_tab_045.%d.hdf5" %(core_idx) # Directory for writable files.
        fname = tmp 

        with h5py.File(fname, 'r+') as f:
            f['Header'].attrs.create("NumPart_Total", numpart_total, dtype = np.int32) # Least significant 32 bits of total number of particles.
            
            for type_idx in range(0, TypeMax):
                if numpart_total[type_idx] > pow(2, 32) - 1:
                    numpart_highword[type_idx] = numpart_total[type_idx] >> 32
            
            f['Header'].attrs.create("NumPart_Total_HighWord", numpart_highword, dtype = np.int32) # Most significant bits of total number of particles (if Number of particles is greater than 2^32 - 1). 
            #print "Done Total stuff for Core %d" %(core_idx)                

def write_fof_groups(snapshot_idx):


    for core_idx in range(0, num_cores):
        numpart_thisfile = np.zeros((TypeMax), dtype = np.int32)
        #tmp = "groups_%03d/fof_tab_%03d.%d.hdf5" %(snapshot_idx, snapshot_idx, core_idx) 
        #fname = filepath + tmp
 
        tmp = "/lustre/projects/p134_swin/jseiler/tinkering/groups_045/fof_tab_045.%d.hdf5" %(core_idx) # Directory for writable files.
        fname = tmp 

        with h5py.File(fname, 'r+') as f:
            try: 
                grouptype = f['IDs']['MemberType']
                grouplen = f['Group']['GroupLenType']
            except KeyError:
                pass
            else:
                for type_idx in range(0, TypeMax):
                    idx = np.where(grouptype[:] == type_idx)[0].tolist() # Need to convert this to a list so h5py can handle slicing.
                    if len(idx) == 0:
                        continue
                   
                    #print len(idx) 
                    pos = f['IDs']['MemberPos'][idx]
                    vel = f['IDs']['MemberVel'][idx]

                    groupID = np.zeros((len(idx)), dtype = np.int32)
                    offset = 0
                    for groupID_idx in range(0, np.shape(grouplen)[0]): # Loops over the number of FoF groups are in this file.
                        groupID[offset : grouplen[groupID_idx, type_idx] + offset] = groupID_idx # This assigns the number of particles of this type a group ID. 
                        #print grouplen[groupID_idx, type_idx]
                        #print groupID_idx
                        #print offset 
                        offset += grouplen[groupID_idx, type_idx]  
                    #print groupID
                    exit()

    

def write_snapshot_header(snapshot_idx):

    '''
    numpart_total = np.zeros((TypeMax), dtype = np.int64)
    for core_idx in range(0, num_cores):    
        numpart_thisfile = np.zeros((TypeMax), dtype = np.int32)

        tmp = "snapdir_%03d/snapshot_%03d.%d.hdf5" %(snapshot_idx, snapshot_idx, core_idx)
        fname = filepath + tmp
 
        with h5py.File(fname, 'r') as f:
            for type_idx in range(0, TypeMax):
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

    
 
def link_fof_snapshot(snapshot_idx):

    tmp = "groups_{0:03d}/fof_tab_{0:03d}.0.hdf5".format(snapshot_idx)
    fname = filepath + tmp 

    ds = yt.load(fname) # Loads in the entire FoF catalogue for this snapshot.
    ad = ds.all_data() # Container for all the data that we can work with.

    num_groups = ad['Group'].shape[0] # Number of groups within this snapshot.
    num_particles_per_group = np.zeros((num_groups), dtype = np.int32)
    print("For snapshot {0} there are {1} groups.".format(snapshot_idx, num_groups))

    ## Now the data format of the snapshot is split into the ParticleType groups. ##
    ## E.g. The IDs of only particle type 1 is stored in a 'PartType1' group. ##
    ## Hence we don't want to try and be searching for the ID of a particle type 2 within the particle type 1 group. ##
    ## So we need to store information about the particles in length 'TypeMax' lists. ##

    groupid = [[] for x in range(TypeMax)] # Array to hold the FoF Group ID for each group.
    fof_partid = [[] for x in range(TypeMax)] # Array to hold the Particle ID for each particle within each group.
    
    for group_idx in range(0, num_groups): # Loop over each group.
        halo = ds.halo('Group', group_idx) # Loads all the information for the specified group.
        num_particles_per_group[group_idx] = halo['member_ids'].shape[0]

        for particle_idx in range(0, int(num_particles_per_group[group_idx])): # Loop over the number of particles within this group.
            parttype = int(halo['MemberType'][particle_idx])
            groupid[parttype].append(group_idx)
            fof_partid[parttype].append(np.int64(halo['member_ids'][particle_idx]))
        print("Done group {0}".format(group_idx))

    # Our cross-matching requires the FoF particle IDs to be sorted. 
    # Need to be careful to sort the FoF group ID properly as well.
    for type_idx in range(0, TypeMax):
        groupid[type_idx] = [x for _,x in sorted(zip(fof_partid[type_idx], groupid[type_idx]))] # Taken from https://stackoverflow.com/questions/6618515/sorting-list-based-on-values-from-another-list
        fof_partid[type_idx] = sorted(fof_partid[type_idx])
    
    ## At this point we now have the particle IDs for each of the groups within the snapshot stored. ##
    ## Now need to load in the snapshot chunk by chunk and search for the particle IDs. ##
    ## Note: Since we want to do the snapshot loading piecewise (so we don't need ~250Gb ram) we won't use yt. ##

    for core_idx in range(0, num_cores): 
        tmp = "snapdir_{0:03d}/snapshot_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx)
        fname = filepath + tmp 


        with h5py.File(fname, 'r') as f:
            for type_idx in range(0, TypeMax):
                tmp = 'PartType{0}'.format(type_idx)
                try:
                    particles = f[tmp]
                except KeyError:
                    pass
                else:
                    snapshot_partid = particles['ParticleIDs']
                    print("Our key list (FoF Particle IDs) has length {0} and our match list (Snapshot Particle IDs) has length {1}".format(len(fof_partid[type_idx]), len(snapshot_partid)))

                    for particle_idx in range(0, len(snapshot_partid)): # Check every particle in this snapshot chunk.
                        particle_found = 0
                        count = 0
                        
                        if (particle_idx % 100000 == 0):
                            print(particle_idx) 
                        search_idx = int(np.ceil(len(fof_partid[type_idx]) / 2))
                        while(particle_found == 0):
                            count += 1 
                            if (snapshot_partid[particle_idx] == fof_partid[type_idx][search_idx]):
                                ids_found.append(snapshot_partid[particle_idx])
                                particle_found = 1
                            elif len(fof_partid[type_idx]) / pow(2, count) < 1.0:
                                break
                            elif (snapshot_partid[particle_idx] > fof_partid[type_idx][search_idx]):
                                search_idx = int(search_idx + np.ceil(len(fof_partid[type_idx]) / pow(2, count + 1)))
                            else:
                                search_idx = int(search_idx - np.ceil(len(fof_partid[type_idx]) / pow(2, count + 1)))                               
                            if search_idx == len(fof_partid[type_idx]):
                                search_idx = len(fof_partid[type_idx])
                        if particle_found == 1:
                            print("Found particle ID {0} after {1} loops".format(snapshot_partid[particle_idx], count)) 
                    #ids_found = np.intersect1d(list(particles['ParticleIDs']), partid[type_idx])
                    print(ids_found)
        exit()

    exit()
if __name__ == '__main__':

    for snapshot_idx in range(snaplow, snaphigh + 1):  
        
        #write_fof_header(snapshot_idx) 
        #write_fof_groups(snapshot_idx)
        #write_snapshot_header(snapshot_idx)
        link_fof_snapshot(snapshot_idx) 
