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

comm= MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

snaplow = 0
snaphigh = 30
num_cores = 256
#filepath = '/lustre/projects/p004_swin/bsmith/1.6Gpc/means/halo_1721228/dm_gadget/data/'
filepath = '/lustre/projects/p004_swin/bsmith/1.6Gpc/means/halo_1721673/dm_gadget/data/'
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

    numpart_allfiles = np.zeros((TypeMax), dtype = np.int64) 

    groupid = [[] for x in range(TypeMax)] # Array to hold the FoF Group ID for each group.
    fof_partid = [[] for x in range(TypeMax)] # Array to hold the Particle ID for each particle within each group.
    particle_groupid = [[] for x in range(TypeMax)] # Array to hold the FoF Group ID for the matched snapshot particles.    
 
    if len(ds.field_list) == 0: # If there aren't any groups in this file, then we skip all the linking.  Note we still need to create files to say there aren't any groups.
        num_groups = 0 
    else:
        num_groups = ad['Group'].shape[0] # Number of groups within this snapshot.
        num_particles_per_group = np.zeros((num_groups), dtype = np.int32)
        print("For snapshot {0} there are {1} groups.".format(snapshot_idx, num_groups))

        ## Now the data format of the snapshot is split into the ParticleType groups. ##
        ## E.g. The IDs of only particle type 1 is stored in a 'PartType1' group. ##
        ## Hence we don't want to try and be searching for the ID of a particle type 2 within the particle type 1 group. ##
        ## So we need to store information about the particles in length 'TypeMax' lists. ##
 
        for group_idx in range(0, num_groups): # Loop over each group.
            halo = ds.halo('Group', group_idx) # Loads all the information for the specified group.
            num_particles_per_group[group_idx] = halo['member_ids'].shape[0]

            for particle_idx in range(0, int(num_particles_per_group[group_idx])): # Loop over the number of particles within this group.
                parttype = int(halo['MemberType'][particle_idx])
                groupid[parttype].append(group_idx)
                fof_partid[parttype].append(np.int64(halo['member_ids'][particle_idx]))
            print("Done group {0}".format(group_idx))

    ## At this point we now have the particle IDs for each of the groups within the snapshot stored. ##
    ## Now need to load in the snapshot chunk by chunk and search for the particle IDs. ##
    ## Note: Since we want to do the snapshot loading piecewise (so we don't need ~250Gb ram) we won't use yt. ##

#    cores = [118]
    for core_idx in range(0 + rank, num_cores, size): 
    #for core_idx in cores: 
        snapshot_groupid = [[] for x in range(TypeMax)] # Array to hold the FoF Group ID for the particles within each snapshot.
        tmp = "snapdir_{0:03d}/snapshot_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx)
        fname = filepath + tmp 

        print("Doing chunk {0}".format(core_idx))
        with h5py.File(fname, 'r') as f:

            NumPart_ThisFile = f['Header'].attrs['NumPart_ThisFile']
 
            for type_idx in range(0, TypeMax):
                print("Doing PartType {0}".format(type_idx))
                tmp = 'PartType{0}'.format(type_idx)


                try:
                    particles = f[tmp]
                except KeyError:
                    pass
                else:
                    if len(fof_partid[type_idx]) > 0:
                        snapshot_partid = particles['ParticleIDs']
                        
                        print("Our key list (FoF Particle IDs) has length {0} and our match list (Snapshot Particle IDs) has length {1}".format(len(fof_partid[type_idx]), len(snapshot_partid)))
     
                        ids_found = np.nonzero(np.in1d(snapshot_partid, fof_partid[type_idx], assume_unique = True)) # np.in1d returns a True if the snapshot particle is found in the FoF list (False otherwise).  np.nonzero then returns only 'Trues'.  
                                                                                                                     # Hence 'ids_found' will be the snapshot particle INDICES for those particles in the FoF group.
                                                                                                                     # Taken from https://stackoverflow.com/questions/11483863/python-intersection-indices-numpy-array.
                        ids_found = ids_found[0].tolist() # Need to explicitly recast as a list for work with h5py.

                        ## We need to grab the correct FoF group ID so we need to know which Snapshot Particle matched with which FoF Particle.  Then we use the index of where the FoF Particle is to grab the correct groupid. ##
                        if len(ids_found) > 0: # Checks to see if any matching IDs were found.
                            print(ids_found)                            
                            particle_id[type_idx].append(particles['ParticleIDs'][ids_found])                         
                            fof_partid_idx = np.int32(np.nonzero(np.in1d(fof_partid[type_idx], snapshot_partid, assume_unique = True))[0]) 

                            groupids_added = []
                            sum_parts = 0
                            for i in fof_partid_idx: 
                                particle_groupid[type_idx].append(groupid[type_idx][i])
                                if (groupid[type_idx][i] in groupids_added) == False:
                                    groupids_added.append(groupid[type_idx][i]) 
                                    sum_parts += num_particles_per_group[groupid[type_idx][i]]
                                    print("Group {0} has {1} particles.".format(groupid[type_idx][i], num_particles_per_group[groupid[type_idx][i]]))

                            print("The number of particles that were found is {0}.  The number of particles in the groups in this chunk (groups {1}) is {2}".format(len(ids_found), groupids_added, sum_parts)) 

                   
                    
                    particle_fof_id = np.full((NumPart_ThisFile[type_idx]), 1<<30, dtype = np.int32) # First initialize every Snapshot Particle to be 'Not in a group' (1<<30 is HBT+'s flag for this).
                    if len(particle_groupid[type_idx]) > 0: # Then loop through each of the matched particles,
                        for i in range(0, len(ids_found)):
                            particle_fof_id[ids_found[i]] = particle_groupid[i] # Update the GroupID value for the Snapshot Particles to be the correct FoF GroupID that we matched previously.

                    snapshot_groupid[type_idx].append(particle_fof_id)
                    snapshot_groupid[type_idx] = snapshot_groupid[type_idx][0] # Fix up the indexing of this array.  The result of this is 'snapshot_partid[type_idx]' returns an array (not a nested one). 

        ## At this point we have matched any particles in the snapshot that are in FoF Groups. ##
        ## We have also constructed an array (snapshot_groupid) that contains the FoF Group for each snapshot particle; using '1<<30' if the particle is not in a group. ##

        fname = "/lustre/projects/p134_swin/jseiler/simulations/my_britton_fof/groups_{0:03d}/my_fof_groupids_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx) 

        with h5py.File(fname, 'w') as f:
            for type_idx in range(0, TypeMax):
                if NumPart_ThisFile[type_idx] == 0:
                    continue
                
                name = 'PartType{0}'.format(type_idx)
                f.create_group(name)

                name = 'PartType{0}/GroupNumber'.format(type_idx) 
                dset = f.create_dataset(name, dtype = np.int32, data = snapshot_groupid[type_idx]) # Velocity of each particle (km/s).

        print("Written data to {0}".format(fname))
        
if __name__ == '__main__':

    for snapshot_idx in range(snaplow, snaphigh + 1):  
        
        #write_fof_header(snapshot_idx) 
        #write_fof_groups(snapshot_idx)
        #write_snapshot_header(snapshot_idx)
        link_fof_snapshot(snapshot_idx) 
