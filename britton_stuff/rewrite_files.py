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

snaplow = 7
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

    if len(ds.field_list) == 0: # If there aren't any groups in this file, then we skip all the linking.  Note we still need to create files to say there aren't any groups.
        num_groups = 0 
    else:
        num_groups = ad['Group'].shape[0] # Number of groups within this snapshot.
        num_particles_per_group = np.zeros((num_groups), dtype = np.int32)
        print("For snapshot {0} there are {1} groups.".format(snapshot_idx, num_groups))

        particle_position = [[] for x in range(TypeMax)]
        particle_velocity = [[] for x in range(TypeMax)]
        particle_id = [[] for x in range(TypeMax)]
        particle_groupid = [[] for x in range(TypeMax)]

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

    ## At this point we now have the particle IDs for each of the groups within the snapshot stored. ##
    ## Now need to load in the snapshot chunk by chunk and search for the particle IDs. ##
    ## Note: Since we want to do the snapshot loading piecewise (so we don't need ~250Gb ram) we won't use yt. ##

#    cores = [118]
    for core_idx in range(0 + rank, num_cores, size): 
    #for core_idx in cores: 
        numpart_thisfile = np.zeros((TypeMax), dtype = np.int32)
 
        tmp = "snapdir_{0:03d}/snapshot_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx)
        fname = filepath + tmp 

        print("Doing chunk {0}".format(core_idx))
        with h5py.File(fname, 'r') as f: 

            ## First grab some attributes that will be used for this snapshot. ##

            BoxSize = f['Header'].attrs['BoxSize'] # Size of the simulation box (in Mpc/h).
            Time = f['Header'].attrs['Time']            
            Omega_m = f['Header'].attrs['Omega0'] # Matter density.
            Omega_l = f['Header'].attrs['OmegaLambda'] # Dark energy density.
            particle_mass = f['Header'].attrs['MassTable'] # Table of particle masses (length equal to TypeMax)

            if num_groups > 0:
                for type_idx in range(0, TypeMax):
                    if len(fof_partid[type_idx]) > 0:
                        print("Doing PartType {0}".format(type_idx))
                        tmp = 'PartType{0}'.format(type_idx)
                        try:
                            particles = f[tmp]
                        except KeyError:
                            pass
                        else:
                            snapshot_partid = particles['ParticleIDs']
                            
                            print("Our key list (FoF Particle IDs) has length {0} and our match list (Snapshot Particle IDs) has length {1}".format(len(fof_partid[type_idx]), len(snapshot_partid)))
         
                            ids_found = np.nonzero(np.in1d(snapshot_partid, fof_partid[type_idx], assume_unique = True)) # np.in1d returns a True if the snapshot particle is found in the FoF list (False otherwise).  np.nonzero then returns only 'Trues'.  
                                                                                                                         # Hence 'ids_found' will be the snapshot particle INDICES for those particles in the FoF group.
                                                                                                                         # Taken from https://stackoverflow.com/questions/11483863/python-intersection-indices-numpy-array.
                            ids_found = ids_found[0].tolist() # Need to explicitly recast as a list for work with h5py.
                            ## Now we've found any matching IDs we need to grab the particles position/velocity. ##
                            if len(ids_found) > 0: # Checks to see if any matching IDs were found.
                                print(ids_found)                            
                                particle_position[type_idx].append(particles['Coordinates'][ids_found]) # For some super weird reason this needs 3 indices to properly reference.  I.e. particle_position[1][0] will return all the coordinates for Particle Type 1.
                                particle_velocity[type_idx].append(particles['Velocities'][ids_found])
                                particle_id[type_idx].append(particles['ParticleIDs'][ids_found])

                                ## Finally we need to grab the correct FoF group ID so we need to know which Snapshot Particle matched with which FoF Particle.  Then we use the index of where the FoF Particle is to grab the correct groupid. ##
                                fof_partid_idx = np.int32(np.nonzero(np.in1d(fof_partid[type_idx], snapshot_partid, assume_unique = True))[0]) 

                                groupids_added = []
                                sum_parts = 0
                                for i in fof_partid_idx: 
                                    particle_groupid[type_idx].append(groupid[type_idx][i])
                                    if (groupid[type_idx][i] in groupids_added) == False:
                                        groupids_added.append(groupid[type_idx][i]) 
                                        sum_parts += num_particles_per_group[groupid[type_idx][i]]
                                        print("Group {0} has {1} particles.".format(groupid[type_idx][i], num_particles_per_group[groupid[type_idx][i]]))

                                numpart_thisfile[type_idx] += len(ids_found)
                                numpart_allfiles[type_idx] += len(ids_found)
     
                                print("The number of particles that were found is {0}.  The number of particles in the groups in this chunk (groups {1}) is {2}".format(len(ids_found), groupids_added, sum_parts)) 
                                #assert(len(ids_found) == sum_parts)
        ## At this point we have successfully linked all the particles in this Snapshot chunk to the FoF Groups (this can, and will at high z, be zero particles). ##
        ## Now we need to construct a hdf5 file for this chunk, WRITE OUT THE HEADER, and then write out all of the particle data. ##
        ## Note: Since an important piece of information for the header is the total number of FoF Particles within the group, we will need to write out the individual parts of the header, keep a running total of the number of particles, then loop through N_cores again to write out the cumulative information. ##

        fname = "/lustre/projects/p134_swin/jseiler/simulations/my_britton_fof/groups_{0:03d}/my_fof_tab_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx)
       
        with h5py.File(fname, 'w') as f:
            f.create_group("Header") 
    
            ## Write local header attributes. ##

            f['Header'].attrs.create("NumFilePerSnapshot", num_cores, dtype = np.int32)
            f['Header'].attrs.create("BoxSize", BoxSize, dtype = np.float32)
            f['Header'].attrs.create("Time", Time, dtype = np.float32)
            f['Header'].attrs.create("Omega0", Omega_m, dtype = np.float32)
            f['Header'].attrs.create("OmegaLambda", Omega_l, dtype = np.float32)
            f['Header'].attrs.create("MassTable", particle_mass, dtype = np.float32)
            
            f['Header'].attrs.create("NumPart_ThisFile", numpart_thisfile, dtype = np.int32) 

            if num_groups > 0:
                for type_idx in range(0, TypeMax):
                    if len(particle_position[type_idx]) == 0:
                        continue
                    name = 'PartType{0}'.format(type_idx)
                    f.create_group(name)
                    
                    name = 'PartType{0}/Coordinates'.format(type_idx)
                    dset = f.create_dataset(name, dtype = np.float32, data = particle_position[type_idx][0]) # Position of each particle (Mpc/h).

                    name = 'PartType{0}/Velocities'.format(type_idx) 
                    dset = f.create_dataset(name, dtype = np.float32, data = particle_velocity[type_idx][0]) # Velocity of each particle (km/s).

                    name = 'PartType{0}/ParticleIDs'.format(type_idx) 
                    dset = f.create_dataset(name, dtype = np.int64, data = particle_id[type_idx][0]) # ID of each particle.

                    name = 'PartType{0}/GroupNumber'.format(type_idx) 
                    dset = f.create_dataset(name, dtype = np.int64, data = particle_groupid[type_idx]) # FoF ID of each particle.
                
            print("Written Data.")

    ## Here we have looped over all the cores. ##  
    ## Now need to open up each file one last time and write the Total Number of Particles attribute. ##

    comm.Barrier()

    if rank == 0:
        numpart_allfiles_total = np.zeros_like(numpart_allfiles)
    else:
        numpart_allfiles_total = None 
        
    comm.Reduce([numpart_allfiles, MPI.DOUBLE], [numpart_allfiles_total, MPI.DOUBLE], op = MPI.SUM, root = 0) # Sum all the stellar mass and pass to Rank 0.

    if rank == 0:
        for core_idx in range(0, num_cores):
#        for core_idx in cores:
                    
            fname = "/lustre/projects/p134_swin/jseiler/simulations/my_britton_fof/groups_{0:03d}/my_fof_tab_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx)
           
            with h5py.File(fname, 'r+') as f:

                f['Header'].attrs.create("NumPart_Total", numpart_allfiles_total, dtype = np.int32) # Least significant 32 bits of total number of particles.
                numpart_highword = np.zeros((TypeMax), dtype = np.int32)                

                for type_idx in range(0, TypeMax):
                    if numpart_allfiles_total[type_idx] > pow(2, 32) - 1:
                        numpart_highword[type_idx] = numpart_allfiles_total[type_idx] >> 32
                
                f['Header'].attrs.create("NumPart_Total_HighWord", numpart_highword, dtype = np.int32) # Most significant bits of total number of particles (if Number of particles is greater than 2^32 - 1). 
    
if __name__ == '__main__':

    for snapshot_idx in range(snaplow, snaphigh + 1):  
        
        #write_fof_header(snapshot_idx) 
        #write_fof_groups(snapshot_idx)
        #write_snapshot_header(snapshot_idx)
        link_fof_snapshot(snapshot_idx) 
