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

from matplotlib.patches import Circle
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

snaplow = 20 
snaphigh = 20
num_cores = 256 
groupdir = '/lustre/projects/p004_swin/bsmith/1.6Gpc/means/halo_1721673/dm_gadget/data/'
snapdir = '/lustre/projects/p134_swin/jseiler/simulations/1.6Gpc/means/halo_1721673/dm_gadget/data/'
linking_outdir = '/lustre/projects/p134_swin/jseiler/simulations/my_britton_fof_full/'
TypeMax = 6

subfind_dir = '/lustre/projects/p134_swin/jseiler/subfind_britton/' 

bin_width = 0.1
 
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

    
 
def link_fof_snapshot_ids(snapshot_idx):
    '''
    This function goes through each snapshot within the simulation and matches the FoF group number to each particle.  
    This matching is done using the 'fof_tab' files. If a particle is not within a FoF group (i.e. it is 'unbound') it is assigned a special value of 1 << 30 = 1073741824.
    Note: This function is different to link_fof_snapshots_full() in that it only saves a file containing the groupID and particleID, no other properties are saved.

    Parameters 
    ----------
    snapshot_idx : int
        Snapshot number that we are doing the linking for.

    Returns
    -------
    No returns.
    
    Units
    -----
    All units are kept in internal units for the simulation.
    ''' 
    

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
 
        for group_idx in trange(num_groups): # Loop over each group.
            halo = ds.halo('Group', group_idx) # Loads all the information for the specified group.
            num_particles_per_group[group_idx] = halo['member_ids'].shape[0]

            for particle_idx in range(0, int(num_particles_per_group[group_idx])): # Loop over the number of particles within this group.
                parttype = int(halo['MemberType'][particle_idx])
                groupid[parttype].append(group_idx)
                fof_partid[parttype].append(np.int64(halo['member_ids'][particle_idx]))
            #print("Done group {0}".format(group_idx))

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
                #print("Doing PartType {0}".format(type_idx))
                tmp = 'PartType{0}'.format(type_idx)

                try:
                    particles = f[tmp]
                except KeyError:
                    pass
                else:
                    snapshot_partid = particles['ParticleIDs']
                    if len(fof_partid[type_idx]) > 0:
                
                        
                        #print("Our key list (FoF Particle IDs) has length {0} and our match list (Snapshot Particle IDs) has length {1}".format(len(fof_partid[type_idx]), len(snapshot_partid)))
     
                        ids_found = np.nonzero(np.in1d(snapshot_partid, fof_partid[type_idx])) # np.in1d returns a True if the snapshot particle is found in the FoF list (False otherwise).  np.nonzero then returns the indices of those with a value 'True'.  
                                                                                                                     # Hence 'ids_found' will be the snapshot particle INDICES for those particles in the FoF group.
                                                                                                                     # Taken from https://stackoverflow.com/questions/11483863/python-intersection-indices-numpy-array.
                        ids_found = ids_found[0].tolist() # Need to explicitly recast as a list for work with h5py.

                        ## We need to grab the correct FoF group ID so we need to know which Snapshot Particle matched with which FoF Particle.  Then we use the index of where the FoF Particle is to grab the correct groupid. ##
                        if len(ids_found) > 0: # Checks to see if any matching IDs were found.                         
                            fof_partid_idx = np.int32(np.nonzero(np.in1d(fof_partid[type_idx], snapshot_partid, assume_unique = True))[0]) 

                            groupids_added = []
                            sum_parts = 0
                            for i in fof_partid_idx: 
                                particle_groupid[type_idx].append(groupid[type_idx][i])
                                if (groupid[type_idx][i] in groupids_added) == False:
                                    groupids_added.append(groupid[type_idx][i]) 
                                    sum_parts += num_particles_per_group[groupid[type_idx][i]]

                            #print("The number of particles that were found is {0}.  The number of particles in the groups in this chunk (groups {1}) is {2}".format(len(ids_found), groupids_added, sum_parts)) 
                            print("The number of particles that were found is {0}.".format(len(ids_found)))

                    particle_fof_id = np.full((NumPart_ThisFile[type_idx]), 1<<30, dtype = np.int32) # First initialize every Snapshot Particle to be 'Not in a group' (1<<30 is HBT+'s flag for this).
                    if len(particle_groupid[type_idx]) > 0: # Then loop through each of the matched particles,
                        for i in range(0, len(ids_found)):

#                            print(ids_found[i])     
#                            print(particle_groupid[type_idx][i])
                            particle_fof_id[ids_found[i]] = particle_groupid[type_idx][i] # Update the GroupID value for the Snapshot Particles to be the correct FoF GroupID that we matched previously.

                    snapshot_groupid[type_idx].append(particle_fof_id)
                    snapshot_groupid[type_idx] = snapshot_groupid[type_idx][0] # Fix up the indexing of this array.  The result of this is 'snapshot_partid[type_idx]' returns an array (not a nested one). 


            ## At this point we have matched any particles in the snapshot that are in FoF Groups. ##
            ## We have also constructed an array (snapshot_groupid) that contains the FoF Group for each snapshot particle; using '1<<30' if the particle is not in a group. ##

            fname = "/lustre/projects/p134_swin/jseiler/simulations/my_britton_fof/groups_{0:03d}/my_fof_groupids_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx) 

            with h5py.File(fname, 'w') as f2:
                for type_idx in range(0, TypeMax):
                    if NumPart_ThisFile[type_idx] == 0:
                        continue
                    
                    name = 'PartType{0}'.format(type_idx)
                    particle_ids = f[name]['ParticleIDs']
                    f2.create_group(name)

                    name = 'PartType{0}/GroupNumber'.format(type_idx) 
                    dset = f2.create_dataset(name, dtype = np.int32, data = snapshot_groupid[type_idx]) # Velocity of each particle (km/s).

                    name = 'PartType{0}/ParticleID'.format(type_idx) 
                    dset = f2.create_dataset(name, dtype = np.int64, data = particle_ids) # Velocity of each particle (km/s).
                
                

        print("Written data to {0}".format(fname))

def link_fof_snapshot_full(snapshot_idx, add_unbound_particles):
    '''
    This function goes through each snapshot within the simulation and matches the FoF particles with those in the original snapshot list.  It then saves the FoF particles with the relevant properties to a separate file. 
    Note: This function is different to link_fof_snapshots_ids() in that it saves a file containing the full information regarding the FoF particle, not merely its ID/GroupID. 

    Parameters 
    ----------
    snapshot_idx : int
        Snapshot number that we are doing the linking for.
    add_unbound_particles : int
        If we wish to add a random 1000 unbound particles this flag should be set to '1'.  Otherwise we only include particles bound in FoFs.

    Returns
    -------
    No returns.
    
    Units
    -----
    All units are kept in internal units for the simulation.
    ''' 

    numpart_allfiles = np.zeros((TypeMax), dtype = np.int64)    

    if rank == 0:
        fof_partid = [[] for x in range(TypeMax)] # Array to hold the Particle ID for each particle within each group.

        for core_idx in trange(num_cores):
            fname = "{1}groups_{0:03d}/fof_tab_{0:03d}.{2}.hdf5".format(snapshot_idx, groupdir, core_idx)

            with h5py.File(fname, "r") as file_fof_tab:
                Ngroups_Total = file_fof_tab['Header'].attrs['Ngroups_Total']
                if(Ngroups_Total == 0):
                    break
               
                Ngroups_ThisFile = file_fof_tab['Header'].attrs['Ngroups_ThisFile']
                if(Ngroups_ThisFile == 0):
                    continue
                part_ids = file_fof_tab['IDs']['ID'][:]
                part_types = file_fof_tab['IDs']['MemberType'][:]

                for part_idx in range(len(part_ids)):
                    parttype = part_types[part_idx]
                    fof_partid[parttype].append(np.int64(part_ids[part_idx])) 


        print("For snapshot {0} there are {1} groups.".format(snapshot_idx, Ngroups_Total))
    else:
        fof_partid = None   
        Ngroups_Total = None

    print("I am about to broadcast everything.") 
    fof_partid = comm.bcast(fof_partid, root = 0)
    Ngroups_Total = comm.bcast(Ngroups_Total, root = 0)    
        
    '''
    ds = yt.load(fname) # Loads in the entire FoF catalogue for this snapshot.
    ad = ds.all_data() # Container for all the data that we can work with.

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

        groupid = [[] for x in range(TypeMax)] # Array to hold the FoF Group ID for each group.
        
        for group_idx in trange(num_groups): # Loop over each group.
            halo = ds.halo('Group', group_idx) # Loads all the information for the specified group.
            num_particles_per_group[group_idx] = halo['member_ids'].shape[0]

            for particle_idx in range(0, int(num_particles_per_group[group_idx])): # Loop over the number of particles within this group.
                parttype = int(halo['MemberType'][particle_idx])
                groupid[parttype].append(group_idx)
                fof_partid[parttype].append(np.int64(halo['member_ids'][particle_idx]))
        
    '''
    ## At this point we now have the particle IDs for each of the groups within the snapshot stored. ##
    ## Now need to load in the snapshot chunk by chunk and search for the particle IDs. ##
    ## Note: Since we want to do the snapshot loading piecewise (so we don't need ~250Gb ram) we won't use yt. ##

#    cores = [118]
    for core_idx in range(0 + rank, num_cores, size): 
    #for core_idx in cores:
 
        particle_position = [[] for x in range(TypeMax)]
        particle_velocity = [[] for x in range(TypeMax)]
        particle_id = [[] for x in range(TypeMax)]
        numpart_thisfile = np.zeros((TypeMax), dtype = np.int32)

        if snapshot_idx < 66: 
            tmp = "snapdir_{0:03d}/snapdir_{0:03d}/snapshot_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx)
        else:
            tmp = "snapdir_{0:03d}/snapshot_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx)
        fname = snapdir + tmp 

        print("Doing chunk {0}".format(core_idx))
        with h5py.File(fname, 'r') as f: 

            ## First grab some attributes that will be used for this snapshot. ##

            BoxSize = f['Header'].attrs['BoxSize'] # Size of the simulation box (in Mpc/h).
            Time = f['Header'].attrs['Time'] # Scale factor. 
            Omega_m = f['Header'].attrs['Omega0'] # Matter density.
            Omega_l = f['Header'].attrs['OmegaLambda'] # Dark energy density.
            particle_mass = f['Header'].attrs['MassTable'] # Table of particle masses (length equal to TypeMax)

            if Ngroups_Total > 0:
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
                            
#                            print("Our key list (FoF Particle IDs) has length {0} and our match list (Snapshot Particle IDs) has length {1}".format(len(fof_partid[type_idx]), len(snapshot_partid)))
         
                            ids_found = (np.nonzero(np.in1d(snapshot_partid, fof_partid[type_idx], assume_unique = True)))[0] # np.in1d returns a True if the snapshot particle is found in the FoF list (False otherwise).  np.nonzero then returns the indices of those values that have a 'True'.  
                                                                                                                         # Hence 'ids_found' will be the snapshot particle INDICES for those particles in the FoF group.
                                                                                                                         # Taken from https://stackoverflow.com/questions/11483863/python-intersection-indices-numpy-array.
                          
                            if(add_unbound_particles == 1 and len(ids_found) > 0):
                                ids_before = len(ids_found)  
                                random_snap_parts = sample(range(len(snapshot_partid)), 1000)
                                ids_found = np.append(ids_found, random_snap_parts)

                                print("We have added 1000 extra particles.")
                                ids_found = np.unique(ids_found) 
                                print("Ensuring that we have only unique particles we therefore have added {0} unbound particles.".format(len(ids_found) - ids_before))
                    

                            ids_found = ids_found.tolist() # Need to explicitly recast as a list for work with h5py.                                                                                             
                            ## Now we've found any matching IDs we need to grab the particles position/velocity. ##
                            if len(ids_found) > 0: # Checks to see if any matching IDs were found.

                                particle_position[type_idx].append(particles['Coordinates'][ids_found]) # For some super weird reason this needs 3 indices to properly reference.  I.e. particle_position[1][0] will return all the coordinates for Particle Type 1.
                                particle_velocity[type_idx].append(particles['Velocities'][ids_found])
                                particle_id[type_idx].append(particles['ParticleIDs'][ids_found])

                                ## Finally we need to grab the correct FoF group ID so we need to know which Snapshot Particle matched with which FoF Particle.  Then we use the index of where the FoF Particle is to grab the correct groupid. ##
                                fof_partid_idx = np.int32(np.nonzero(np.in1d(fof_partid[type_idx], snapshot_partid, assume_unique = True))[0]) 

                                numpart_thisfile[type_idx] += len(ids_found)
                                numpart_allfiles[type_idx] += len(ids_found)
                                print("I am rank {0} and I found {1} particles from snapshot {2} in chunk {3}".format(rank, len(ids_found), snapshot_idx, core_idx))   
 
        ## At this point we have successfully linked all the particles in this Snapshot chunk to the FoF Groups (this can, and will at high z, be zero particles). ##
        ## Now we need to construct a hdf5 file for this chunk, WRITE OUT THE HEADER, and then write out all of the particle data. ##
        ## Note: Since an important piece of information for the header is the total number of FoF Particles within the group, we will need to write out the individual parts of the header, keep a running total of the number of particles, then loop through N_cores again to write out the cumulative information. ##

        fname = "{2}groups_{0:03d}/my_fof_tab_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx, linking_outdir)
       
        with h5py.File(fname, 'w') as f:
            f.create_group("Header") 
    
            ## Write local header attributes. ##

            f['Header'].attrs.create("NumFilesPerSnapshot", num_cores, dtype = np.int32)
            f['Header'].attrs.create("BoxSize", BoxSize, dtype = np.float32)
            f['Header'].attrs.create("Time", Time, dtype = np.float32)
            f['Header'].attrs.create("Omega0", Omega_m, dtype = np.float32)
            f['Header'].attrs.create("OmegaLambda", Omega_l, dtype = np.float32)
            f['Header'].attrs.create("MassTable", particle_mass, dtype = np.float32)
            
            f['Header'].attrs.create("NumPart_ThisFile", numpart_thisfile, dtype = np.int32) 

            if Ngroups_Total > 0:
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

            print("Written Data to {0}.".format(fname))

    ## Here we have looped over all the cores. ## 
 
    ## Now need to open up each file one last time and write the Total Number of Particles attribute. ##

    print("I am Task {0} and I have finished my writing.  Now waiting on other tasks to finish.".format(rank))
    comm.Barrier()

    print("Now going through all the chunks and summing up the number of particles.")
    if rank == 0:
        numpart_allfiles_total = np.zeros_like(numpart_allfiles)
    else:
        numpart_allfiles_total = None 
        
    comm.Reduce([numpart_allfiles, MPI.DOUBLE], [numpart_allfiles_total, MPI.DOUBLE], op = MPI.SUM, root = 0) # Sum all the stellar mass and pass to Rank 0.

    if rank == 0:
        for core_idx in range(0, num_cores):
#        for core_idx in cores:
                    
            fname = "{2}groups_{0:03d}/my_fof_tab_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx, linking_outdir)
           
            with h5py.File(fname, 'r+') as f:

                f['Header'].attrs.create("NumPart_Total", numpart_allfiles_total, dtype = np.int32) # Least significant 32 bits of total number of particles.
                numpart_highword = np.zeros((TypeMax), dtype = np.int32)                

                for type_idx in range(0, TypeMax):
                    if numpart_allfiles_total[type_idx] > pow(2, 32) - 1:
                        numpart_highword[type_idx] = numpart_allfiles_total[type_idx] >> 32
                
                f['Header'].attrs.create("NumPart_Total_HighWord", numpart_highword, dtype = np.int32) # Most significant bits of total number of particles (if Number of particles is greater than 2^32 - 1). 

        print("Fully finished writing to snapshot {0}".format(snapshot_idx))

def check_linking_ids(snapshot_idx):

    for core_idx in range(0, num_cores):
        print("Doing chunk {0}".format(core_idx))

        fname_fof = "/lustre/projects/p134_swin/jseiler/simulations/my_britton_fof/groups_{0:03d}/my_fof_groupids_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx)        

        tmp = "snapdir_{0:03d}/snapshot_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx)
        fname_snap = filepath + tmp 
        with h5py.File(fname_fof, "r") as file_fof, h5py.File(fname_snap, "r") as file_snap:
            NumPart_ThisFile = file_snap['Header'].attrs['NumPart_ThisFile']
            for type_idx in range(TypeMax):
                if(NumPart_ThisFile[type_idx] == 0):
                    continue

                tmp = 'PartType{0}'.format(type_idx)
                group_ids = file_fof[tmp]['GroupNumber']
                fof_particleids = file_fof[tmp]['ParticleID']

                snap_particleids = file_snap[tmp]['ParticleIDs']

                for part_idx in range(len(snap_particleids)):
                    if(part_idx % 1000000 == 0):
                        print(part_idx)
                    if(snap_particleids[part_idx] != fof_particleids[part_idx]):
                        print("Snapshot ID = {0} \t FoF ID = {1}".format(snap_particleids[part_idx], fof_particleids[part_idx]))
                        assert(snap_particleids[part_idx] == fof_particleids[part_idx])


def check_linking_full(snapshot_idx, full_debug): 
    '''
    In this function we want to make sure that our list of linked particles has been correctly created.  We do this in the SIMPLEST and DUMBEST way possible and so there will be a 2 step process.
    1: Ensure that every single particle in the FoF tables have been included in my linked list output.
    2: Ensure that the properties of the particles in the linked list output is the same as the Snapshot properties.

    Only comments for new code has been included here.  Refer to other functions for comments on lines such as 'np.nonzero(...)'.

    Parameters 
    ----------
    snapshot_idx : int
        Snapshot number that we are doing the linking for.
    full_debug : int
        Flag for whether we want a full debug run which involves a lot of printing to screen.

    Returns
    -------
    No returns. The function will run until all particles are checked.
    
    Units
    -----
    All units are kept in internal units for the simulation.
    
    '''

    ## First checking all FoF particles have been included in the linked list. ##

    numpart_allfiles = np.zeros((TypeMax), dtype = np.int64)    

    if rank == 0:
        fof_tab_ids = []

        for core_idx in trange(num_cores):
            fname = "{1}groups_{0:03d}/fof_tab_{0:03d}.{2}.hdf5".format(snapshot_idx, groupdir, core_idx)

            with h5py.File(fname, "r") as file_fof_tab:
                Ngroups_Total = file_fof_tab['Header'].attrs['Ngroups_Total']
                if(Ngroups_Total == 0):
                    break
               
                Ngroups_ThisFile = file_fof_tab['Header'].attrs['Ngroups_ThisFile']
                if(Ngroups_ThisFile == 0):
                    continue
                part_ids = file_fof_tab['IDs']['ID'][:]

                for part_idx in range(len(part_ids)):
                    fof_tab_ids.append(np.int64(part_ids[part_idx])) 


        print("For snapshot {0} there are {1} groups.".format(snapshot_idx, Ngroups_Total))
    else:
        fof_partid = None   
        Ngroups_Total = None

    print("I am about to broadcast everything.") 
#    fof_partid = comm.bcast(fof_tab_ids, root = 0)
#    Ngroups_Total = comm.bcast(Ngroups_Total, root = 0)    


    '''
    if(ds == None):
        fname_fof_tab_ids = "{1}groups_{0:03d}/fof_tab_{0:03d}.0.hdf5".format(snapshot_idx, groupdir)
        print("Reading from file {0}".format(fname_fof_tab_ids)) 
        ds = yt.load(fname_fof_tab_ids) # Loads in the entire FoF catalogue for this snapshot.
        ad = ds.all_data() # Container for all the data that we can work with.

    fof_tab_ids = [] 

    # First load all the FoF particle IDs from Britton's FoF_tab files. #
    if len(ds.field_list) == 0: 
        num_groups = 0 
    else:
        num_groups = ad['Group'].shape[0] # Number of groups within this snapshot.
        num_particles_per_group = np.zeros((num_groups), dtype = np.int32)
        print("For snapshot {0} there are {1} groups.".format(snapshot_idx, num_groups))
 
        for group_idx in trange(num_groups): # Loop over each group.
            halo = ds.halo('Group', group_idx) # Loads all the information for the specified group.
            num_particles_per_group[group_idx] = halo['member_ids'].shape[0]

            for particle_idx in range(0, int(num_particles_per_group[group_idx])): # Loop over the number of particles within this group.
                parttype = int(halo['MemberType'][particle_idx])
                fof_tab_ids.append(np.int64(halo['member_ids'][particle_idx]))
    '''

    # Next we load in all the particle IDs within my linked lists. #
    print("Loaded in all the FoF groups.  Now loading the linked list.")
    linked_list_ids = []
    header_sum = np.zeros((TypeMax))
    header_maxpart = np.zeros((TypeMax))

    for core_idx in range(num_cores): 
        fname_linked_list_ids = "{2}groups_{0:03d}/my_fof_tab_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx, linking_outdir)
        print("Chunk {0}".format(core_idx))       
 
        with h5py.File(fname_linked_list_ids, "r") as file_linked_list: # Open up the linked list.

            if(Ngroups_Total == 0): # If there aren't any FoF groups at all for this snapshot, 
                assert(all(N == 0 for N in (file_linked_list['Header'].attrs['NumPart_ThisFile']))) # Confirm that there aren't any in the linked list,
                continue # Then move onto the next chunk (we will still check the rest of the chunks to make 100% sure).

            for type_idx in range(TypeMax):
                if(core_idx > 0):
                    if(full_debug == 1):
                        print("NumPart_Total from previous chunk = {0} \t NumPart_Total from current chunk = {1}".format(header_maxpart[type_idx], file_linked_list['Header'].attrs['NumPart_Total']))
                    assert(header_maxpart[type_idx] == file_linked_list['Header'].attrs['NumPart_Total'][type_idx])
                header_maxpart[type_idx] = file_linked_list['Header'].attrs['NumPart_Total'][type_idx] 

                name = 'PartType{0}'.format(type_idx)
                try:
                    linked_list_ids_tmp = file_linked_list[name]['ParticleIDs'][:]
                except KeyError:
                    pass
                else:
                    header_sum[type_idx] += file_linked_list['Header'].attrs['NumPart_ThisFile'][type_idx]
                                        
                    for i in range(len(linked_list_ids_tmp)):
                        linked_list_ids.append(linked_list_ids_tmp[i])

    for type_idx in range(TypeMax):
        assert(header_sum[type_idx] == header_maxpart[type_idx])
    
    print("The header of my linked_list has a total of {0} particles and the fof_tab_ids has length {1}".format(len(fof_tab_ids), sum(header_maxpart)))
    assert(sum(header_maxpart) == len(fof_tab_ids))
   
    # Now that we have both the particle IDs within the FoF groups (from Brittons) and the IDs of the particles within my linked list, ensure that ALL entries in Britton's groups are in my FoF groups. #
    print("The length of the linked list IDs is {0} and the length of the fof tab IDs is {1}".format(len(linked_list_ids), len(fof_tab_ids)))
    assert(len(linked_list_ids) == len(fof_tab_ids)) # Make sure the arrays are the same size.

    ids_found = (np.nonzero(np.in1d(fof_tab_ids, linked_list_ids)))[0] # Returns the indices of the matches between the FoF tab and linked list.
    assert(len(ids_found) == len(fof_tab_ids)) # Check to see if all the FoF tab ids can be found in the linked list.
        
    ids_found = (np.nonzero(np.in1d(linked_list_ids, fof_tab_ids)))[0] # Also want to check the reverse match is satisfied.
    assert(len(ids_found) == len(fof_tab_ids)) # Check to see if all the FoF tab ids can be found in the linked list.

    assert(len(np.unique(fof_tab_ids)) == len(fof_tab_ids)) # Ensure that each ID is only repeated once.
    assert(len(np.unique(linked_list_ids)) == len(linked_list_ids))

    print("Checked that all the IDs within the FoF_tab file are present within the output linked list.")

    ## We have now ensured that all the particles from the FoF tab are found in the linked output list. ##
    ## We now must check to ensure that all the properties of the linked output list are IDENTICAL to those found within the snapshot list. ##
    ## I do this in the stupidest way possible. ##

    # For every single chunk, I first find those particles from the FoF groups that are present within this chunk (i.e. identical to what I did for the linking procedure). # 

    # Check verification exactly once. #
    print("We are now checking the properties of each particle in the FoF linked list to ensure they match those found in the Snapshot.")
    for core_idx in range(0 + rank, num_cores, size):     
        print("Checking chunk {0}.".format(core_idx))
        if snapshot_idx < 66: 
            fname_snapshot = "{2}snapdir_{0:03d}/snapdir_{0:03d}/snapshot_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx, snapdir)
        else:
            fname_snapshot = "{2}snapdir_{0:03d}/snapshot_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx, snapdir)
#        fname_snapshot = "{2}snapdir_{0:03d}/snapshot_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx, snapdir) 
        
        with h5py.File(fname_snapshot, "r") as file_snapshot: # Open up the linked list. 
            if(Ngroups_Total == 0): # If there aren't any FoF groups at all for this snapshot, 
                break # Just keep going to the next snapshot then.
            for type_idx in range(TypeMax):
                print("Doing PartType {0}".format(type_idx))
                name = 'PartType{0}'.format(type_idx)
                try:
                    snapshot_partid = file_snapshot[name]['ParticleIDs'][:]
                except KeyError:
                    pass
                else:

                    if (full_debug == 1):
                        print("Length of snapshot IDs in this chunk is {0}.".format(len(snapshot_partid)))

                    ids_found = (np.nonzero(np.in1d(snapshot_partid, fof_tab_ids)))[0] 
                    snapshot_ids_found = snapshot_partid[ids_found]

                    print("We found {0} particles from this chunk that were also present in the FoF Linked List.".format(len(ids_found))) 
                    ## Now that we have matched any FoF particles in this Snapshot, let's ensure the stored velocities/positions are identical. ## 

                    if len(ids_found) > 0:
                        # snapshot_part_idx is ambiguous. #
                        # Zip + Enumerate # 
#                        for snapshot_part_idx, partid enumerate(snapshot_partid[ids_found]):

                        fname_linked_list = "{2}groups_{0:03d}/my_fof_tab_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx, linking_outdir)
                        with h5py.File(fname_linked_list, "r") as file_linked_list:
                            try:
                                linked_list_partid = file_linked_list[name]['ParticleIDs'][:]
                            except KeyError:
                                pass
                            else:
                                
                                linked_list_ids_found = (np.nonzero(np.in1d(linked_list_partid, snapshot_ids_found)))[0]
                                assert(len(linked_list_ids_found) == len(linked_list_partid))
                                assert(len(linked_list_ids_found) == len(snapshot_ids_found))

                                #print("Linked list z = {0}, Snapshot z = {1}".format(file_linked_list[name]['Coordinates'][linked_list_part_idx, 2], file_snapshot[name]['Coordinates'][ids_found[snapshot_part_idx]][2]))
                                for i in range(len(linked_list_ids_found)):                               
                                    if(i % 10000 == 0):
                                        print("Finished checking particle {0}".format(i)) 
                                    assert(file_linked_list[name]['Coordinates'][linked_list_ids_found[i], 0] == file_snapshot[name]['Coordinates'][ids_found[i], 0]) 
                                    assert(file_linked_list[name]['Coordinates'][linked_list_ids_found[i], 1] == file_snapshot[name]['Coordinates'][ids_found[i], 1])
                                    assert(file_linked_list[name]['Coordinates'][linked_list_ids_found[i], 2] == file_snapshot[name]['Coordinates'][ids_found[i], 2])

                                    assert(file_linked_list[name]['Velocities'][linked_list_ids_found[i], 0] == file_snapshot[name]['Velocities'][ids_found[i], 0])
                                    assert(file_linked_list[name]['Velocities'][linked_list_ids_found[i], 1] == file_snapshot[name]['Velocities'][ids_found[i], 1])
                                    assert(file_linked_list[name]['Velocities'][linked_list_ids_found[i], 2] == file_snapshot[name]['Velocities'][ids_found[i], 2])
 
    print("All particles have been accounted for.  Good job!!!")

def check_subfind_results(snapshot_idx):

    '''
    This function will be used to check if the results from SUBFIND are sensible/correct.  We will do a number of checks here.

    Since the 'snapshots' that we pass into SUBFIND contain only the bound particles, we can ensure that SUBFIND correctly put them all into groups.
    1: The 'halos/subfind_###.catalog_groups' should contain the number of groups given by the FoF tab.
    2: The 'halos/subfind_###.catalog_particles' should contain the number of particles given by the FoF tab.
    3: The 'catalogs/subfind_###.catalog_groups_properties' contains the number of particles within each FoF halo.  These particle numbers should match those given by the FoF tab.
    '''

    ## First let's load in all the FoF tab properties. ##

    '''
    fname_fof_tab_ids = "{1}groups_{0:03d}/fof_tab_{0:03d}.0.hdf5".format(snapshot_idx, groupdir)
    
    ds = yt.load(fname_fof_tab_ids) # Loads in the entire FoF catalogue for this snapshot.
    ad = ds.all_data() # Container for all the data that we can work with.

    # First load all the FoF particle IDs from Britton's FoF_tab files. #
    if len(ds.field_list) == 0: 
        num_groups = 0 
    else:
        num_groups = ad['Group'].shape[0] # Number of groups within this snapshot.
        num_particles_per_group = np.zeros((num_groups), dtype = np.int32)
        print("For snapshot {0} there are {1} groups.".format(snapshot_idx, num_groups))

        fof_tab_ids = [] 
        
        for group_idx in trange(num_groups): # Loop over each group.
            halo = ds.halo('Group', group_idx) # Loads all the information for the specified group.
            num_particles_per_group[group_idx] = halo['member_ids'].shape[0]

            for particle_idx in range(0, int(num_particles_per_group[group_idx])): # Loop over the number of particles within this group.
                parttype = int(halo['MemberType'][particle_idx])
                fof_tab_ids.append(np.int64(halo['member_ids'][particle_idx]))


    ## Check number of groups found by SUBFIND matches the FoF tab value. ##

    fname_subfind_groups = "{0}halos/subfind_{1:03d}.catalog_groups".format(subfind_dir, snapshot_idx)
   
    print("Reading from file {0}".format(fname_subfind_groups)) 
    with open(fname_subfind_groups, 'rb') as file_subfind_groups:
        N_groups = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1) 

    print("Number of groups in the FoF Tab is {0}.  Number of groups from SUBFIND is {1}".format(num_groups, N_groups))
#    assert(N_groups == num_groups)
    '''

    ## Do (Sub)Halo Mass Function ## 
    Halo_Desc_full = [
    ('id_MBP',              np.int64),
    ('M_vir',               np.float64),
    ('n_particles',         np.int16),
    ('position_COM',        (np.float32, 3)),
    ('position_MBP',        (np.float32, 3)),
    ('velocity_COM',        (np.float32, 3)),
    ('velocity_MBP',        (np.float32, 3)),
    ('R_vir',               np.float32),
    ('R_halo',              np.float32),
    ('R_max',               np.float32),
    ('V_max',               np.float32),
    ('sigma_v',             np.float32),
    ('spin',                (np.float32, 3)),
    ('q_triaxial',          np.float32),
    ('s_triaxial',          np.float32),
    ('shape_eigen_vectors', (np.float32, (3,3))),
    ('padding',             (np.int16, 2))
                     ] # Note that there are also a padding of 8 bytes following this array. 

    names = [Halo_Desc_full[i][0] for i in range(len(Halo_Desc_full))]
    formats = [Halo_Desc_full[i][1] for i in range(len(Halo_Desc_full))]
    Halo_Desc = np.dtype({'names':names, 'formats':formats}, align=True)

    # Reading SUBFIND Output #
    # Halos #
    fname_subfind_groups = "{0}catalogs/subfind_{1:03d}.catalog_groups_properties/subfind_{1:03d}.catalog_groups_properties.0".format(subfind_dir, snapshot_idx)

    print("Reading from file {0}".format(fname_subfind_groups)) 
    with open(fname_subfind_groups, 'rb') as file_subfind_groups:
        file_number = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)
        n_files = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)
        N_groups_thisfile = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)[0]
        N_groups_allfile = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)       
        
        Halos = np.fromfile(file_subfind_groups, Halo_Desc, N_groups_thisfile)  # Read in the galaxy structures
        
    halo_mass_subfind = np.log10(Halos['M_vir'])
    npart_subfind = Halos['n_particles']
    (halo_counts_subfind, bin_edges, bin_middle) = AllVars.Calculate_Histogram(halo_mass_subfind, bin_width, 0, 6, 11)
    (npart_binned_subfind_mean, npart_binned_subfind_std, npart_binned_subfind_N, bin_middle) = AllVars.Calculate_2D_Mean(halo_mass_subfind, npart_subfind, bin_width, 6, 11) 

    # Subhalos #
    fname_subfind_groups = "{0}catalogs/subfind_{1:03d}.catalog_subgroups_properties/subfind_{1:03d}.catalog_subgroups_properties.0".format(subfind_dir, snapshot_idx)

    print("Reading from file {0}".format(fname_subfind_groups)) 
    with open(fname_subfind_groups, 'rb') as file_subfind_groups:
        file_number = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)
        n_files = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)
        N_groups_thisfile = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)[0]
        N_groups_allfile = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)
       
        Subhalos = np.fromfile(file_subfind_groups, Halo_Desc, N_groups_thisfile)  # Read in the galaxy structures
        
    subhalo_mass_subfind = np.log10(Subhalos['M_vir'])
    (subhalo_counts_subfind, bin_edges, bin_middle) = AllVars.Calculate_Histogram(subhalo_mass_subfind, bin_width, 0, 6, 11)

    # Reading the FoF Tab #

    mass_fof_tab = []
    npart_fof_tab = []
    for core_idx in range(num_cores):
        fname_fof_tab_ids = "{1}groups_{0:03d}/fof_tab_{0:03d}.{2}.hdf5".format(snapshot_idx, groupdir, core_idx)

        with h5py.File(fname_fof_tab_ids, "r") as file_fof_tab:
            for group_idx in range(file_fof_tab['Group']['GroupMass'].shape[0]):
                mass_fof_tab.append(np.log10(file_fof_tab['Group']['GroupMass'][group_idx] * 1.0e10))
                npart_fof_tab.append(file_fof_tab['Group']['GroupLen'][group_idx])

    print("The FoF Tab had {0} groups whereas SUBFIND has {1} groups.".format(len(mass_fof_tab), len(halo_mass_subfind))) 
    (counts_fof_tab, bin_edges, bin_middle) = AllVars.Calculate_Histogram(mass_fof_tab, bin_width, 0, 6, 11) 
    (npart_binned_fof_tab_mean, npart_binned_fof_tab_std, npart_binned_fof_tab_N, bin_middle) = AllVars.Calculate_2D_Mean(mass_fof_tab, npart_fof_tab, bin_width, 6, 11) 

    ## Plotting ## 
    # Halo Mass Function #       
    ax1 = plt.subplot(211)

    label = "Subfind" 
    ax1.plot(bin_middle, halo_counts_subfind, label = label) 

    label = "FoF Tab"
    ax1.plot(bin_middle, counts_fof_tab, label = label)

    ax1.set_ylabel(r'Count', fontsize = PlotScripts.global_fontsize) 

    ax1.set_yscale('log', nonposy='clip')
    ax1.set_xlim([7, 11])
    ax1.set_ylim([0, max(counts_fof_tab) + 1000]) 

    leg = ax1.legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(PlotScripts.global_legendsize)

    # Subhalo Mass Function #

    ax2 = plt.subplot(212)

    ax2.plot(bin_middle, np.divide(npart_binned_subfind_mean, npart_binned_fof_tab_mean), 'r')

    ax2.set_xlim([7, 11])


    ax2.set_xlabel(r'$\log_{10}\ M_{H} \:[M_{\odot}]$', fontsize = PlotScripts.global_fontsize)
    ax2.set_ylabel(r'$N_\mathrm{Subfind} / N_\mathrm{FoF Tab}$', fontsize = PlotScripts.global_fontsize) 

    plt.tight_layout()

    outputFile = "./HaloCounts_z{0:.3f}.png".format(AllVars.SnapZ[snapshot_idx])
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))
    plt.close()


    # Halos and Subhalos Spatial Plot #
    ax1 = plt.subplot(111)

    w_halos = np.where((Halos['position_COM'][:,2] > 775.0) & (Halos['position_COM'][:,2] < 776.0))[0]
    w_subhalos = np.where((Subhalos['position_COM'][:,2] > 775.0) & (Subhalos['position_COM'][:,2] < 776.0))[0]

    for halo_idx in w_halos:
        log_mass = np.log10(Halos[halo_idx]['M_vir'])
        print("Halo {0}: x = {1:.4f} y = {2:.4f} z = {3:.4f} M_vir = {4:.4e}".format(halo_idx, Halos[halo_idx]['position_COM'][0], Halos[halo_idx]['position_COM'][1], Halos[halo_idx]['position_COM'][2], Halos[halo_idx]['M_vir']))      
        circ = Circle((Halos[halo_idx]['position_COM'][0], Halos[halo_idx]['position_COM'][1]), pow(1.2, log_mass), fill = False)
        ax1.add_patch(circ)

    for halo_idx in w_subhalos:
        log_mass = np.log10(Subhalos[halo_idx]['M_vir'])
        print("Halo {0}: x = {1:.4f} y = {2:.4f} z = {3:.4f} M_vir = {4:.4e}".format(halo_idx, Subhalos[halo_idx]['position_COM'][0], Subhalos[halo_idx]['position_COM'][1], Subhalos[halo_idx]['position_COM'][2], Subhalos[halo_idx]['M_vir']))      
        circ = Circle((Subhalos[halo_idx]['position_COM'][0], Subhalos[halo_idx]['position_COM'][1]), pow(1.1, log_mass), fill = False, color = 'b')
        ax1.add_patch(circ)



    ax1.set_xlim([770, 830])
    ax1.set_ylim([770, 830])

    outputFile = "./Halos_Subhalos_z{0:.3f}.png".format(AllVars.SnapZ[snapshot_idx])
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))
    plt.close()

def find_mass():

    fname = "{0}/snapdir_000/snapdir_000/snapshot_000.0.hdf5".format(snapdir)
    mass_table = [ 0.0, 0.0001150956450801459112886354629878837841, 0.0009207651606411672903090837039030702726, 0.0073661212851293383224726696312245621812, 0.0589289702810347065797813570497964974493, 0.4714317619800567626953125000000000000000] # For some reason the Mass Table in the header does not contain the mass for PartType5.
    mass_table = [ 0.0, 0.0001150956450801459112886354629878837841, 0.0009207651606411672903090837039030702726, 0.0073661212851293383224726696312245621812, 0.0589289702810347065797813570497964974493, 0.0] # For some reason the Mass Table in the header does not contain the mass for PartType5.

    mass_tot = 0
    numpart_tot = 0
    with h5py.File(fname, "r") as f:
        for type_idx in range(TypeMax):
            numpart_thistype = f['Header'].attrs['NumPart_Total'][type_idx]
            numpart_thistype += f['Header'].attrs['NumPart_Total_HighWord'][type_idx] << 32 
            

            numpart_tot += numpart_thistype
            mass_tot += mass_table[type_idx]*numpart_thistype 

            print("For PartType {0} there is {1} particles.".format(type_idx, numpart_thistype))

    print("The total number of particles within the simulation is {0}".format(numpart_tot))
    print("The total mass within the simulation is {0}".format(mass_tot))
if __name__ == '__main__':

    if (len(sys.argv) != 3):
        print("Usage: python3 rewrite_files.py <snaplow> <snaphigh>")
        exit()

    AllVars.Set_Params_Britton()
    PlotScripts.Set_Params_Plot()
    #find_mass()
   
    snaplow = int(sys.argv[1])
    snaphigh = int(sys.argv[2])

    print("Snaplow = {0}, snaphigh = {1}".format(snaplow, snaphigh))
 
    for snapshot_idx in range(snaplow, snaphigh + 1):  
        
        #write_fof_header(snapshot_idx) 
        #write_fof_groups(snapshot_idx)
        #write_snapshot_header(snapshot_idx)

        #link_fof_snapshot_ids(snapshot_idx) 
        #check_linking_ids(snapshot_idx)

        #link_fof_snapshot_full(snapshot_idx, 0)
        check_linking_full(snapshot_idx, 0)
        
        #check_subfind_results(snapshot_idx) 
