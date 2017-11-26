#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')

import astropy

import h5py
import numpy as np
import pylab as plt
from hmf import MassFunction
from hmf import cosmo
from astropy.cosmology import FlatLambdaCDM
from random import sample
import matplotlib.patches as patches

import sys
sys.path.append('/home/jseiler/SAGE-stuff/output/')

import AllVars
import PlotScripts
import ReadScripts

output_format = ".png"
dilute = 500000
bin_width = 0.1

num_cores = 1
TypeMax = 6
snaplow = 0
snaphigh = 0

groupdir = '/lustre/projects/p004_swin/bsmith/1.6Gpc/means/halo_1721673/dm_gadget/data/'
snapdir = '/lustre/projects/p134_swin/jseiler/simulations/1.6Gpc/means/halo_1721673/dm_gadget/data/'
linking_outdir = '/lustre/projects/p134_swin/jseiler/simulations/my_britton_fof_full/'


def plot_snapshot(filepath, snapshot_idx, do_slice, do_flip):

    fig = plt.figure()

    #ax1 = fig.add_subplot(121) 
    #ax2 = plt.axes([0.58, 0.58, 0.29, 0.29])
    ax2 = fig.add_subplot(111) 

    for core_idx in range(0, num_cores):
        print("Plotting chunk {0}".format(core_idx))
        if snapshot_idx < 66:
            fname = "{2}snapdir_{0:03d}/snapdir_{0:03d}/snapshot_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx, snapdir)
        else:
            fname = "{2}snapdir_{0:03d}/snapshot_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx, snapdir)
        fname = "/lustre/projects/p004_swin/bsmith/1.6Gpc/means/halo_1721673/dm_gadget/data/snapdir_020/snapshot_020.0.hdf5"

        print("Reading from file {0}".format(fname)) 
        with h5py.File(fname, 'r') as f:
            for type_idx in range(TypeMax): 
                print("plotting ParticleType {0}".format(type_idx))
                tmp = "PartType{0}".format(type_idx)
                try:
                    part = f[tmp]['Coordinates'][:]
                except KeyError:
                    pass
                else:


                    if (do_slice == 1):
                        #w = np.where((part[:,2] > 800.0) & (part[:,2] < 804.0))[0]
                        w = np.where((part[:,0] > 817.0) & (part[:,1] < 779.0))[0]
                        print("There is {0} particles within this slice.".format(len(w)))
                    else:
                        if(len(part) > dilute):    
                            print("There are {0} particles and we are diluting it down to {1}".format(len(part), dilute)) 
                            w = sample(list(np.arange(0, len(part))), dilute) 
    #                        w = np.sort(w)
                        else:
                            w = np.arange(0, len(part))  

#                    print("The x coords range from [{0:.4f}, {1:.4f}]".format(min(part[:,0]), max(part[:,0])))
#                    print("The y coords range from [{0:.4f}, {1:.4f}]".format(min(part[:,1]), max(part[:,1])))
#                    print("The z coords range from [{0:.4f}, {1:.4f}]".format(min(part[:,2]), max(part[:,2])))
   
#                    ax1.scatter(part[w,0], part[w,1], marker = 'o', alpha = 0.5, color = PlotScripts.colors[type_idx-1], s = 1)
                    
                    if(do_flip == 1):
                        #ax2.scatter(part[w,1], part[w,0], marker = 'o', alpha = 0.5, color = PlotScripts.colors[type_idx-1], s = 1)                
                        ax2.scatter((part[w,1] - 775.0), (part[w,0] - 775.0), marker = 'o', alpha = 0.5, color = PlotScripts.colors[type_idx-1], s = 1)                
                    else:
                        #ax2.scatter(part[w,0] - 775.0, part[w,1] - 775.0, marker = 'o', alpha = 0.5, color = PlotScripts.colors[type_idx-1], s = 1)                
                        ax2.scatter(part[w,2] - 775.0, part[w,0] - 775.0, marker = 'o', alpha = 0.5, color = PlotScripts.colors[type_idx-1], s = 1)                
#                        ax2.scatter((part[w,0] - 775.0), (part[w,1] - 775.0), marker = 'o', alpha = 0.5, color = PlotScripts.colors[type_idx-1], s = 1)                

    #ax2.add_patch(patches.Rectangle((817, 779), 1.0, 1.0, fill=False))
    #ax2.add_patch(patches.Rectangle((45.0, 45.0), 5.0, 5.0, fill=False))
    #ax2.add_patch(patches.Rectangle((10.0, 45.0), 5.0, 5.0, fill=False))

    ax2.add_patch(patches.Rectangle((45.0, 45.0), 5.0, 5.0, fill=False))
    ax2.add_patch(patches.Rectangle((45.0, 0.0), 5.0, 10.0, fill=False))


    for type_idx in range(1, TypeMax): # PartType0 has no particles ever.
        label = 'PartType{0}'.format(type_idx)
        ax2.scatter(np.nan, np.nan, color = PlotScripts.colors[type_idx-1], label = label)

    #ax1.set_xlim([0, 1600])
    #ax1.set_ylim([0, 1600])

    #ax2.set_xlim([774, 827])
    #ax2.set_ylim([774, 827])

    ax2.set_xlim([0, 50])
    ax2.set_ylim([0, 50])


    #ax1.set_xlabel("")
    #ax1.set_ylabel("y [Mpc/h]")

    ax2.set_xlabel("z [Mpc/h]")    
    ax2.set_ylabel("x [Mpc/h]")
    
    leg = ax2.legend(loc='upper left', numpoints=1,labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(PlotScripts.global_legendsize)
        t.set_alpha(1)

    plt.tight_layout()
    
   
    if do_flip == 1:
        flip_tag = "flip"
    else:
        flip_tag = "noflip"

    if do_slice == 1:
        outputFile = './Slice_Core{2}_{0}{3}{1}'.format(snapshot_idx, output_format, core_idx, flip_tag) 
    else:
        outputFile = './Grid_AllPart_Core{2}_{0}{3}{1}'.format(snapshot_idx, output_format, core_idx, flip_tag)
     
    outputFile = './Grid_AllPart_Snapshot20_chunk0_zx.png'
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))
    plt.close()

def plot_density_grid(filepath, snapshot_idx, GridSize):

    fig = plt.figure()

    ax1 = fig.add_subplot(111) 

    density_grid = np.zeros((GridSize, GridSize, GridSize), np.float64)

    bound_low = 775.0
    bound_high = 825.0
    print("Creating a density grid.")
    for core_idx in range(10, num_cores):
        print("Summing chunk {0}".format(core_idx))
        if snapshot_idx < 66:
            fname = "{2}snapdir_{0:03d}/snapdir_{0:03d}/snapshot_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx, snapdir)
        else:
            fname = "{2}snapdir_{0:03d}/snapshot_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx, snapdir)
    
        with h5py.File(fname, 'r') as f:
            for type_idx in range(TypeMax): 
                print("plotting ParticleType {0}".format(type_idx))
                tmp = "PartType{0}".format(type_idx)
                try:
                    part = f[tmp]['Coordinates'][:]
                except KeyError:
                    pass
                else:
                    for part_idx in range(len(part)):
                        if(part_idx % 1e5 == 0):
                            print("Particle {0}".format(part_idx))
                        x_pos = part[part_idx,0]
                        y_pos = part[part_idx,1]
                        z_pos = part[part_idx,2]

                        if (x_pos < bound_low or x_pos > bound_high or y_pos < bound_low or y_pos > bound_high or z_pos < bound_low or z_pos > bound_high):
                            continue

                        x_grid = int((x_pos - bound_low) * GridSize/ AllVars.BoxSize) 
                        y_grid = int((y_pos - bound_low) * GridSize/ AllVars.BoxSize) 
                        z_grid = int((z_pos - bound_low) * GridSize/ AllVars.BoxSize) 

                        if(x_grid >= GridSize):
                            x_grid = int(GridSize -1)
                        if(x_grid < 0):
                            x_grid = 0

                        if(y_grid >= GridSize):
                            y_grid = int(GridSize -1)
                        if(y_grid < 0):
                            y_grid = 0

                        if(z_grid >= GridSize):
                            z_grid = int(GridSize -1)
                        if(z_grid < 0):
                            z_grid = 0

                        if(type_idx != 5):  
                            mass = f['Header'].attrs['MassTable'][type_idx]
                        else:
                            mass = 0.471431761 # For some reason Britton's Sim doesn't have a proper mass for PartType5 in the MassTable.

                        density_grid[x_grid, y_grid, z_grid] += mass
 
    total_mass = np.sum(density_grid)
    norm_density = total_mass / pow(GridSize, 3)
    print("The total mass is {0} and the norm_density is {1}".format(total_mass, norm_density))
    
    density_grid /= norm_density

    cut_slice = 0
    im = ax1.imshow(density_grid[:,:,cut_slice:cut_slice+127].mean(axis = -1), interpolation='bilinear', origin='low', extent =[0,AllVars.BoxSize,0,AllVars.BoxSize], cmap = 'Purples', vmin = 0.12, vmax = 25.0)

    cbar = plt.colorbar(im, ax = ax1)
    cbar.set_label(r'$\rho/\langle \rho \rangle$')

    ax1.set_xlabel(r'$\mathrm{x}  (h^{-1}Mpc)$')
    ax1.set_ylabel(r'$\mathrm{y}  (h^{-1}Mpc)$')

    ax1.set_xlim([0.0, AllVars.BoxSize])
    ax1.set_ylim([0.0, AllVars.BoxSize])
    
    outputFile = './density_grid_snapshot{0}_core{1}'.format(snapshot_idx, num_cores - 1) 
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))

    plt.close()




def plot_halos(filepath, snapshot_idx):

    #ax1 = plt.subplot(211)
    ax2 = plt.subplot(111)
    
    mass = [] 

    have_fof = 1
    for core_idx in range(num_cores):
        
        tmp = "groups_%03d/fof_tab_%03d.%d.hdf5" %(snapshot_idx, snapshot_idx, core_idx)
        fname = filepath + tmp

        with h5py.File(fname, 'r') as f:

            try:
                position = f['Group']['GroupPos'] 
            except KeyError:
                pass
            else:
                have_fof = 1
                print("Found {0} groups for snapshot {1}, Core {2}".format(len(position), snapshot_idx, core_idx))

                #ax1.scatter(position[:,0], position[:,1], marker = 'o', alpha = 0.5, color = 'r')
                for i in range(len(position)):
                    mass.append(np.log10(f['Group']['GroupMass'][i] * 1.0e10))

    if (have_fof == 1):
        (counts, bin_edges, bin_middle) = AllVars.Calculate_Histogram(mass, bin_width, 0)
        
        # Box is [774, 831]
        ax2.plot(bin_middle,  counts / pow(831.0-774.0,3) / bin_width / pow(10,bin_middle), ls = '-', color = 'r', linewidth = PlotScripts.global_linewidth)
    ax2.set_yscale('log', nonposy='clip')
#    ax2.set_ylim([1e-4, 2e-1])
    ax2.set_xlim([6.0, 11.0])
    ax2.set_xlabel(r"$\log_{10}\mathrm{M}_\mathrm{H} \: [\mathrm{M}_\odot]$", fontsize = PlotScripts.global_fontsize)
    ax2.set_ylabel(r'$\left(\frac{dn}{dM}\right) [\mathrm{Mpc}^{-3}\: \mathrm{M}_\odot^{-1}]$', fontsize = PlotScripts.global_fontsize)

    #ax1.set_xlim([775, 830])
    #ax1.set_ylim([765, 830])

    #ax1.set_xlabel(r"$\mathrm{x} \: [\mathrm{Mpc}]$")
    #ax1.set_ylabel(r"$\mathrm{y} \: [\mathrm{Mpc}]$")

    plt.tight_layout()
    outputFile = './halo_1721228/%d%s' %(snapshot_idx, output_format) 
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))
    plt.close()

def check_bounds(filepath, snapshot_idx):

    
    min_x = np.full(TypeMax, 1e10)
    max_x = np.full(TypeMax, -1e10)

    min_y = np.full(TypeMax, 1e10)
    max_y = np.full(TypeMax, -1e10)

    min_z = np.full(TypeMax, 1e10)
    max_z = np.full(TypeMax, -1e10)

    for core_idx in range(num_cores):

        tmp = "snapdir_%03d/snapshot_%03d.%d.hdf5" %(snapshot_idx, snapshot_idx, core_idx)
        fname = filepath + tmp
        print("Chunk {0}".format(core_idx))
        with h5py.File(fname, 'r') as f:

            for part_idx in range(TypeMax):
                tmp = "PartType%d" %(part_idx)
            
                try:    
                    part = f[tmp]['Coordinates']
                except KeyError:
                    pass
                else:

                    local_x_min = min(part[:,0])
                    local_x_max = max(part[:,0])
                    if(local_x_min < min_x[part_idx]):
                        min_x[part_idx ] = local_x_min
                        print("New x min = {0} for PartType {1}".format(local_x_min, part_idx))
                    if(local_x_max > max_x[part_idx]):
                        max_x[part_idx] = local_x_max
                        print("New x max = {0} for PartType {1}".format(local_x_max, part_idx))

                    local_y_min = min(part[:,1])
                    local_y_max = max(part[:,1])
                    if(local_y_min < min_y[part_idx]):
                        min_y[part_idx] = local_y_min
                        print("New y min = {0} for PartType {1}".format(local_y_min, part_idx))
                    if(local_y_max > max_y[part_idx]):
                        max_y[part_idx] = local_y_max
                        print("New y max = {0} for PartType {1}".format(local_y_max, part_idx))

                    local_z_min = min(part[:,2])
                    local_z_max = max(part[:,2])
                    if(local_z_min < min_z[part_idx]):
                        min_z[part_idx] = local_z_min
                        print("New z min = {0} for PartTzpe {1}".format(local_z_min, part_idx))
                    if(local_z_max > max_z[part_idx]):
                        max_z[part_idx] = local_z_max
                        print("New z max = {0} for PartTzpe {1}".format(local_z_max, part_idx))


    print("For snapshot {0} the smallest x coorindate is {1} [Mpc/h] and maximum is {2} [Mpc/h] (for each particle type)".format(snapshot_idx, min_x, max_x)) 

def plot_hmf(filepath, snapshot_idx, hmf):
        
    tmp = "snapdir_{0:03d}/snapshot_{0:03d}.0.hdf5".format(snapshot_idx)
    fname = filepath + tmp

    mass = []
    pos_x = []
    pos_y = []
    pos_z = []
    have_fof_this_snap = 0

    cosmol = AllVars.Set_Params_Britton() # Set the parameters for Britton's model.
    with h5py.File(fname, 'r') as f:       
        z = f['Header'].attrs['Redshift']
        hmf.update(z = z)
        print("Doing snapshot {1} at redshift {0}".format(z, snapshot_idx))  
        
        massfunc = hmf.dndlog10m
        hmf_bins = np.linspace(6.0, 11.0, num = (11.0 - 6.0) / 0.01)
   
#        print(massfunc * pow(AllVars.Hubble_h, 3)  

    for core_idx in range(num_cores):
        if core_idx % 100 == 0:
            print("Doing core {0}".format(core_idx))
        tmp = "groups_{0:03d}/fof_tab_{0:03d}.{1}.hdf5".format(snapshot_idx, core_idx)
        fname = filepath + tmp
        
        with h5py.File(fname, 'r') as f:       
            Ngroups_ThisFile = f['Header'].attrs['Ngroups_ThisFile']
            Ngroups_Total = f['Header'].attrs['Ngroups_Total']
            if(Ngroups_ThisFile > 0):               
                have_fof_this_snap = 1 
                for group_idx in range(Ngroups_ThisFile):
                    mass.append(np.log10(f['Group']['GroupMass'][group_idx] * 1.0e10 / AllVars.Hubble_h))
#                    pos_x.append(f['Group']['GroupPos'][group_idx][0]) 
#                    pos_y.append(f['Group']['GroupPos'][group_idx][1]) 
#                    pos_z.append(f['Group']['GroupPos'][group_idx][2]) 

    assert(len(mass) == Ngroups_Total)

#    print("Minimum x position is {0}, maximum x position is {1}".format(min(pos_x), max(pos_x)))
#    print("Minimum y position is {0}, maximum y position is {1}".format(min(pos_y), max(pos_y)))
#    print("Minimum z position is {0}, maximum z position is {1}".format(min(pos_z), max(pos_z)))

    if (have_fof_this_snap == 1):
        (counts, bin_edges, bin_middle) = AllVars.Calculate_Histogram(mass, bin_width, 0, 6, 11)
     
    ax1 = plt.subplot(111)

    label = "Britton : z = {0:.2f}".format(z)
    ax1.plot(bin_middle, np.multiply(counts / pow(AllVars.BoxSize,3) * pow(AllVars.Hubble_h, 3) / bin_width, bin_middle), label = label) 
    
    label = "HMF : z = {0:.2f}".format(z)
    ax1.plot(hmf_bins - np.log10(AllVars.Hubble_h), massfunc * pow(AllVars.Hubble_h,3), label = label)

    ax1.set_xlabel(r'$\log_{10}\ M_{H} \:[M_{\odot}]$', fontsize = PlotScripts.global_fontsize)
    ax1.set_ylabel(r'$\left(\frac{dn}{d\log{M}}\right) \ [\mathrm{Mpc}^{-3}]$', fontsize = PlotScripts.global_fontsize)

    plt.yscale('log', nonposy='clip')
    plt.axis([6, 11.5, 1e-6, 5e0])

    leg = plt.legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(PlotScripts.global_legendsize)

    outputFile = "hmf_z{0:.2f}.png".format(z)
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))
    plt.close()
 
if __name__ == '__main__':

    PlotScripts.Set_Params_Plot()
    filepath = '/lustre/projects/p004_swin/bsmith/1.6Gpc/means/halo_1721673/dm_gadget/data/'   

    cosmol = AllVars.Set_Params_Britton() # Set the parameters for Britton's model.
    my_cosmo = cosmo.Cosmology(cosmo_model = cosmol) # Update the hmf cosmology.   
    britton_cosmo = FlatLambdaCDM(H0 = 69.5, Om0 = 0.285, Ob0 = 0.04845) 
    hmf = MassFunction()
    hmf.update(cosmo_params = {"H0" : 69.5, "Om0" : 0.285}, Mmax = 11, Mmin = 6)
 
    for snapshot_idx in range(snaplow, snaphigh + 1):
        #plot_snapshot(filepath, snapshot_idx, 0, 0)
        #plot_density_grid(filepath, snapshot_idx, 128)
        #plot_halos(filepath, snapshot_idx)
        #check_bounds(filepath, snapshot_idx)
        #plot_hmf(filepath, snapshot_idx, hmf)         
