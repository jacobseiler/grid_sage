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

import sys
sys.path.append('/home/jseiler/SAGE-stuff/output/')

import AllVars
import PlotScripts
import ReadScripts

output_format = ".png"
dilute = 10000
bin_width = 0.1

num_cores = 256
TypeMax = 6
snaplow = 20
snaphigh = 100

def plot_snapshot(filepath, snapshot_idx):

    ax1 = plt.subplot(111)

    tmp = "snapdir_%03d/snapshot_%03d.129.hdf5" %(snapshot_idx, snapshot_idx)
    fname = filepath + tmp

    with h5py.File(fname, 'r') as f:

        for i in range(TypeMax):
            print("plotting ParticleType {0}".format(i))
            tmp = "PartType%d" %(i)
            try:
                part = f[tmp]['Coordinates']
            except KeyError:
                pass
            else:
                if(len(part) > dilute): 
                    w = sample(np.arange(0, len(part)), dilute) 
                    w = np.sort(w)
                else:
                    w = np.arange(0, len(part))  
               
                print("The x coords range from [{0:.4f}, {1:.4f}]".format(min(part[:,0]), max(part[:,0])))
                print("The y coords range from [{0:.4f}, {1:.4f}]".format(min(part[:,1]), max(part[:,1])))
                print("The z coords range from [{0:.4f}, {1:.4f}]".format(min(part[:,2]), max(part[:,2])))
                
                ax1.scatter(part[w,0], part[w,1], marker = 'o', alpha = 0.5, color = PlotScripts.colors[i-1])


    ax1.set_xlim([760, 830])
    ax1.set_ylim([760, 830])

    ax1.set_xlabel("x [mpc]")
    ax1.set_ylabel("y [mpc]")

    outputFile = './AllPart_%d%s' %(snapshot_idx, output_format) 
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
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

    min_y = [1e10, 1e10, 1e10, 1e10]
    max_y = [-1e10, -1e10, -1e10, -1e10]

    min_z = [1e10, 1e10, 1e10, 1e10]
    max_z = [-1e10, -1e10, -1e10, -1e10] 

    for core_idx in range(100):

        tmp = "snapdir_%03d/snapshot_%03d.%d.hdf5" %(snapshot_idx, snapshot_idx, core_idx)
        fname = filepath + tmp
        print(fname)
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
                    if(local_x_min < min_x[part_idx - 1]):
                        min_x[part_idx - 1] = local_x_min
                        print("New x min = {0} for PartType {1}".format(local_x_min, part_idx))
                    if(local_x_max > max_x[part_idx -1]):
                        max_x[part_idx - 1] = local_x_max
                        print("New x max = {0} for PartType {1}".format(local_x_max, part_idx))

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
        #plot_snapshot(filepath, snapshot_idx)
        #plot_halos(filepath, snapshot_idx)
        #check_bounds(filepath, snapshot_idx)
        plot_hmf(filepath, snapshot_idx, hmf) 
