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

output_format = ".png"
dilute = 10000
bin_width = 0.1

def plot_snapshot(filepath, snapshot_idx):

    ax1 = plt.subplot(111)

    tmp = "snapdir_%03d/snapshot_%03d.129.hdf5" %(snapshot_idx, snapshot_idx)
    fname = filepath + tmp

    with h5py.File(fname, 'r') as f:

        for i in xrange(1, 5):
            tmp = "PartType%d" %(i)
            part = f[tmp]['Coordinates']

            if(len(part) > dilute): 
                w = sample(np.arange(0, len(part)), dilute) 
                w = np.sort(w)
            else:
                w = np.arange(0, len(part))  
           
            print "The x coords range from [%.4f, %.4f]" %(min(part[:,0]), max(part[:,0]))
            print "The y coords range from [%.4f, %.4f]" %(min(part[:,1]), max(part[:,1]))
            print "The z coords range from [%.4f, %.4f]" %(min(part[:,2]), max(part[:,2]))
            
            ax1.scatter(part[w,0], part[w,1], marker = 'o', alpha = 0.5, color = PlotScripts.colors[i-1])
            print "plotting ParticleType %d" %(i)

    ax1.set_xlim([760, 830])
    ax1.set_ylim([760, 830])

    ax1.set_xlabel("x [mpc]")
    ax1.set_ylabel("y [mpc]")

    outputFile = './AllPart_%d%s' %(snapshot_idx, output_format) 
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

def plot_halos(filepath, snapshot_idx):

    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)
    
    mass = [] 

    have_fof = 0
    for core_idx in xrange(0, 256):
        
        tmp = "groups_%03d/fof_tab_%03d.%d.hdf5" %(snapshot_idx, snapshot_idx, core_idx)
        fname = filepath + tmp

        with h5py.File(fname, 'r') as f:

            try:
                position = f['Group']['GroupPos'] 
            except KeyError:
                pass
            else:
                have_fof = 1
                print "Found %d groups for snapshot %d, Core %d" %(len(position), snapshot_idx, core_idx)

                ax1.scatter(position[:,0], position[:,1], marker = 'o', alpha = 0.5, color = 'r')
                for i in xrange(0, len(position)):
                    mass.append(np.log10(f['Group']['GroupMass'][i] * 1.0e10))

    if (have_fof == 1):
        (counts, bin_edges, bin_middle) = AllVars.Calculate_Histogram(mass, bin_width, 0)
        
        # Box is [770, 825]
        ax2.plot(bin_middle,  counts / pow(825.0-770.0,3) / bin_width / pow(10,bin_middle), ls = '-', color = 'r', linewidth = PlotScripts.global_linewidth)
    ax2.set_yscale('log', nonposy='clip')
#    ax2.set_ylim([1e-4, 2e-1])
    ax2.set_xlim([6.0, 11.0])
    ax2.set_xlabel(r"$\log_{10}\mathrm{M}_\mathrm{H} \: [\mathrm{M}_\odot]$", fontsize = PlotScripts.global_fontsize)
    ax2.set_ylabel(r'$\left(\frac{dn}{dM}\right) [\mathrm{Mpc}^{-3}\: \mathrm{M}_\odot^{-1}]$', fontsize = PlotScripts.global_fontsize)

    ax1.set_xlim([775, 830])
    ax1.set_ylim([765, 830])

    ax1.set_xlabel(r"$\mathrm{x} \: [\mathrm{Mpc}]$")
    ax1.set_ylabel(r"$\mathrm{y} \: [\mathrm{Mpc}]$")

    plt.tight_layout()
    outputFile = './halo_1721228/%d%s' %(snapshot_idx, output_format) 
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

if __name__ == '__main__':

    PlotScripts.Set_Params_Plot()
    filepath = '/lustre/projects/p004_swin/bsmith/1.6Gpc/means/halo_1721228/dm_gadget/data/'

    for snapshot_idx in xrange(0, 131):
        #plot_snapshot(filepath, snapshot_idx)
        plot_halos(filepath, snapshot_idx) 
