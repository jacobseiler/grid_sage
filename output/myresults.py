import matplotlib
matplotlib.use('Agg')

import os
import heapq
import h5py as h5
import numpy as np
import pylab as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from numpy import *
from random import sample, seed, randint
from os.path import getsize as getFileSize
import math
import random
import csv
from cycler import cycler
from io import StringIO
#np.set_printoptions(threshold=np.nan)
from collections import Counter
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import AxesGrid
from astropy import units as u
from astropy import cosmology


import matplotlib.ticker as mtick
import PlotScripts
import ReadScripts
import AllVars

matplotlib.rcdefaults()
plt.rc('axes', color_cycle=['k', 'b', 'r', 'g', 'm', 'c','y', '0.5'], labelsize='x-large')
plt.rc('lines', linewidth='2.0')
# plt.rc('font', variant='monospace')
plt.rc('legend', numpoints=1, fontsize='x-large')
plt.rc('text', usetex=True)

label_size = 20
extra_size = 2
talk_fontsize = 20
talk_legendsize = 18
plt.rc('xtick', labelsize=talk_fontsize)
plt.rc('ytick', labelsize=talk_fontsize)
plt.rc('text', usetex=True)
tick_interval = 0.25
np.set_printoptions(formatter={'float': lambda x: "{0:0.10e}".format(x)})

colors = ['r', 'b', 'g', 'c', 'm', 'k']
markers = ['x', 'o', '^', 's', 'D']
linestyles = ['-', '--', '-.', ':']

AllVars.Set_Constants()
AllVars.Set_Params_Mysim()
cosmo = cosmology.FlatLambdaCDM(H0 = AllVars.Hubble_h*100, Om0 = AllVars.Omega_m) 
t_BigBang = cosmo.lookback_time(100000).value # Lookback time to the Big Bang in Gyr.
z_plot = np.arange(6, 14)  #Range of redshift we wish to plot.
time_xlim = [315, 930]
time_tick_interval = 25

Output_Format = ".png"

def calculate_beta(MUV, z):
	
## Calculation of the dust attenuation parameter Beta ##
## Fit values are from Bouwens (2015) ApJ 793, 115 ##
## For z = 5 and 6, Bouwens uses a piece-wise linear relationship and a linear relationship for higher redshift. ##

## INPUT ##
# MUV: A value of the absolute magnitude in the UV (generally M1600) in the AB magnitude system.
# z: Redshift value.

## OUTPUT ##
# beta: Value of the UV continuum paramaeter beta. 

	if (z >= 4.5 and z < 5.5): # z = 5 fits.
		if (MUV > -18.8):
			dB = -0.08
		else:
			dB = -0.17
		B = -2.05
		offset = 18.8
	elif (z >= 5.5 and z < 6.5): # z = 6 fits.
		if (MUV > -18.8):
			dB = -0.08
		else:
			dB = -0.24
		B = -2.22
		offset = 18.8

	elif (z >= 6.5 and z < 7.5): # z = 7 fits.
		dB = -0.20
		B = -2.05
		offset = 19.5
	elif (z >= 7.5 and z < 8.5): # z = 8 fits.
		dB = -0.15
		B = -2.13
		offset = 19.5
	elif (z >= 8.5 and z < 9.5): # z = 9 fits.
		dB = -0.16
		B = -2.19
		offset = 19.5
	elif (z >= 9.5 and z < 10.5): # z = 10 fits.
		dB = -0.16
		B = -2.16
		offset = 19.5

	beta = dB * (MUV + offset) + B

	return beta
		
def multiply(n):
	total = 1
	for i in xrange(0, len(n)):
		total *= n[i]
		print "%.4e" %(total)
	return total

def Calculate_Histogram(Data, Bin_Width, Weights, min_hist=None, max_hist=None):

# This calculates the counts and Bin_Edges for a given set of data.

## Input ##
# Data is an array containing the data to be binned.
# Bin_Width is the width of the bins.
# Weights is either 0 or 1.  0: Denotes that the histogram should be a frequency (count) histogram. 1: Denotes that the histogram should be a probability histogram.
# min_hist: Minimum value that will be binned.  OPTIONAL: If not specified will be given by the minimum of the data - 10 times the binwidth.
# max_hist: Minimum value that will be binned.  OPTIONAL: If not specified will be given by the maximum of the data - 10 times the binwidth.


## Output ##
# Counts: The count (either frequency or probability) in each bin.
# Bin_Edges: The location of the edges of the bins.
# Bin_Middle: The middle of the bins.

    if (min_hist == None): 
	mi = np.floor(min(Data)) - 10*Bin_Width
	ma = np.floor(max(Data)) + 10*Bin_Width
    else:
	print "EQHTQEIT"
	mi = min_hist 
	ma = max_hist 

    print "The minimum of the data being binned is %.4e" %(mi)
    print "The maximum of the data being binned is %.4e" %(ma) 
	    
    NB = (ma - mi) / Bin_Width 

#    print "The total of the data being binned is %.4e" %(sum(Data[Data >= mi and Data <= ma]))

    if (Weights == 0):
        (counts, Bin_Edges) = np.histogram(Data, range=(mi, ma), bins=NB)
    else:
        weights = np.ones_like(Data)/len(Data)
        (counts, Bin_Edges) = np.histogram(Data, range=(mi, ma), bins=NB, weights = weights)

    Bin_Middle = Bin_Edges[:-1] + 0.5 * Bin_Width

    return (counts, Bin_Edges, Bin_Middle)

def Sum_Log(Array):

    Sum = 0
    for i in xrange(0, len(Array)):
    	Sum += 10**Array[i]

    return Sum


def Std_Log(Array, Mean):

    Sum = 0
    for i in xrange(0, len(Array)):
	Sum += (10**Array[i] - Mean)**2

    Sum *= 1.0/len(Array)

    Std = np.sqrt(Sum)
    return Std

def Plot_Spectrum():

    SFR = [0.1, 0.3, 0.9, 2.7, 8.1, 24.3]
    Lambda_tot = []
    Total_tot = []
    Stellar_tot = []
    Labels = []
    for i in xrange(0, len(SFR)):

    	fname = "/Users/100921091/Desktop/sage-forked/sage-1/STARBURST/output/output_files/SFR_%.1f_Z_33.spectrum1" %(SFR[i])
	fd = open(fname, 'rb')
	t, Lambda, Total, Stellar, Nebular = np.loadtxt(fd, dtype = np.float32, skiprows=6, unpack=True)  
    
	Lambda_tot.append(Lambda[1222:2442])
	Total_tot.append(Total[1222:2442])
	Stellar_tot.append(Stellar[1222:2442])
	Labels.append("%.1f" %(SFR[i]))


    print Stellar_tot[0]
    print Stellar_tot[1]
    print Stellar_tot[1]/Stellar_tot[0]

    PlotScripts.Plot_Line(Lambda_tot, Stellar_tot, [0], [1,0, ["Galaxies"]], Labels, [80, 2e6, 28, 42], [1, 0], ['Lambda', 'L'], 1, 'Stellar2', '.png') 

def Plot_Ngamma_Ratio():

    SFR = np.arange(0.01, 1, 0.01)
    Ngamma = []
    delta_SFR = np.empty(len(SFR)-1)
    delta_Ngamma = np.empty(len(SFR)-1)

    for i in xrange(0, len(SFR)):
        fname = "/Users/100921091/Desktop/sage-forked/sage-1/STARBURST/output/output_files/SFR_%.2f_Z_35.quanta1" %(SFR[i])
        fd = open(fname, 'rb')
        t, N_HI, L_HI, N_HeI, L_HeI, N_HeII, L_HeII, L = np.loadtxt(fd, dtype = np.float32, skiprows=7, unpack=True)	

        Ngamma.append(N_HI[-1])

    for i in xrange(1, len(SFR)):
        delta_SFR[i-1] = SFR[i]/float(SFR[i-1])
        delta_Ngamma[i-1] = Ngamma[i]/float(Ngamma[i-1])

    SFR = np.log10(SFR)

    PlotScripts.Fit_Linear([SFR], [Ngamma], ["Galaxies"], [min(SFR), max(SFR), min(Ngamma), max(Ngamma)], [r'log(SFR)', r'log(N$\gamma$)'], 1, 'SFR_Ngamma_Z35', '.png') 
#    PlotScripts.Plot_Scatter([SFR], [Ngamma], ["Galaxies"], [min(SFR), max(SFR), min(Ngamma), max(Ngamma)], [1,0], ['Delta SFR', 'Delta Ngamma'], 1, 'deltaSFR_deltaNgamma', '.png') 


def Photon_Totals(Simulation, Redshift, Photons, Mysim_Len): 


    title = ["MySim 512 Halo", "MySim 512 Galaxy", "Millennium 1024 Halo", "Millennium 1024 Galaxy"]

    if Simulation == 0:
        OutputFile = 'TotalPhotons_Millennium'
    elif Simulation == 1:
        OutputFile = 'TotalPhotons_MySim'
    elif Simulation == 2:
        OutputFile = 'TotalPhotons_Both'

    print Redshift
    print Photons

    PlotScripts.Plot_Scatter(Redshift, Photons, title, [4.5, 15.5, 53, 58.3], [0,0], ['z', r'Log Ionizing Photons [s$^{-1}$]'], 1, 2, OutputFile, '.png')  

def StellarMassFunction(Simulation, Redshift, Mass, HaloPartStellarMass, MySim_Len, Model_Tags, Output_Tag):

    title = []
    Normalization = []

    ## Plot Parameters ##
    binwidth = 0.1
    Observations = 2

    Frequency = 0 # 0 for a frequency (count) histogram, 1 for a probability histogram.
    errorwidth = 2
    delta = 0.05
    caps = 5
    ##

    ## Normalization and Titles for MySim ##
    for i in xrange(0, MySim_Len): 

        #tmp = r'Recursive $\texttt{SAGE}, f_\mathrm{esc} \: \propto \: M_H^\beta$: z = %.2f' %(Redshift[i])
        tmp = r'z = %.2f' %(Redshift[i])

        title.append(tmp)

        #AllVars.Set_Params_Mysim()
      	AllVars.Set_Params_MiniMill()
      
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
    
        Normalization.append(norm)


    ## Normalization and Titles for MySim2 ## 
    
    for i in xrange(0, MySim_Len): 
        #tmp = r'Recursive $\texttt{SAGE}, f_\mathrm{esc} \: \propto \: M_H^{-\beta}$: z = %.2f' %(Redshift[i])
        tmp = r'Delayed: z = %.2f' %(Redshift[i])
        title.append(tmp)

        #AllVars.Set_Params_Mysim()
      	AllVars.Set_Params_MiniMill()
       
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
    
        Normalization.append(norm)
   
### Plotting ###

    plt.figure()  
    ax = plt.subplot(111)  

### Plots the Histograms ###

    for i in xrange(0, len(Mass)):

        (counts, Bin_Edges, Bin_Middle) = Calculate_Histogram(Mass[i], binwidth, Frequency)
        if (i < MySim_Len):
            ls = linestyles[0]
	    label = title[i]
        else:
            ls = linestyles[1]
	    label = ''
        plt.plot(Bin_Middle, counts / Normalization[i], color = colors[i % MySim_Len], linestyle = ls, label = label)

##

### Manually enter in the labels for the different simulations. ###
    for i in xrange(0, 2):
	    plt.plot(1e100, 1e100, color = 'k', ls = linestyles[i], label = Model_Tags[i])

### Draws a vertical line to denote lower bounds for what is an 'acceptable' Stellar Mass ### 

    plt.axvline(x = np.log10(HaloPartStellarMass), ymin = 0, ymax = 10, linestyle = '-.')

## 

    plt.yscale('log', nonposy='clip')

    plt.axis([8, 12.0, 1e-6, 1e-1])

    ax.set_xlabel(r'$\log_{10}\ m_{\mathrm{*}} \:[M_{\odot}]$', fontsize = talk_fontsize)
    ax.set_ylabel(r'$\Phi\ [\mathrm{Mpc}^{-3}\: \mathrm{dex}^{-1}]$', fontsize = talk_fontsize)
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

### If we want to put observations on the figure ###

    if Observations == 1:
    	Gonzalez_z6 = np.array([[7.77 + delta, -2.0956, -1.8596, -2.3539],
                                [8.27 + delta, -2.1742, -1.9494, -2.4101],
                                [8.77 + delta, -2.5674, -2.3876, -2.7921],
                                [9.27 + delta, -2.8483, -2.6573, -3.0843],
                                [9.77 + delta, -3.5787, -3.3764, -3.8258],
                                [10.27 + delta, -4.3202, -4.0281, -4.5674]], dtype = np.float32)

                                #plt.errorbar(Gonzalez_z6[:,0], 10**Gonzalez_z6[:,1], yerr= (10**Gonzalez_z6[:,3], 10**Gonzalez_z6[:,2]) , alpha=0.75, lw=1.0, marker='o', ls='none', label = 'Gonzalez 2011, z = 6', color = 'cyan')


        Gonzalez_z7 = np.array([[7.75, -2.1828, -1.7463, -2.5858],
                                [8.26, -2.25, -1.8694, -2.2631],
                                [8.77, -2.7425, -2.3731, -3.1231],
                                [9.27, -3.0672, -2.6753, -3.4142],
                                [9.76, -3.8731, -3.4831, -4.2537]], dtype = np.float32)

#plt.errorbar(Gonzalez_z7[:,0], 10**Gonzalez_z7[:,1], yerr= (10**Gonzalez_z7[:,3], 10**Gonzalez_z7[:,2]) , alpha=0.75, lw=1.0, marker='o', ls='none'    , label = 'Gonzalez 2011, z = 7', color = 'magenta')

        Song_z6 = np.array([[7.25 - delta, -1.47, -1.47 + 0.35, -1.47 - 0.23],
                            [7.75 - delta, -1.81, -1.81 + 0.23, -1.81 - 0.28],
                            [8.25 - delta, -2.26, -2.26 + 0.21, -2.26 - 0.16],
                            [8.75 - delta, -2.65, -2.65 + 0.15, -2.65 - 0.15],
                            [9.25 - delta, -3.14, -3.14 + 0.12, -3.14 - 0.11],
                            [9.75 - delta, -3.69, -3.69 + 0.12, -3.69 - 0.13],
                            [10.25 - delta, -4.27, -4.27 + 0.38, -4.27 - 0.86]], dtype = np.float32)

        plt.errorbar(Song_z6[:,0], 10**Song_z6[:,1], yerr= (10**Song_z6[:,1] - 10**Song_z6[:,3], 10**Song_z6[:,2] - 10**Song_z6[:,1]), xerr = 0.25, capsize = caps, elinewidth = errorwidth, alpha = 1.0, lw=2.0, marker='o', ls='none', label = 'Song 2015, z = 6', color = '#bd0026')

        Song_z7 = np.array([[7.25, -1.63, -1.63 + 0.54, -1.63 - 0.54],
                            [7.75, -2.07, -2.07 + 0.45, -2.07 - 0.41],
                            [8.25, -2.49, -2.49 + 0.38, -2.49 - 0.32],
                            [8.75, -2.96, -2.96 + 0.32, -2.96 - 0.30],
                            [9.25, -3.47, -3.47 + 0.32, -3.47 - 0.35],
                            [9.75, -4.11, -4.11 + 0.41, -4.11 - 0.57],
                            [10.25, -4.61, -4.61 + 0.72, -4.61 - 0.82],
                            [10.75, -5.24, -5.24 + 0.90, -5.25 - 0.57]], dtype = np.float32)

        plt.errorbar(Song_z7[:,0], 10**Song_z7[:,1], yerr= (10**Song_z7[:,1] - 10**Song_z7[:,3], 10**Song_z7[:,2] - 10**Song_z7[:,1]), xerr = 0.25, capsize = caps, alpha=0.75, elinewidth = errorwidth, lw=1.0, marker='o', ls='none', label = 'Song 2015, z = 7', color = 'blue')

        Song_z8 = np.array([[7.25, -1.73, -1.73 + 1.01, -1.73 - 0.84],
                            [7.75, -2.28, -2.28 + 0.84, -2.28 - 0.64],
                            [8.25, -2.88, -2.88 + 0.75, -2.88 - 0.57],
                            [8.75, -3.45, -3.45 + 0.57, -3.45 - 0.60],
                            [9.25, -4.21, -4.21 + 0.63, -4.21 - 0.78],
                            [9.75, -5.31, -5.31 + 1.01, -5.31 - 1.64]], dtype = np.float32)


        plt.errorbar(Song_z8[:,0], 10**Song_z8[:,1], yerr= (10**Song_z8[:,1] - 10**Song_z8[:,3], 10**Song_z8[:,2] - 10**Song_z8[:,1]), xerr = 0.25, capsize = caps, alpha=0.75, elinewidth = errorwidth, lw=1.0, marker='o', ls='none', label = 'Song 2015, z = 8', color = 'green')

    if (Observations == 2):
		Baldry = np.array([
		    [7.05, 1.3531e-01, 6.0741e-02],
		    [7.15, 1.3474e-01, 6.0109e-02],
		    [7.25, 2.0971e-01, 7.7965e-02],
		    [7.35, 1.7161e-01, 3.1841e-02],
		    [7.45, 2.1648e-01, 5.7832e-02],
		    [7.55, 2.1645e-01, 3.9988e-02],
		    [7.65, 2.0837e-01, 4.8713e-02],
		    [7.75, 2.0402e-01, 7.0061e-02],
		    [7.85, 1.5536e-01, 3.9182e-02],
		    [7.95, 1.5232e-01, 2.6824e-02],
		    [8.05, 1.5067e-01, 4.8824e-02],
		    [8.15, 1.3032e-01, 2.1892e-02],
		    [8.25, 1.2545e-01, 3.5526e-02],
		    [8.35, 9.8472e-02, 2.7181e-02],
		    [8.45, 8.7194e-02, 2.8345e-02],
		    [8.55, 7.0758e-02, 2.0808e-02],
		    [8.65, 5.8190e-02, 1.3359e-02],
		    [8.75, 5.6057e-02, 1.3512e-02],
		    [8.85, 5.1380e-02, 1.2815e-02],
		    [8.95, 4.4206e-02, 9.6866e-03],
		    [9.05, 4.1149e-02, 1.0169e-02],
		    [9.15, 3.4959e-02, 6.7898e-03],
		    [9.25, 3.3111e-02, 8.3704e-03],
		    [9.35, 3.0138e-02, 4.7741e-03],
		    [9.45, 2.6692e-02, 5.5029e-03],
		    [9.55, 2.4656e-02, 4.4359e-03],
		    [9.65, 2.2885e-02, 3.7915e-03],
		    [9.75, 2.1849e-02, 3.9812e-03],
		    [9.85, 2.0383e-02, 3.2930e-03],
		    [9.95, 1.9929e-02, 2.9370e-03],
		    [10.05, 1.8865e-02, 2.4624e-03],
		    [10.15, 1.8136e-02, 2.5208e-03],
		    [10.25, 1.7657e-02, 2.4217e-03],
		    [10.35, 1.6616e-02, 2.2784e-03],
		    [10.45, 1.6114e-02, 2.1783e-03],
		    [10.55, 1.4366e-02, 1.8819e-03],
		    [10.65, 1.2588e-02, 1.8249e-03],
		    [10.75, 1.1372e-02, 1.4436e-03],
		    [10.85, 9.1213e-03, 1.5816e-03],
		    [10.95, 6.1125e-03, 9.6735e-04],
		    [11.05, 4.3923e-03, 9.6254e-04],
		    [11.15, 2.5463e-03, 5.0038e-04],
		    [11.25, 1.4298e-03, 4.2816e-04],
		    [11.35, 6.4867e-04, 1.6439e-04],
		    [11.45, 2.8294e-04, 9.9799e-05],
		    [11.55, 1.0617e-04, 4.9085e-05],
		    [11.65, 3.2702e-05, 2.4546e-05],
		    [11.75, 1.2571e-05, 1.2571e-05],
		    [11.85, 8.4589e-06, 8.4589e-06],
		    [11.95, 7.4764e-06, 7.4764e-06],
		    ], dtype=np.float32)
		

		Baldry_xval = np.log10(10 ** Baldry[:, 0]  /AllVars.Hubble_h/AllVars.Hubble_h)
		Baldry_xval = Baldry_xval - 0.26  # convert back to Chabrier IMF
		Baldry_yvalU = (Baldry[:, 1]+Baldry[:, 2]) * AllVars.Hubble_h*AllVars.Hubble_h*AllVars.Hubble_h
		Baldry_yvalL = (Baldry[:, 1]-Baldry[:, 2]) * AllVars.Hubble_h*AllVars.Hubble_h*AllVars.Hubble_h

		plt.fill_between(Baldry_xval, Baldry_yvalU, Baldry_yvalL, 
		    facecolor='purple', alpha=0.25, label='Baldry et al. 2008 (z=0.1)')



    leg = plt.legend(loc='lower left', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(talk_legendsize)

    plt.tight_layout()

    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()
 

#####

def CentralGalaxy_Comparison(Simulation, Redshift, Mass, Photons):

    min_mass_H = 9
    max_mass_H = 12
    min_mass_G = 9
    max_mass_G = 12
    min_Photons = 49
    max_Photons = 53

    if Simulation == 0:
        OutputFile = 'CentralStellar_Photons_Fiducial_Millennium'
    elif Simulation == 1:
        OutputFile = 'CentralStellar_Photons_Fiducial_MySim_NoReion'

    PlotScripts.Plot_Scatter_SixPanel(Mass, Photons, 0, 1, ['Halos', 'Galaxies'], [min_mass_H, max_mass_H, min_Photons, max_Photons], [r'Log Halo Mass [$M_{\odot}$]', r'Log Ionizing Photons [s$^{-1}$]'], 2, Redshift, OutputFile, '.png')
    #PlotScripts.Plot_HeatMap_SixPanel(Mass, Photons, 0, 1, ['Halos', 'Galaxies'], [min_mass_H, max_mass_H, min_Photons, max_Photons], [r'Log Halo Mass [$M_{\odot}$]', r'Log Ionizing Photons [s$^{-1}$]'], 2, Redshift, OutputFile, '.png')

def CentralGalaxy_Comparison_Difference(Simulation,Redshift, Mass, Photons):

    min_mass_H = 9
    max_mass_H = 12
    min_mass_G = 9
    max_mass_G = 12

    if Simulation == 0:
        OutputFile = 'CentralStellar_Photons_Fiducial_Millennium_Difference'
    elif Simulation == 1:
        OutputFile = 'CentralStellar_Photons_Fiducial_MySim_Difference'

    PlotScripts.Plot_Scatter_SixPanel(Mass, Photons, 1, 1, ['Halos', 'Galaxies'], [min_mass_H, max_mass_H, -2, 1], [r'Log Halo Mass [$M_{\odot}$]', r'Log (Halo/Galaxy Ionizing Photons [s$^{-1}$]'], 2, Redshift, OutputFile, '.png')

##

def CentralGalaxy_Projection(Simulation, Redshift, Mass, Photons):

    if Simulation == 0:
        OutputFile = 'CentralStellar_Photons_Fiducial_Millennium_HistProjection'
    elif Simulation == 1:
        OutputFile = 'CentralStellar_Photons_Fiducial_MySim_HistProjection'

    title = []

    for i in xrange(0, len(Redshift)):
        tmp = 'z = %.2f: 512' %(Redshift[i])
        title.append(tmp)

    for i in xrange(0, len(Redshift)):
        tmp = 'z = %.2f: 1024' %(Redshift[i])
        title.append(tmp)

    PlotScripts.Plot_Scatter_TwoHists(Mass, Photons, title, [9, 12, 51, 54], [0, 0.30, 0, 0.30], 0.2, 1, 1, [r'Log Halo Mass [$M_{\odot}$]', r'Log Ionizing Photons [s$^{-1}$]'], 0, OutputFile, '.png') 

##

def Metallicity(Simulation, Redshift, Mass, Z):

    title = []
    for i in xrange(0, len(Redshift)):
        tmp = 'z = %.2f' %(Redshift[i])
        title.append(tmp)

    PlotScripts.Plot_Scatter(Mass, Z, title, [6, 11, 7.0, 9.5], [0, 0], [r'$\log_{10}\ M_{\mathrm{*}}\ (M_{\odot})$', r'$12\ +\ \log_{10}[\mathrm{O/H}]$'], 4, 1, "Z", ".png")
    PlotScripts.Plot_Scatter_TwoHists(Mass, Z, title, [6, 11, 7.0, 9.5], [0, 0.1, 0, 0.15], 0.1, 1, 1, [r'$\log_{10}\ M_{\mathrm{*}}\ (M_{\odot})$', r'$12\ +\ \log_{10}[\mathrm{O/H}]$'], 1, 'Z_Hist', '.png')

##

def HaloMassFunction(Simulation, Redshift, Mass, Mysim_Len): 

    title = []
    normalization = []
    binwidth = 0.1

    
    for i in xrange(0, Mysim_Len): 
        tmp = 'Mysim 512: z = %.2f' %(Redshift[i])
        title.append(tmp)

        AllVars.Set_Params_Mysim()
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
    
        normalization.append(norm)
 
    for i in xrange(0, Mysim_Len): 
        tmp = 'Mysim 1024: z = %.2f' %(Redshift[i])
        title.append(tmp)

        AllVars.Set_Params_Mysim()
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
    
        normalization.append(norm)

    '''    
    for i in xrange(Mysim_Len, len(Redshift)):

        tmp = 'Millennium: z = %.2f' %(Redshift[i])
        title.append(tmp)

        AllVars.Set_Params_MiniMill()
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
    
        normalization.append(norm)
    ''' 
    if Simulation == 0:
        OutputFile = 'HaloMassFunction_Millennium'
    elif Simulation == 1:
        OutputFile = 'HaloMassFunction_Comparison'
    elif Simulation ==2:   
        OutputFile = 'HaloMassFunction_Comparison'

    PlotScripts.Plot_Histogram(Mass, title, [7, 13.5, 1e-6, 1e-1], binwidth, normalization, 0, [0, 1], [r'$\log_{10}\ M_{\mathrm{Vir}}\ (M_{\odot})$', r'$\Phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$'], 1, 0, OutputFile, '.png')

##

def Central_Galaxy_Projection(Simulation, Redshift, Mass, Photons):

    if Simulation == 0:
        OutputFile = 'CentralStellar_Photons_Fiducial_Millennium_HistProjection'
    elif Simulation == 1:
        OutputFile = 'CentralStellar_Photons_Fiducial_MySim_HistProjection'



    PlotScripts.Plot_Scatter_TwoHists(Mass, Photons, title, [9, 12, 52, 55], [0, 0.25, 0, 0.25], 0.2, 1, 1, [r'Log Halo Mass [$M_{\odot}$]', r'Log Ionizing Photons [s$^{-1}$]'], 0, OutputFile, '.png') 



##

def FilteringMass(Simulation, SnapListZ, Mass, Sobacchi_idx, MySim_Len, model_tags):

    

 
    Normalization = []
 

    ## Plot Parameters ##
    binwidth = 0.1
    Observations = 0
    Output_Tag = "Filtering_MHpos_Fiducial"

    ##

    Mean_Array = []
    Std_Array = []

    Gnedin_Mean = np.zeros((len(SnapListZ)))
    Gnedin_Std = np.zeros((len(SnapListZ)))

    Sobacchi_Mean = np.zeros((len(SnapListZ)))
    Sobacchi_Std = np.zeros((len(SnapListZ)))

    Sobacchi_Mean2 = np.zeros((len(SnapListZ)))
    Sobacchi_Std2 = np.zeros((len(SnapListZ)))

    Mvir_Mean = np.zeros((len(SnapListZ)))
    Mvir_Std = np.zeros((len(SnapListZ)))

    Mvir_Mean2 = np.zeros((len(SnapListZ)))
    Mvir_Std2 = np.zeros((len(SnapListZ)))

    print "SNAPSHOTS", SnapListZ
    for i in xrange(0, len(Mass)):

	Mean = np.zeros((len(SnapListZ)))
    	for j in xrange(0, len(SnapListZ)):
		if(i in Sobacchi_idx == True): 
			Mean[j] = np.mean(Mass[i][j])
			Std[j] = np.std(Mass[i][j])
		else:
			y = 10
	ax1.plot(SnapListZ, np.log10(Mean), color = colors[i], linestyle = linestyles[i], label = model_tags[i], lw = 3)
	ax1.fill_between(SnapListZ, np.log10(Mean) - 0.434*Std/Mean, np.log10(Mean) + 0.434*Std/Mean, alpha = 0.3, color = colors[i])






	w = np.where(Sobacchi[i] > 1.0)[0]

	Sobacchi_Mean[i] = np.mean(Sobacchi[i][w])
	Sobacchi_Std[i] = np.std(Sobacchi[i][w])


	w2 = np.where(Sobacchi2[i] > 1.0)[0]

	print "w2", w2
	Sobacchi_Mean2[i] = np.mean(Sobacchi2[i][w2])
	Sobacchi_Std2[i] = np.std(Sobacchi2[i][w2])

	print "Mvir min", min(Mvir[i])

    fig = plt.figure()
    ax = plt.subplot(111)
          
    ax.plot(SnapListZ, np.log10(Gnedin_Mean), color = 'r', linestyle = '--', label = model_tags[0], lw = 3)
    print "Gnedin_Mean", Gnedin_Mean

    print "Sobacchi_Mean", np.log10(Sobacchi_Mean)
    print "Sobacchi_Std", np.log10(Sobacchi_Std)
    ax.plot(SnapListZ, np.log10(Sobacchi_Mean), color = 'b', linestyle = '-', label = model_tags[1], lw = 3)
    ax.fill_between(SnapListZ, np.log10(Sobacchi_Mean)-0.434*Sobacchi_Std/Sobacchi_Mean, np.log10(Sobacchi_Mean)+0.434*Sobacchi_Std/Sobacchi_Mean, alpha = 0.3, color = 'b')	


    ax.plot(SnapListZ, np.log10(Sobacchi_Mean2), color = 'g', linestyle = '-', label = model_tags[2])
    ax.fill_between(SnapListZ, np.log10(Sobacchi_Mean2)-0.434*Sobacchi_Std2/Sobacchi_Mean2, np.log10(Sobacchi_Mean2)+0.434*Sobacchi_Std2/Sobacchi_Mean2, alpha = 0.3, color = 'g')	

    print "Mvir_Mean", np.log10(Mvir_Mean)
    print "Mvir_Std", np.log10(Mvir_Std)
    ax.plot(SnapListZ, np.log10(Mvir_Mean), color = 'k', linestyle = '-', label = model_tags[3])
    ax.fill_between(SnapListZ, np.log10(Mvir_Mean)-0.434*Mvir_Std/Mvir_Mean, np.log10(Mvir_Mean)+0.434*Mvir_Std/Mvir_Mean, alpha = 0.2, color = 'k')

    #ax.plot(SnapListZ, np.log10(Mvir_Mean2), color = 'c', linestyle = '-', label = model_tags[4])
    #ax.fill_between(SnapListZ, np.log10(Mvir_Mean2)-0.434*Mvir_Std2/Mvir_Mean2, np.log10(Mvir_Mean2)+0.434*Mvir_Std2/Mvir_Mean2, alpha = 0.2, color = 'c')
 
    ax.set_xlabel(r"z", fontsize = label_size)
    ax.set_ylabel(r"$\mathrm{log}_{10} \: M_\odot$", fontsize = label_size)

    ax.set_xlim([5.5,13])
    ax.set_ylim([6,13])

    ax.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    ax.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))

    leg = plt.legend(loc= 1, numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')


    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()


##

def StarFormationRate_Redshift(Simulation, SnapListZ, SFR1, SFR2):

    title = [r"$\mathrm{Noreion Reionmine}$", r"$\mathrm{Noreion}$"]
 
    ## Plot Parameters ##
    binwidth = 0.1
    Observations = 0
    Output_Tag = "SFR"

    ##

    SFR1_Mean = np.zeros((len(SnapListZ))) 
    SFR1_Std = np.zeros((len(SnapListZ)))

    SFR2_Mean = np.zeros((len(SnapListZ)))
    SFR2_Std = np.zeros((len(SnapListZ)))

    print "SFR1", SFR1
    print "SFR2", SFR2

    for i in xrange(0, len(SnapListZ)):

	SFR1_Mean[i] = np.mean(10**SFR1[i])
	SFR1_Std[i] = np.std(10**SFR1[i])

	SFR2_Mean[i] = np.mean(10**SFR2[i])
	SFR2_Std[i] = np.std(10**SFR2[i])

 
    fig = plt.figure()
    ax = plt.subplot(111)
          
    ax.plot(SnapListZ, np.log10(SFR1_Mean), color = 'r', linestyle = '-', label = title[0])
    ax.fill_between(SnapListZ, np.log10(SFR1_Mean)-0.434*SFR1_Std/SFR1_Mean, np.log10(SFR1_Mean)+0.434*SFR1_Std[i]/SFR1_Mean, alpha = 0.5, color = 'r')	

    print "SFR1_Mean", np.log10(SFR1_Mean)
    print "SFR1_Std", np.log10(SFR1_Std)

    ax.plot(SnapListZ, np.log10(SFR2_Mean), color = 'b', linestyle = '-', label = title[1])
    ax.fill_between(SnapListZ, np.log10(SFR2_Mean)-0.434*SFR2_Std/SFR2_Mean, np.log10(SFR2_Mean)+0.434*SFR2_Std[i]/SFR2_Mean, alpha = 0.5, color = 'b')	

    print "SFR2_Mean", np.log10(SFR2_Mean)
    print "SFR2_Std", np.log10(SFR2_Std)

    ax.set_xlabel(r"z", fontsize = label_size + extra_size)
    ax.set_ylabel(r"$\mathrm{SFR} \: [M_\odot\mathrm{yr}^{-1}]$", fontsize = label_size + extra_size)

    ax.set_xlim([5.5,13])
    ax.set_ylim([-3,2.5])

    ax.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    ax.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))

    leg = plt.legend(loc= 1, numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')


    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()


##

def PhotonsStellarMass(Simulation, SnapListZ, Mass, Photons):

	Output_Tag = "PhotonsStellarMass"

	bound1 = 7 
	bound2 = 8
	bound3 = 9
	bound4 = 10
	bound5 = 11  

	bounds = [7, 8, 9, 10, 100]
	labels = [r"$7 \leq M_* \le 8$", r"$8 \leq M_* \le 9$", r"$9 \leq M_* \le 10$", r"$M_* \geq 10$"]

	photons_mean = []
	photons_std = []

	for i in xrange(0, len(bounds)-1):
		for j in xrange(0, len(SnapListZ)):
			w = np.where((Mass[j] >= bounds[i]) & (Mass[j] < bounds[i+1]))[0]
#			print "Redshift %.4f, bounds %d-%d" %(SnapListZ[j], bounds[i], bounds[i+1])
#			print "len(w)", len(w)
#			print "W", w
#			print "Photons[j][w]", Photons[j][w]

			if (len(w) > 0):
				photons_mean.append(Sum_Log(Photons[j][w]))
#				photons_std.append(Std_Log(Photons[j][w], Sum_Log(Photons[j][w])/len(w)))
			else:
				photons_mean.append(nan)
				photons_std.append(nan)

			print "photons_mean", photons_mean[i*len(SnapListZ) + j]

	## Photons are now in an array of [Bin1, Bin1, Bin1, Bin2, Bin2, Bin2, Bin3, Bin3, Bin3]
	##                                SNAP1,SNAP2,SNAP3,SNAP1,SNAP2,SNAP3,SNAP1,SNAP2,SNAP3

	ax = plt.subplot(111)

	for i in xrange(0, len(bounds)-1):
		print "For bounds ", bounds[i], "-", bounds[i+1], "photons_mean = ", photons_mean[i*len(SnapListZ):(i+1)*len(SnapListZ)]
		low_index = i*len(SnapListZ)
		high_index = (i+1)*len(SnapListZ)

		ax.plot(SnapListZ, np.log10(photons_mean[low_index:high_index]), color = colours[i], label = labels[i]) 
#    		ax.fill_between(SnapListZ, np.log10(photons_mean[low_index:high_index])-0.434*photons_std[low_index:high_index]/photons_mean[low_index:high_index], np.log10(photons_mean[low_index:high_index])+0.434*photons_std[low_index:high_index]/photons_mean[low_index:high_index], alpha = 0.5, color = colours[i])	


	ax.set_xlabel(r"z", fontsize = label_size + extra_size)
	ax.set_ylabel(r'Total $\dot{N}_{\gamma,\mathrm{HI}} \: [\mathrm{s}^{-1}]$', fontsize = label_size + extra_size)
    
	ax.set_xlim([5.5,13])
	ax.set_ylim([50, 60])

	
        leg = plt.legend(loc= 1, numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
        	t.set_fontsize('medium')


        outputFile = './' + Output_Tag + Output_Format
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()


##

def LymanAlphaLF(Simulation, Redshift, LymanAlpha, MySim_Len):

    title = []
    Normalization = []

    print LymanAlpha
    ## Plot Parameters ##
    binwidth = 0.05
    Observations = 1
    Output_Tag = "Test_LymanAlpha"
    Frequency = 0 # 0 for a frequency (count) histogram, 1 for a probaility histogram.
    ##

    ## Normalization and Titles for MySim ##
    for i in xrange(0, MySim_Len): 
        tmp = '512: z = %.2f' %(Redshift[i])
        title.append(tmp)

        AllVars.Set_Params_Mysim()
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
    
        Normalization.append(norm)


    ## Normalization and Titles for MySim2 ## 
    for i in xrange(0, MySim_Len): 
        tmp = '1024: z = %.2f' %(Redshift[i])
        title.append(tmp)

        AllVars.Set_Params_Mysim()
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
    
        Normalization.append(norm)

### Plotting ###

    plt.figure()  
    ax = plt.subplot(111)  

### Plots the Histograms ###

    for i in xrange(0, len(LymanAlpha)):

        (counts, Bin_Edges, Bin_Middle) = Calculate_Histogram(LymanAlpha[i], binwidth, Frequency)
        if (i < MySim_Len):
            ls = '-'
        else:
            ls = '--'
        plt.plot(Bin_Middle, counts / Normalization[i], colours[i], linestyle = ls, label = title[i])

##


## 

    plt.yscale('log', nonposy='clip')

    plt.axis([42, 44, 1e-6, 1e-1])

    ax.set_xlabel(r'$\log_{10}\ L_{\mathrm{Ly}\alpha}\ [\mathrm{erg} \: \mathrm{s}^{-1}]$')
    ax.set_ylabel(r'$\Phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$')
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

### If we want to put observations on the figure ###

    if Observations == 1:
    	Konno_2014_z73 = np.array([[42.5, 2.5e-4, 5e-4, 1e-4],
                               	[42.65, 9.5e-5, 2e-4, 5e-5],  
                               	[42.85, 2e-5, 3e-6, 7.5e-5]], dtype = np.float32) 

        plt.errorbar(Konno_2014_z73[:,0], Konno_2014_z73[:,1], yerr= (Konno_2014_z73[:,1] - Konno_2014_z73[:,3], Konno_2014_z73[:,2] - Konno_2014_z73[:,1]), xerr = 0.0, alpha = 1.0, lw=2.0, marker='o', ls='none', label = 'Konno et al. (2014), z = 7.3', color = '#bd0026')
   
    	Ouchi_2010_z66 = np.array([[42.5, 9.5e-4, 1.1e-3, 7e-4],
                               	[42.65, 5e-4, 5.5e-4, 4.0e-4],  
                               	[42.85, 2.5e-4, 3.5e-4, 2e-4],
				[43.05, 5.5e-5, 9.0e-5, 1.5e-5]], dtype = np.float32)
#				[43.50, 9.0e-6, 0.8e-5, 5.0e-5]], dtype = np.float32) 

        plt.errorbar(Ouchi_2010_z66[:,0], Ouchi_2010_z66[:,1], yerr= (Ouchi_2010_z66[:,1] - Ouchi_2010_z66[:,3], Ouchi_2010_z66[:,2] - Ouchi_2010_z66[:,1]), xerr = 0.0, alpha = 1.0, lw=2.0, marker='o', ls='none', label = 'Ouchi et al. (2010), z = 6.6', color = 'b')

 
    leg = plt.legend(loc=1, numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')

    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()
 
##

def PhotonsVsStellarMass_OLD(Simulation, SnapListZ, Mass, Photons):
   
    print "Plotting a Six panel Ngamma vs Stellar Mass plot."
 
    nrows = 2
    ncols = 3

    binwidth = 0.1
    low_mass = 7
    high_mass = 12
   
    Output_Tag = "PhotonsVsStellarMass"

    Limits = [low_mass, high_mass, 53, 57]

    bins = np.arange(low_mass,high_mass, binwidth)
    bins_mid = bins + binwidth/2.0
    bins_mid = bins_mid[:-1] # len(bins_mid) should be 1 less than len(bins) as the last bin doesn't have a midpoint.


    fig = plt.figure()
    for i in xrange(1, nrows*ncols + 1):
	print "On panel %d" %(i)
        ax = fig.add_subplot(nrows, ncols, i)
    
	Photons_Sum = []

	for j in xrange(0, len(bins)-1):
		w = np.where((Mass[i-1] >= bins[j]) & (Mass[i-1] < bins[j+1]))[0]
		print "w", w	
		if (len(w) != 0):
			Photons_Sum.append(Sum_Log(Photons[i-1][w]))
		else:
			Photons_Sum.append(nan)
		print "Photons_Sum[-1]", Photons_Sum[-1]
	ax.plot(bins_mid, np.log10(Photons_Sum))
		
        ax.set_xticks(np.arange(Limits[0], Limits[1] + 1.0, 1))
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f'))
        ax.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
        
        ax.set_yticks(np.arange(Limits[2], Limits[3] + 0.5, 0.5))
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
        ax.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))

        if i < 4:
            ax.tick_params(axis='both', labelbottom = 'off')
        else:
            ax.xaxis.tick_bottom()
        if i > 3:
            ax.xaxis.get_major_ticks()[-1].set_visible(False)

        if i != 1 and i != 4:
            ax.tick_params(axis = 'y', labelleft = 'off')
        if i == 1:
            ax.yaxis.get_major_ticks()[0].set_visible(True)
        if i == 4:
            ax.yaxis.get_major_ticks()[-1].set_visible(False)

        if i == 6:
            ax.xaxis.get_major_ticks()[-1].set_visible(True)

        ax.set_xlim(Limits[0:2])
        ax.set_ylim(Limits[2:4])
        
        label = "z = %.2f" %(SnapListZ[i-1])
        ax.text(11.25, 56.0, label, fontsize = label_size)


    plt.tight_layout()
    fig.text(0.5, 0.01, "Mass", ha = 'center')
    fig.text(0.001, 0.5, "Photons", va = 'center', rotation = 'vertical')
    fig.subplots_adjust(hspace = 0.0, wspace = 0.0)
    
    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

##

def PhotonsVsStellarMass(Simulation, SnapListZ, Mass, Photons):
   
    print "Plotting histogram for Ngamma vs M*" 
 
    binwidth = 0.1
    low_mass = 7
    high_mass = 12
   
    Output_Tag = "PhotonsVsStellarMassTQhiehtibt"

    Limits = [low_mass, high_mass, 53, 57]

    bins = np.arange(low_mass,high_mass, binwidth)
    bins_mid = bins + binwidth/2.0
    bins_mid = bins_mid[:-1] # len(bins_mid) should be 1 less than len(bins) as the last bin doesn't have a midpoint.


    ax = plt.subplot(111) 
    for i in xrange(0, len(SnapListZ)):
    
	Photons_Sum = []

	print "REDSHIFT %.4f" %(SnapListZ[i])
	for j in xrange(0, len(bins)-1):
		w = np.where((Mass[i] >= bins[j]) & (Mass[i] < bins[j+1]))[0]
		if (len(w) != 0):
			Photons_Sum.append(Sum_Log(Photons[i][w]))
			print "Photons_Sum[-1]", Photons_Sum[-1], "Bins %.1f-%.1f" %(bins[j], bins[j+1])
		else:
			Photons_Sum.append(nan)


	label = "z = %.2f" %(SnapListZ[i])
	ax.plot(bins_mid, np.log10(Photons_Sum), label = label)
 
    ax.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval)) 
    ax.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))

    ax.set_xlabel(r'$\log_{10}\ M_{\mathrm{*}}\ (M_{\odot})$') 
    ax.set_ylabel(r'$\dot{N}_\gamma dM_\mathrm{*}$')
            
    ax.set_xlim(Limits[0:2])
    ax.set_ylim(Limits[2:4])

    leg = plt.legend(loc=1, numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')

       
    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()


##

def fesc(simulation, SnapListZ, fesc, model_tags, Output_Tag):

    ax1 = plt.subplot(111)

    t = np.empty(len(SnapListZ))
        
    for i in xrange(0, len(SnapListZ)):
	t[i] = (t_BigBang - cosmo.lookback_time(SnapListZ[i]).value) * 1.0e3   

    print SnapListZ
    print t
    print len(t)

    for j in xrange(0, len(fesc)):
	avg = []
	std = []
	
    	for i in xrange(0, len(SnapListZ)):
	    	avg.append(np.mean(fesc[j][i]))
	    	std.append(np.std(fesc[j][i]))
    	ax1.plot(t, avg, color = colors[j], ls = linestyles[j], label = model_tags[j], lw = 3)  
    	ax1.fill_between(t, np.subtract(avg,std), np.add(avg,std), color = colors[j], alpha = 0.25)

    ax1.xaxis.set_minor_locator(mtick.MultipleLocator(time_tick_interval))
    ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.025))
    ax1.set_xlim(time_xlim)
    ax1.set_ylim([-0.05, 0.7])

    ax2 = ax1.twiny()

    t_plot = (t_BigBang - cosmo.lookback_time(z_plot).value) * 1.0e3 # Corresponding Time values on the bottom.
    z_labels = ["$%d$" % x for x in z_plot] # Properly Latex-ize the labels.

    ax2.set_xlabel(r"$z$", size = label_size)
    ax2.set_xlim(time_xlim)
    ax2.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
    ax2.set_xticklabels(z_labels) # But label them as redshifts.

    ax1.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = label_size)
    ax1.set_ylabel(r'$f_\mathrm{esc}$', fontsize = talk_fontsize)

    leg = ax1.legend(loc=1, numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(talk_legendsize)

    plt.tight_layout()
    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

##

def UVLF(Simulation, Redshift, MUV, MUV_Obs, MySim_Len, tag, Output_Tag): 

    title = []
    Normalization = []

    ## Plot Parameters ##
    binwidth = 0.2
    Observations = 1
    Frequency = 0 # 0 for a frequency (count) histogram, 1 for a probaility histogram.
    errorwidth = 2
    delta = 0.00

    caps = 5
    ##

    ## Normalization and Titles for MySim ##
    for i in xrange(0, MySim_Len): 
        tmp = '%s: z = %.2f' %(tag, Redshift[i])
        title.append(tmp)

        AllVars.Set_Params_Mysim()
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
    
        Normalization.append(norm)

### Plotting ###

    nrows = 2
    ncols = 3

    fig = plt.figure()  
    Limits = [-23, -16, 1e-6, 1e-1]
### Plots the Histograms ###

    for i in xrange(1, 7):
	print "On panel %d" %(i)
        ax = fig.add_subplot(nrows, ncols, i)

        (counts, Bin_Edges, Bin_Middle) = Calculate_Histogram(MUV[i-1], binwidth, Frequency) 
        ax.plot(Bin_Middle, counts / Normalization[i-1], color = 'r', linestyle = '-') 

        (counts, Bin_Edges, Bin_Middle) = Calculate_Histogram(MUV_Obs[i-1], binwidth, Frequency)
        ax.plot(Bin_Middle, counts / Normalization[i-1], color = 'r', linestyle = '--', lw = 2)
	
        ax.set_xticks(np.arange(Limits[0], Limits[1] + 1.0, 1))
        #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f'))
        ax.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
        
#        ax.set_yticks(np.arange(Limits[2], Limits[3] + 0.5, 0.5))
#        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
#        ax.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
	ax.set_yscale('log', nonposy='clip')

        if i < 4:
            ax.tick_params(axis='both', labelbottom = 'off')
        else:
            ax.xaxis.tick_bottom()
        if i > 3:
            ax.xaxis.get_major_ticks()[-1].set_visible(False)

        if i != 1 and i != 4:
            ax.tick_params(axis = 'y', labelleft = 'off')
        if i == 1:
            ax.yaxis.get_major_ticks()[0].set_visible(False)
        if i == 4:
            ax.yaxis.get_major_ticks()[-1].set_visible(False)

        if i == 6:
            ax.xaxis.get_major_ticks()[-1].set_visible(True)

        ax.set_xlim(Limits[0:2])
        ax.set_ylim(Limits[2:4])
        
        label = "z = %.2f" %(SnapListZ[i-1])
        ax.text(-22.8, 5e-2, label, fontsize = label_size)

	if (i == 1):

    		## z = 5 ##	
		Obs = np.array([[-23.11 + delta, 0.000002, 0.000002 + 0.000002, 0.000000000001],
                                [-22.61 + delta, 0.000006, 0.000006 + 0.000003, 0.000006 - 0.000003],
                                [-22.11 + delta, 0.000034, 0.000034 + 0.000008, 0.000034 - 0.000008],
                                [-21.61 + delta, 0.000101, 0.000101 + 0.000014, 0.000101 - 0.000014],
                                [-21.11 + delta, 0.000265, 0.000265 + 0.000025, 0.000265 - 0.000025],
                                [-20.61 + delta, 0.000676, 0.000676 + 0.000046, 0.000676 - 0.000046],
				[-20.11 + delta, 0.001029, 0.001029 + 0.000067, 0.001029 - 0.000067],
				[-19.61 + delta, 0.001329, 0.001329 + 0.000094, 0.001329 - 0.000094],
				[-18.36 + delta, 0.004460, 0.004460 + 0.000540, 0.004460 - 0.000540],
				[-17.36 + delta, 0.008600, 0.008600 + 0.001760, 0.008600 - 0.001760],
				[-16.36 + delta, 0.024400, 0.024400 + 0.007160, 0.024400 - 0.007160]], dtype = np.float32)


		## Intrinsic ##
    		fname = "/home/jseiler/SAGE-stuff/output/lfintrinsic/lfint_z5.dat"
		fd = open(fname, 'rb')
		M_i, Phi_i = np.loadtxt(fd, dtype = np.float32, skiprows=1, unpack=True)  


		## Observed ##
    		fname = "/home/jseiler/SAGE-stuff/output/lfobserved/lf_z5.dat"
		fd = open(fname, 'rb')
		M_O, Phi_O = np.loadtxt(fd, dtype = np.float32, skiprows=1, unpack=True)  
	
	if (i == 2):

    		## z = 6 ##	
		Obs = np.array([[-22.52 + delta, 0.000002, 0.000002 + 0.000002, 0.000000000001],
                                [-22.02 + delta, 0.000015, 0.000015 + 0.000006, 0.000015 - 0.000006],
                                [-21.52 + delta, 0.000043, 0.000043 + 0.000012, 0.000043 - 0.000012],
                                [-21.02 + delta, 0.000176, 0.000176 + 0.000025, 0.000176 - 0.000025],
                                [-20.52 + delta, 0.000320, 0.000320 + 0.000041, 0.000320 - 0.000041],
                                [-19.52 + delta, 0.000698, 0.000698 + 0.000083, 0.000698 - 0.000083],
				[-18.77 + delta, 0.001900, 0.001900 + 0.000320, 0.001900 - 0.000320],
				[-17.77 + delta, 0.006680, 0.006680 + 0.001380, 0.006680 - 0.001380],
				[-16.77 + delta, 0.013640, 0.013640 + 0.004200, 0.013640 - 0.004200]], dtype = np.float32)		

		## Intrinsic ##
    		fname = "/home/jseiler/SAGE-stuff/output/lfintrinsic/lfint_z6.dat"
		fd = open(fname, 'rb')
		M_i, Phi_i = np.loadtxt(fd, dtype = np.float32, skiprows=1, unpack=True)  


		## Observed ##
    		fname = "/home/jseiler/SAGE-stuff/output/lfobserved/lf_z6.dat"
		fd = open(fname, 'rb')
		M_O, Phi_O = np.loadtxt(fd, dtype = np.float32, skiprows=1, unpack=True)  	


	if (i == 3):
		## z = 7 ##

		NoBound = np.array([[-22.66 + delta, 0.000002]], dtype = np.float32) 
		if (Observations == 1):
			ax.errorbar(NoBound[:,0], NoBound[:,1], yerr = 0.000001, uplims = True, alpha=0.75, lw=1.0, marker='o', linestyle='none', color = 'cyan')

    		Obs = np.array([[-22.16 + delta, 0.000001, 0.000001 + 0.000002, 0.000002 - 0.000001],
                                [-21.66 + delta, 0.000033, 0.000033 + 0.000009, 0.000033 - 0.000009],
                                [-21.16 + delta, 0.000048, 0.000048 + 0.000015, 0.000048 - 0.000015],
                                [-20.66 + delta, 0.000193, 0.000193 + 0.000034, 0.000193 - 0.000034],
				[-20.16 + delta, 0.000309, 0.000309 + 0.000061, 0.000309 - 0.000061],
                                [-19.66 + delta, 0.000654, 0.000654 + 0.000100, 0.000654 - 0.000100],
				[-19.16 + delta, 0.000907, 0.000907 + 0.000117, 0.000907 - 0.000177],
				[-18.66 + delta, 0.001717, 0.001717 + 0.000478, 0.001717 - 0.000478],
				[-17.91 + delta, 0.005840, 0.005840 + 0.001460, 0.005840 - 0.001460],
				[-16.91 + delta, 0.008500, 0.008500 + 0.002940, 0.008500 - 0.002940]], dtype = np.float32)

		## Intrinsic ##
    		fname = "/home/jseiler/SAGE-stuff/output/lfintrinsic/lfint_z7.dat"
		fd = open(fname, 'rb')
		M_i, Phi_i = np.loadtxt(fd, dtype = np.float32, skiprows=1, unpack=True)  


		## Observed ##
    		fname = "/home/jseiler/SAGE-stuff/output/lfobserved/lf_z7.dat"
		fd = open(fname, 'rb')
		M_O, Phi_O = np.loadtxt(fd, dtype = np.float32, skiprows=1, unpack=True)  	

	if (i == 4):
    		## z = 8 ##

		
		NoBound = np.array([[-22.87 + delta, 0.000002],
				    [-22.37 + delta, 0.000002]], dtype = np.float32) 
		if (Observations == 1):
			ax.errorbar(NoBound[:,0], NoBound[:,1], yerr = 0.000001, uplims = True, alpha=0.75, lw=1.0, marker='o', linestyle='none', color = 'cyan')
	
		Obs = np.array([[-21.87 + delta, 0.000005, 0.000005 + 0.000003, 0.000005 - 0.000003],
                                [-21.37 + delta, 0.000013, 0.000013 + 0.000005, 0.000013 - 0.000005],
                                [-20.87 + delta, 0.000058, 0.000058 + 0.000015, 0.000058 - 0.000015],
                                [-20.37 + delta, 0.000060, 0.000060 + 0.000026, 0.000060 - 0.000026],
                                [-19.87 + delta, 0.000331, 0.000331 + 0.000104, 0.000331 - 0.000104],
                                [-19.37 + delta, 0.000676, 0.000676 + 0.000046, 0.000676 - 0.000046],
				[-18.62 + delta, 0.001060, 0.001060 + 0.000340, 0.001060 - 0.000340],
				[-17.62 + delta, 0.002740, 0.002740 + 0.001040, 0.002740 - 0.001040]], dtype = np.float32)


		## Intrinsic ##
    		fname = "/home/jseiler/SAGE-stuff/output/lfintrinsic/lfint_z8.dat"
		fd = open(fname, 'rb')
		M_i, Phi_i = np.loadtxt(fd, dtype = np.float32, skiprows=1, unpack=True)  


		## Observed ##
    		fname = "/home/jseiler/SAGE-stuff/output/lfobserved/lf_z8.dat"
		fd = open(fname, 'rb')
		M_O, Phi_O = np.loadtxt(fd, dtype = np.float32, skiprows=1, unpack=True)  


	if (i == 5):
    		## z = 9 ##

		
		NoBound = np.array([[-21.94 + delta, 1.0e-3*0.0024]], dtype = np.float32) 
		if (Observations == 1):
			ax.errorbar(NoBound[:,0], NoBound[:,1], yerr = 1.0e-3*0.001, uplims = True, alpha=0.75, lw=1.0, marker='o', linestyle='none', color = 'cyan')
	
		Obs = np.array([[-21.14 + delta, 1.0e-3*0.004400, 1.0e-3*(0.004400 + 0.004200), 1.0e-3*(0.004400 - 0.002400)],
                                [-20.34 + delta, 1.0e-3*0.032200, 1.0e-3*(0.032200 + 0.021700), 1.0e-3*(0.032200 - 0.013800)]], dtype = np.float32)

		## Intrinsic ##
    		fname = "/home/jseiler/SAGE-stuff/output/lfintrinsic/lfint_z9.dat"
		fd = open(fname, 'rb')
		M_i, Phi_i = np.loadtxt(fd, dtype = np.float32, skiprows=1, unpack=True)  


		## Observed ##
    		fname = "/home/jseiler/SAGE-stuff/output/lfobserved/lf_z9.dat"
		fd = open(fname, 'rb')
		M_O, Phi_O = np.loadtxt(fd, dtype = np.float32, skiprows=1, unpack=True)  

	if (i == 6):
    		## z = 10 ##

		
		NoBound = np.array([[-22.23 + delta, 0.000001],
				    [-19.23 + delta, 0.000049]], dtype = np.float32) 
		if (Observations == 1):
			ax.errorbar(NoBound[:,0], NoBound[:,1], yerr = [0.000001, 0.00001], uplims = True, alpha=0.75, lw=1.0, marker='o', linestyle='none', color = 'cyan')
	
		Obs = np.array([[-21.23 + delta, 0.000001, 0.000001 + 0.000001, 0.000001 - 0.000001],
                                [-20.23 + delta, 0.000010, 0.000010 + 0.000005, 0.000010 - 0.000005],
                                [-18.23 + delta, 0.000266, 0.000266 + 0.000171, 0.000266 - 0.000171]], dtype = np.float32)

		## Intrinsic ##
    		fname = "/home/jseiler/SAGE-stuff/output/lfintrinsic/lfint_z10.dat"
		fd = open(fname, 'rb')
		M_i, Phi_i = np.loadtxt(fd, dtype = np.float32, skiprows=1, unpack=True)  


		## Observed ##
    		fname = "/home/jseiler/SAGE-stuff/output/lfobserved/lf_z10.dat"
		fd = open(fname, 'rb')
		M_O, Phi_O = np.loadtxt(fd, dtype = np.float32, skiprows=1, unpack=True)  



	if (Observations == 1):

		ax.plot(M_i, 10**Phi_i, color = 'k', linestyle = '-')
		ax.plot(M_O, 10**Phi_O, color = 'k', linestyle = '--', lw = 2)
		ax.errorbar(Obs[:,0], Obs[:,1], yerr= (Obs[:,1] - Obs[:,3], Obs[:,2] - Obs[:,1]) , alpha=0.75, lw=1.0, marker='o', linestyle='none', color = 'cyan')
		if (i == 1):

			ax.plot(1e10, 1e10, color = 'r', linestyle = '-', label = 'SAGE Intrinsic')
			ax.plot(1e10, 1e10, color = 'r', linestyle = '--', lw = 2, label = 'SAGE Observed')
			
			ax.plot(1e10, 1e10, color = 'k', linestyle = '-', label = 'DRAGONS Intrinsic')
			ax.plot(1e10, 1e10, color = 'k', linestyle = '--', lw = 2, label = 'DRAGONS Observed')

			ax.errorbar(1e10, 1e10, color = 'c', xerr = 1, yerr = 1, marker = 'o', label = "Bouwens et. al (2015)")
			
        		leg = ax.legend(loc = 4, scatterpoints=1, labelspacing=0.0)
    			leg.draw_frame(False)  # Don't want a box frame
			for t in leg.get_texts():  # Reduce the size of the text
			    	t.set_fontsize(label_size - 6)
##

    plt.tight_layout()
    fig.text(0.5, 0.01, r"$M_{1600}$", ha = 'center')
    fig.text(0.001, 0.5, r'$\Phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$', va = 'center', rotation = 'vertical')
    fig.subplots_adjust(hspace = 0.0, wspace = 0.0)
   
    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()
##

def SFR_Hist(Simulation, Redshift, SFR, MySim_Len):

    title = []
    Normalization = []

    ## Plot Parameters ##
    binwidth = 0.1
    Observations = 1
    Output_Tag = "SFR_Hist_reionmine"
    Frequency = 0 # 0 for a frequency (count) histogram, 1 for a probaility histogram.
    errorwidth = 2
    delta = 0.05
    caps = 5
    ##

    ## Normalization and Titles for MySim ##
    for i in xrange(0, MySim_Len): 
        tmp = 'Noreion reionmine MH Beta: z = %.2f' %(Redshift[i])
        title.append(tmp)

        AllVars.Set_Params_Mysim()
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
    
        Normalization.append(norm)


    ## Normalization and Titles for MySim2 ## 
    for i in xrange(0, MySim_Len): 
        tmp = 'Noreion reionmine fesc = 0.15: z = %.2f' %(Redshift[i])
        title.append(tmp)

        AllVars.Set_Params_Mysim()
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
    
        Normalization.append(norm)

### Plotting ###

    plt.figure()  
    ax = plt.subplot(111)  

### Plots the Histograms ###

    for i in xrange(0, len(SFR)):

        (counts, Bin_Edges, Bin_Middle) = Calculate_Histogram(SFR[i], binwidth, Frequency)
        if (i < MySim_Len):
            ls = '-'
        else:
            ls = '--'
        plt.plot(Bin_Middle, np.log10(counts / Normalization[i]), colours[i], linestyle = ls, label = title[i])

##

    #plt.yscale('log', nonposy='clip')

    plt.axis([-4, 4, -6, 0.2])

    ax.set_xlabel(r'$\log_{10}\ SFR [(M_{\odot}) \mathrm{yr}^{-1}$')
    ax.set_ylabel(r'$\Phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$')
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

    leg = plt.legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')

    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()
##

def SFRVsStellarMass(Simulation, SnapListZ, Mass, SFR):
   
    print "Plotting histogram for SFR vs M*" 
 
    binwidth = 0.1
    low_mass = 7
    high_mass = 12
   
    Output_Tag = "SFRVsStellarMassFraction_verylowSF"

    Limits = [low_mass, high_mass, 1e-3, 1.1]

    bins = np.arange(low_mass,high_mass, binwidth)
    bins_mid = bins + binwidth/2.0
    bins_mid = bins_mid[:-1] # len(bins_mid) should be 1 less than len(bins) as the last bin doesn't have a midpoint.


    ax = plt.subplot(111) 
    for i in xrange(0, len(SnapListZ)):
    
	SFR_Sum = []

	print "REDSHIFT %.4f" %(SnapListZ[i])
	pleb_sum = 0.0
	for j in xrange(0, len(bins)-1):
		w = np.where((Mass[i] >= bins[j]) & (Mass[i] < bins[j+1]))[0]
		if (len(w) != 0):
			SFR_Sum.append(Sum_Log(SFR[i][w]))
			pleb_sum += Sum_Log(SFR[i][w])
		else:
			SFR_Sum.append(np.nan)


	label = "z = %.2f" %(SnapListZ[i])
	ax.plot(bins_mid, SFR_Sum/pleb_sum, label = label)
 
	print "SFR_Sum = ", SFR_Sum, "pleb_sum", pleb_sum, "SFR_Sum/pleb_sum", np.log10(SFR_Sum/pleb_sum), "log(SFR_Sum) - log(pleb_sum)", np.log10(SFR_Sum) - np.log10(pleb_sum) 

 
    ax.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval)) 
    ax.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))

    ax.set_yscale('log', nonposy='clip')
    ax.set_xlabel(r'$\log_{10}\ M_{\mathrm{*}}\ [M_{\odot}]$') 
    #ax.set_ylabel(r'$\log_{10}\ \mathrm{SFR} \: [M_{\odot} \mathrm{yr}^{-1}]$')
    ax.set_ylabel(r'$\mathrm{SFR}/\mathrm{SFR}_\mathrm{Tot}$')
    
    ax.set_xlim(Limits[0:2])
    ax.set_ylim(Limits[2:4])

    leg = plt.legend(loc=1, numpoints=1, labelspacing=0.1)

    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')

    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()




##

def sSFR_Hist(Simulation, Redshift, sSFR, MySim_Len):

    title = []
    Normalization = []

    ## Plot Parameters ##
    binwidth = 0.1
    Observations = 1
    Output_Tag = "sSFR_Hist"
    Frequency = 0 # 0 for a frequency (count) histogram, 1 for a probaility histogram.
    errorwidth = 2
    delta = 0.05
    caps = 5
    ##

    ## Normalization and Titles for MySim ##
    for i in xrange(0, MySim_Len): 
        tmp = 'Noreion: z = %.2f' %(Redshift[i])
        title.append(tmp)

        AllVars.Set_Params_Mysim()
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
    
        Normalization.append(norm)


    ## Normalization and Titles for MySim2 ## 
    for i in xrange(0, MySim_Len): 
        tmp = 'Reionmine: z = %.2f' %(Redshift[i])
        title.append(tmp)

        AllVars.Set_Params_Mysim()
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
    
        Normalization.append(norm)

### Plotting ###

    plt.figure()  
    ax = plt.subplot(111)  

### Plots the Histograms ###

    for i in xrange(0, len(sSFR)):

        (counts, Bin_Edges, Bin_Middle) = Calculate_Histogram(sSFR[i], binwidth, Frequency)
        if (i < MySim_Len):
            ls = '-'
        else:
            ls = '--'
        plt.plot(Bin_Middle, np.log10(counts / Normalization[i]), colours[i], linestyle = ls, label = title[i])

##

    #plt.yscale('log', nonposy='clip')

#    plt.axis([-4, 4, -6, 0.2])

    ax.set_xlabel(r'$\log_{10}\ sSFR [\mathrm{yr}^{-1}$')
    ax.set_ylabel(r'$\Phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$')
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

    leg = plt.legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')

    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()
##

def EjectedFracVsStellarMass(Simulation, Redshift, mass, EjectedFraction, MySim_Len, Output_Tag):

    print "Ejected Mass vs Stellar Mass"
    title = []

    fig = plt.figure()

    ax1 = plt.subplot(211)

    binwidth = 0.1
    low_mass = 3
    high_mass = 14

    bins = np.arange(low_mass,high_mass, binwidth)
    bins_mid = bins + binwidth/2.0
    bins_mid = bins_mid[:-1] # len(bins_mid) should be 1 less than len(bins) as the last bin doesn't have a midpoint.

    for i in xrange(0, len(Redshift)):
	    ejected_sum = [] 
	    for j in xrange(0, len(bins)-1):
		w = np.where((mass[i] >= bins[j]) & (mass[i] < bins[j+1]))[0]
		if (len(w) != 0):
			ejected_sum.append(np.sum(EjectedFraction[i][w])/len(w))
		else:
			ejected_sum.append(nan)
	
	    ax1.scatter(np.mean(mass[i]), 0.3, s = 30, color = colors[i])

	    tmp = 'z = %.2f' %(Redshift[i])
	    ax1.plot(bins_mid, ejected_sum, color = colors[i], label = tmp) 
	
#    ax1.set_xlabel(r'$\log_{10}\ M_{\mathrm{H}}\ [M_{\odot}]$') 
    ax1.set_ylabel(r'$\mathrm{Ejected \: Fraction}$', size = label_size)
    ax1.set_xlim([8.5, 12])
    ax1.set_ylim([-0.05, 0.6])   

    ax1.xaxis.set_minor_locator(mtick.MultipleLocator(0.1))
    ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.025))
    ax1.set_xticklabels([])
   

    leg = ax1.legend(loc=1, numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')

    ax2 = plt.subplot(212)

    hb = ax2.hexbin(mass[3], EjectedFraction[3], cmap = 'inferno', gridsize = 50, bins = 'log')
    counts = hb.get_array()

    cax = fig.add_axes([0.95, 0.14, 0.03, 0.38])
    cbar = fig.colorbar(hb, cax=cax, ticks = np.arange(0, max(counts)+0.2, 0.4))
    cbar.ax.set_ylabel(r"$\log_{10} N$",  rotation = 90, fontsize = label_size) 
    cbar.ax.tick_params(labelsize = label_size - 2) 

    ax2.set_xlabel(r'$\log_{10}\ M_{\mathrm{H}}\ [M_{\odot}]$', size = label_size) 
    ax2.set_ylabel(r'$\mathrm{Ejected \: Fraction}$', size = label_size)

    ax2.set_xlim([8.5, 12])
    ax2.xaxis.set_minor_locator(mtick.MultipleLocator(0.1))
    ax2.yaxis.set_minor_locator(mtick.MultipleLocator(0.025)) 

    fig.subplots_adjust(right = None, hspace = 0.0, wspace = 0.0)
    plt.tight_layout()

    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

##

def HaloPartCount(Simulation, Redshift, HaloCount, MySim_Len): 

    ax = plt.subplot(111)
    binwidth = 1
    Frequency = 0

    Output_Tag = "HaloPartCount"

    for i in xrange(0, len(Redshift)):

        (counts, Bin_Edges, Bin_Middle) = Calculate_Histogram(HaloCount[i], binwidth, Frequency, 0, 200)
        if (i < MySim_Len):
            ls = '-'
        else:
            ls = '--'
	label = r"$z = %.2f$" %(Redshift[i])
        ax.plot(Bin_Middle, counts / 1, colours[i], linestyle = ls, label = label)

    ax.set_yscale('log', nonposy='clip')

    ax.set_xlabel(r'$\mathrm{Number \: Halo \: Particles}$')
    ax.set_ylabel(r'$\mathrm{Count}$')
    ax.xaxis.set_minor_locator(plt.MultipleLocator(10))


    leg = plt.legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
	t.set_fontsize('medium')

    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

##

def SFR(Simulation, Redshift, Mass, SFR, MySim_Len, Output_Tag):

    ax = plt.subplot(111)

    '''
    binwidth = 0.1
    Frequency = 0  

    for i in xrange(0, len(Redshift)):
 
        (counts, Bin_Edges, Bin_Middle) = Calculate_Histogram(SFR[i], binwidth, Frequency, -7, 3)
        if (i < MySim_Len):
            ls = '-'
        else:
            ls = '--'
	label = r"$z = %.2f$" %(Redshift[i])
        ax.plot(Bin_Middle, counts / 1, colors[i], linestyle = ls, label = label)

    ax.set_yscale('log', nonposy='clip')

    ax.set_xlabel(r'$\mathrm{SFR}$')
    ax.set_ylabel(r'$\mathrm{Count}$')
    ax.xaxis.set_minor_locator(plt.MultipleLocator(10))
    '''

    for i in xrange(0, len(Redshift)):
	label = r"$z = %.2f$" %(Redshift[i])
	ax.scatter(Mass, SFR, color = colors[i], alpha = 0.2, label = label)

    ax.set_xlabel(r'$m_*$', size = label_size)
    ax.set_ylabel(r'$\mathrm{SFR}$', size = label_size) 

    leg = plt.legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
	t.set_fontsize('medium')

    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

##

def BaryonFraction(G, G2):


	print 'Plotting the average baryon fraction vs halo mass'

        seed(2222)
    
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
        
        HaloMass = np.log10(G.Mvir * 1.0e10 / AllVars.Hubble_h)
        Baryons = G.StellarMass + G.ColdGas + G.HotGas + G.EjectedMass + G.IntraClusterStars + G.BlackHoleMass

        HaloMass2 = np.log10(G2.Mvir * 1.0e10 / AllVars.Hubble_h)
        Baryons2 = G2.StellarMass + G2.ColdGas + G2.HotGas + G2.EjectedMass + G2.IntraClusterStars + G2.BlackHoleMass


        MinHalo = 11.0 
        MaxHalo = 16.0
        Interval = 0.1
        Nbins = int((MaxHalo-MinHalo)/Interval)
        HaloRange = np.arange(MinHalo, MaxHalo, Interval)
        
        MeanCentralHaloMass = []
        MeanBaryonFraction = []
        MeanBaryonFractionU = []
        MeanBaryonFractionL = []

        MeanStars = []
        MeanCold = []
        MeanHot = []
        MeanEjected = []
        MeanICS = []
        MeanBH = []

        MeanCentralHaloMass2 = []
        MeanBaryonFraction2 = []
        MeanBaryonFractionU2 = []
        MeanBaryonFractionL2 = []

        MeanStars2 = []
        MeanCold2 = []
        MeanHot2 = []
        MeanEjected2 = []
        MeanICS2 = []
        MeanBH2 = []


        for i in xrange(Nbins-1):
            
            w1 = np.where((G.Type == 0) & (HaloMass >= HaloRange[i]) & (HaloMass < HaloRange[i+1]))[0]
            HalosFound = len(w1)
            
            if HalosFound > 2:  
                
                BaryonFraction = []
                CentralHaloMass = []
                
                Stars = []
                Cold = []
                Hot = []
                Ejected = []
                ICS = []
                BH = []
                
                for j in xrange(HalosFound):
                    
                    w2 = np.where(G.CentralGalaxyIndex == G.CentralGalaxyIndex[w1[j]])[0]
                    CentralAndSatellitesFound = len(w2)
                    
                    if CentralAndSatellitesFound > 0:
                        BaryonFraction.append(sum(Baryons[w2]) / G.Mvir[w1[j]])
                        CentralHaloMass.append(np.log10(G.Mvir[w1[j]] * 1.0e10 / AllVars.Hubble_h))

                        Stars.append(sum(G.StellarMass[w2]) / G.Mvir[w1[j]])
                        Cold.append(sum(G.ColdGas[w2]) / G.Mvir[w1[j]])
                        Hot.append(sum(G.HotGas[w2]) / G.Mvir[w1[j]])
                        Ejected.append(sum(G.EjectedMass[w2]) / G.Mvir[w1[j]])
                        ICS.append(sum(G.IntraClusterStars[w2]) / G.Mvir[w1[j]])
                        BH.append(sum(G.BlackHoleMass[w2]) / G.Mvir[w1[j]])                        
                                
                MeanCentralHaloMass.append(np.mean(CentralHaloMass))
                MeanBaryonFraction.append(np.mean(BaryonFraction))
                MeanBaryonFractionU.append(np.mean(Ejected) + np.var(Ejected))
                MeanBaryonFractionL.append(np.mean(Ejected) - np.var(Ejected))
                
                MeanStars.append(np.mean(Stars))
                MeanCold.append(np.mean(Cold))
                MeanHot.append(np.mean(Hot))
                MeanEjected.append(np.mean(Ejected))
                MeanICS.append(np.mean(ICS))
                MeanBH.append(np.mean(BH))

	    w1 = np.where((G2.Type == 0) & (HaloMass2 >= HaloRange[i]) & (HaloMass2 < HaloRange[i+1]))[0]
	    HalosFound = len(w1)

	    if HalosFound > 2:  
		    
	        BaryonFraction2 = []
		CentralHaloMass2 = []
		    
		Stars2 = []
		Cold2 = []
		Hot2 = []
		Ejected2 = []
		ICS2 = []
		BH2 = []
		    
		for j in xrange(HalosFound):
			
		    w2 = np.where(G2.CentralGalaxyIndex == G2.CentralGalaxyIndex[w1[j]])[0]
		    CentralAndSatellitesFound = len(w2)
			
		    if CentralAndSatellitesFound > 0:
			BaryonFraction2.append(sum(Baryons2[w2]) / G2.Mvir[w1[j]])
			CentralHaloMass2.append(np.log10(G2.Mvir[w1[j]] * 1.0e10 / AllVars.Hubble_h))

			Stars2.append(sum(G2.StellarMass[w2]) / G2.Mvir[w1[j]])
			Cold2.append(sum(G2.ColdGas[w2]) / G2.Mvir[w1[j]])
			Hot2.append(sum(G2.HotGas[w2]) / G2.Mvir[w1[j]])
			Ejected2.append(sum(G2.EjectedMass[w2]) / G2.Mvir[w1[j]])
			ICS2.append(sum(G2.IntraClusterStars[w2]) / G2.Mvir[w1[j]])
			BH2.append(sum(G2.BlackHoleMass[w2]) / G2.Mvir[w1[j]])                        
				    
		MeanCentralHaloMass2.append(np.mean(CentralHaloMass2))
		MeanBaryonFraction2.append(np.mean(BaryonFraction2))
		MeanBaryonFractionU2.append(np.mean(Ejected2) + np.var(Ejected2))
		MeanBaryonFractionL2.append(np.mean(Ejected2) - np.var(Ejected2))
		    
		MeanStars2.append(np.mean(Stars2))
		MeanCold2.append(np.mean(Cold2))
		MeanHot2.append(np.mean(Hot2))
		MeanEjected2.append(np.mean(Ejected2))
		MeanICS2.append(np.mean(ICS2))
		MeanBH.append(np.mean(BH))

	    print '  ', i, HaloRange[i], HalosFound
        
        plt.plot(MeanCentralHaloMass, MeanBaryonFraction, 'k-', label='TOTAL')#, color='purple', alpha=0.3)
        plt.fill_between(MeanCentralHaloMass, MeanBaryonFractionU, MeanBaryonFractionL, facecolor='purple', alpha=0.25)
        plt.fill_between(MeanCentralHaloMass2, MeanBaryonFractionU2, MeanBaryonFractionL2, facecolor='cyan', alpha=0.25)
        
        plt.plot(MeanCentralHaloMass, MeanStars, 'm', label='Stars')
        plt.plot(MeanCentralHaloMass, MeanCold, label='Cold', color='blue')
        plt.plot(MeanCentralHaloMass, MeanHot, label='Hot', color='red')
        plt.plot(MeanCentralHaloMass, MeanEjected, label='Ejected', color='green')
        plt.plot(MeanCentralHaloMass, MeanICS, label='ICS', color='yellow')

        plt.plot(MeanCentralHaloMass2, MeanStars2, 'm--')
        plt.plot(MeanCentralHaloMass2, MeanCold2, color = 'blue', ls = '--')
        plt.plot(MeanCentralHaloMass2, MeanHot2, color = 'red', ls = '--')
        plt.plot(MeanCentralHaloMass2, MeanEjected2, color = 'green', ls = '--')
        plt.plot(MeanCentralHaloMass2, MeanICS2, color = 'yellow', ls = '--')

	plt.plot(1e100, 1e100, label = "IRA", color = 'black', ls = '-')
	plt.plot(1e100, 1e100, label = "Delayed", color = 'black', ls = '--')
        
        plt.xlabel(r'$\mathrm{Central}\ \log_{10} M_{\mathrm{vir}}\ (M_{\odot})$')  # Set the x-axis label
        plt.ylabel(r'$\mathrm{Baryon\ Fraction}$')  # Set the y-axis label
            
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
            
        plt.axis([MinHalo - Interval, MaxHalo + Interval, 0.0, 0.23])
            
        leg = plt.legend(bbox_to_anchor=[0.99, 0.6])
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
            
        outputFile = './BaryonFrac_MiniMill' + Output_Format
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()


##

def calculate_delta_ejected():

	SnapList = [11.302, 9.457, 7.227, 5.898, 6.490, 5.000]	

	plt.figure()  # New figure
	ax = plt.subplot(111)  # 1 plot on the figure
		
	for k in xrange(0, len(SnapList)):
      
		print "Doing z = %.3f" %(SnapList[k]) 
		label = "/lustre/projects/p004_swin/jseiler/SAGE_output/1024/clean/DelayedSN/IRA_EtaSN0.5_z%.3f" %(SnapList[k])
    		G, Gal_Desc = ReadScripts.ReadGals_SAGE_NoGrid(label, 0, 124, 101)

		label = "/lustre/projects/p004_swin/jseiler/SAGE_output/1024/clean/DelayedSN/SFH4_EtaSN0.5_z%.3f" %(SnapList[k])
    		G2, Gal_Desc = ReadScripts.ReadGals_SAGE_NoGrid(label, 0, 124, 101)
    
		cut = 10000	
		if len(G) > cut:
			G = G[0:cut]
		
		HaloMass = np.log10(G.Mvir * 1.0e10 / AllVars.Hubble_h)
		Baryons = G.StellarMass + G.ColdGas + G.HotGas + G.EjectedMass + G.IntraClusterStars + G.BlackHoleMass

		HaloMass2 = np.log10(G2.Mvir * 1.0e10 / AllVars.Hubble_h)
		Baryons2 = G2.StellarMass + G2.ColdGas + G2.HotGas + G2.EjectedMass + G2.IntraClusterStars + G2.BlackHoleMass

		MinHalo = 8.0
		MaxHalo = 12.0
		Interval = 0.1
		Nbins = int((MaxHalo-MinHalo)/Interval)
		HaloRange = np.arange(MinHalo, MaxHalo, Interval)
		
		MeanCentralHaloMass = []
		MeanEjected = []
		MeanEjectedU = []
		MeanEjectedL = []

		MeanCentralHaloMass2 = []
		MeanEjected2 = []
		MeanEjectedU2 = []
		MeanEjectedL2 = []

		DeltaEjected = []
    		HaloMassPlot = []

		for i in xrange(Nbins-1):
		    
		    w1 = np.where((G.Type == 0) & (HaloMass >= HaloRange[i]) & (HaloMass < HaloRange[i+1]))[0]
		    HalosFound = len(w1)
		    

		    if HalosFound > 2:  
			
			BaryonFraction = []
			CentralHaloMass = []               
			Ejected = []
			
			for j in xrange(HalosFound):
			    
			    w2 = np.where(G.CentralGalaxyIndex == G.CentralGalaxyIndex[w1[j]])[0]
			    CentralAndSatellitesFound = len(w2)
			    
			    if CentralAndSatellitesFound > 0:
				CentralHaloMass.append(np.log10(G.Mvir[w1[j]] * 1.0e10 / AllVars.Hubble_h))
				Ejected.append(sum(G.EjectedMass[w2]) / G.Mvir[w1[j]])

					
			MeanCentralHaloMass.append(np.mean(CentralHaloMass))
			MeanEjected.append(np.mean(Ejected))
			MeanEjectedU.append(np.mean(Ejected) - np.var(Ejected))
			MeanEjectedL.append(np.mean(Ejected) + np.var(Ejected))

		    w1 = np.where((G2.Type == 0) & (HaloMass2 >= HaloRange[i]) & (HaloMass2 < HaloRange[i+1]))[0]
		    HalosFound2 = len(w1)

		    if HalosFound2 > 2:  
			    
			BaryonFraction2 = []
			CentralHaloMass2 = []
			    
			Stars2 = []
			Cold2 = []
			Hot2 = []
			Ejected2 = []
			ICS2 = []
			BH2 = []
			    
			for j in xrange(HalosFound):
				
			    w2 = np.where(G2.CentralGalaxyIndex == G2.CentralGalaxyIndex[w1[j]])[0]
			    CentralAndSatellitesFound = len(w2)
				
			    if CentralAndSatellitesFound > 0:
			
				CentralHaloMass2.append(np.log10(G2.Mvir[w1[j]] * 1.0e10 / AllVars.Hubble_h))
				Ejected2.append(sum(G2.EjectedMass[w2]) / G2.Mvir[w1[j]])

			MeanCentralHaloMass2.append(np.mean(CentralHaloMass2))
			MeanEjected2.append(np.mean(Ejected2))
			MeanEjectedU2.append(np.mean(Ejected2) - np.var(Ejected2))
			MeanEjectedL2.append(np.mean(Ejected2) + np.var(Ejected2))
				
			if(HalosFound > 2 and HalosFound2 > 2):   
				DeltaEjected.append(np.mean(Ejected) - np.mean(Ejected2))
				HaloMassPlot.append(np.mean(CentralHaloMass2))

		    print '  ', i, HaloRange[i], HalosFound

		    #print "Ejected", Ejected
		    #print "Ejected2", Ejected2
		print "HaloMassPlot", HaloMassPlot
	        print "DeltaEjected", DeltaEjected

		label = "z = %.2f" %(SnapList[k])
        	plt.plot(HaloMassPlot, DeltaEjected, color = colors[k], label = label)
	


        plt.xlabel(r'$\mathrm{Central}\ \log_{10} M_{\mathrm{vir}}\ (M_{\odot})$')  # Set the x-axis label
        #plt.ylabel(r'$\mathrm{Ejected Baryon Fraction, IRA - Delayed SN}$')  # Set the y-axis label
        plt.ylabel(r'$\mathrm{Ejected Baryon Fraction, Instant - Delayed SN}$')  # Set the y-axis label
            
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
           
        plt.axis([MinHalo - Interval, MaxHalo + Interval, -0.05, 0.05])
            
        leg = plt.legend(bbox_to_anchor=[0.99, 0.6])
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
            
        outputFile = './BaryonFraction_Comp_Instant_Delayed' + Output_Format
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()


##

def BlackHoleBulge():

      	AllVars.Set_Params_MiniMill()

	dilute = 1000    
    	plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure
   
	label = "/lustre/projects/p004_swin/jseiler/mini_millennium/Delayed_SN1Myr_z0.000"
    	G, Gal_Desc = ReadScripts.ReadGals_SAGE_DelayedSN(label, 0, 7, 64)

        w = np.where((G.BulgeMass > 0.01) & (G.BlackHoleMass > 0.00001))[0]
        if(len(w) > dilute): w = sample(w, dilute)
    
        bh = np.log10(G.BlackHoleMass[w] * 1.0e10 / AllVars.Hubble_h)
        bulge = np.log10(G.BulgeMass[w] * 1.0e10 / AllVars.Hubble_h)
                    
        plt.scatter(bulge, bh, marker='o', s=1, c='k', alpha=0.5, label='Model galaxies')
                
        # overplot Haring & Rix 2004
        w = 10. ** np.arange(20)
        BHdata = 10. ** (8.2 + 1.12 * np.log10(w / 1.0e11))
        plt.plot(np.log10(w), np.log10(BHdata), 'b-', label="Haring \& Rix 2004")

        plt.ylabel(r'$\log\ M_{\mathrm{BH}}\ (M_{\odot})$')  # Set the y...
        plt.xlabel(r'$\log\ M_{\mathrm{bulge}}\ (M_{\odot})$')  # and the x-axis labels
            
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
            
        plt.axis([8.0, 12.0, 6.0, 10.0])
            
        leg = plt.legend(loc='upper left')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')
            
        outputFile = './BlackHoleBulge_Delayed' + Output_Format 
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()


##

def GasMassFunction():

	sSFRcut = -11.0

      	AllVars.Set_Params_MiniMill()
        plt.figure()  # New figure
        ax = plt.subplot(111)  # 1 plot on the figure

        binwidth = 0.1  # mass function histogram bin width

	label = "/lustre/projects/p004_swin/jseiler/mini_millennium/IRA_z0.000"
    	G, Gal_Desc = ReadScripts.ReadGals_SAGE_DelayedSN(label, 0, 7, 64)

        # calculate all
        w = np.where(G.ColdGas > 0.0)[0]
        mass = np.log10(G.ColdGas[w] * 1.0e10 / AllVars.Hubble_h)
        sSFR = (G.SfrDisk[w] + G.SfrBulge[w]) / (G.StellarMass[w] * 1.0e10 / AllVars.Hubble_h)
        mi = np.floor(min(mass)) - 2
        ma = np.floor(max(mass)) + 2
        NB = (ma - mi) / binwidth

        (counts, binedges) = np.histogram(mass, range=(mi, ma), bins=NB)

        # Set the x-axis values to be the centre of the bins
        xaxeshisto = binedges[:-1] + 0.5 * binwidth
        
        # additionally calculate red
        w = np.where(sSFR < 10.0**sSFRcut)[0]
        massRED = mass[w]
        (countsRED, binedges) = np.histogram(massRED, range=(mi, ma), bins=NB)

        # additionally calculate blue
        w = np.where(sSFR > 10.0**sSFRcut)[0]
        massBLU = mass[w]
        (countsBLU, binedges) = np.histogram(massBLU, range=(mi, ma), bins=NB)

        # Baldry+ 2008 modified data used for the MCMC fitting
        Zwaan = np.array([[6.933,   -0.333],
            [7.057,   -0.490],
            [7.209,   -0.698],
            [7.365,   -0.667],
            [7.528,   -0.823],
            [7.647,   -0.958],
            [7.809,   -0.917],
            [7.971,   -0.948],
            [8.112,   -0.927],
            [8.263,   -0.917],
            [8.404,   -1.062],
            [8.566,   -1.177],
            [8.707,   -1.177],
            [8.853,   -1.312],
            [9.010,   -1.344],
            [9.161,   -1.448],
            [9.302,   -1.604],
            [9.448,   -1.792],
            [9.599,   -2.021],
            [9.740,   -2.406],
            [9.897,   -2.615],
            [10.053,  -3.031],
            [10.178,  -3.677],
            [10.335,  -4.448],
            [10.492,  -5.083]        ], dtype=np.float32)
        
        ObrRaw = np.array([
            [7.300,   -1.104],
            [7.576,   -1.302],
            [7.847,   -1.250],
            [8.133,   -1.240],
            [8.409,   -1.344],
            [8.691,   -1.479],
            [8.956,   -1.792],
            [9.231,   -2.271],
            [9.507,   -3.198],
            [9.788,   -5.062 ]        ], dtype=np.float32)

        ObrCold = np.array([
            [8.009,   -1.042],
            [8.215,   -1.156],
            [8.409,   -0.990],
            [8.604,   -1.156],
            [8.799,   -1.208],
            [9.020,   -1.333],
            [9.194,   -1.385],
            [9.404,   -1.552],
            [9.599,   -1.677],
            [9.788,   -1.812],
            [9.999,   -2.312],
            [10.172,  -2.656],
            [10.362,  -3.500],
            [10.551,  -3.635],
            [10.740,  -5.010]        ], dtype=np.float32)

        ObrCold_xval = np.log10(10**(ObrCold[:, 0])  /AllVars.Hubble_h/AllVars.Hubble_h)
        ObrCold_yval = (10**(ObrCold[:, 1]) * AllVars.Hubble_h*AllVars.Hubble_h*AllVars.Hubble_h)
        Zwaan_xval = np.log10(10**(Zwaan[:, 0]) /AllVars.Hubble_h/AllVars.Hubble_h)
        Zwaan_yval = (10**(Zwaan[:, 1]) * AllVars.Hubble_h*AllVars.Hubble_h*AllVars.Hubble_h)
        ObrRaw_xval = np.log10(10**(ObrRaw[:, 0])  /AllVars.Hubble_h/AllVars.Hubble_h)
        ObrRaw_yval = (10**(ObrRaw[:, 1]) * AllVars.Hubble_h*AllVars.Hubble_h*AllVars.Hubble_h)

        plt.plot(ObrCold_xval, ObrCold_yval, color='black', lw = 7, alpha=0.25, label='Obr. \& Raw. 2009 (Cold Gas)')
        plt.plot(Zwaan_xval, Zwaan_yval, color='cyan', lw = 7, alpha=0.25, label='Zwaan et al. 2005 (HI)')
        plt.plot(ObrRaw_xval, ObrRaw_yval, color='magenta', lw = 7, alpha=0.25, label='Obr. \& Raw. 2009 (H2)')

        
        # Overplot the model histograms
        plt.plot(xaxeshisto, counts    / AllVars.Volume * AllVars.Hubble_h*AllVars.Hubble_h*AllVars.Hubble_h / binwidth, 'k-', label='Model - Cold Gas')

        plt.yscale('log', nonposy='clip')
        plt.axis([8.0, 11.5, 1.0e-6, 1.0e-1])

        # Set the x-axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

        plt.ylabel(r'$\phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$')  # Set the y...
        plt.xlabel(r'$\log_{10} M_{\mathrm{X}}\ (M_{\odot})$')  # and the x-axis labels

        leg = plt.legend(loc='lower left', numpoints=1,
                         labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')

        outputFile = './GasMassFunction_IRA' + Output_Format 
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()


##


def OutflowRate(SnapList, mass1, mass2, outflow1, outflow2, Model_Tags, Output_Tag):


    ax = plt.subplot(111)  # 1 plot on the figure

    binwidth = 0.1
    low_mass = 8.0 
    high_mass = 12.0 

    bins = np.arange(low_mass,high_mass, binwidth)
    bins_mid = bins + binwidth/2.0
    bins_mid = bins_mid[:-1] # len(bins_mid) should be 1 less than len(bins) as the last bin doesn't have a midpoint.

    for i in xrange(0, len(SnapList)):
	Outflow_Sum = []
	Outflow_Sum2 = []	
	for j in xrange(0, len(bins)-1):	
		w = np.where((mass1[i] >= bins[j]) & (mass1[i] < bins[j+1]))[0]

		if (len(w) != 0):
			Outflow_Sum.append(np.sum(outflow1[i][w])/len(w))
		else:
			Outflow_Sum.append(nan)

		w = np.where((mass2[i] >= bins[j]) & (mass2[i] < bins[j+1]))[0]
		if (len(w) != 0):
			Outflow_Sum2.append(np.sum(outflow2[i][w])/len(w))
		else:
			Outflow_Sum2.append(nan)
	

	print Outflow_Sum
	print Outflow_Sum2
#    	ax.plot(bins_mid, Outflow_Sum, color = colors[i], ls = '-')
#    	ax.plot(bins_mid, Outflow_Sum2, color = colors[i], ls = '--')

	tmp_percentage = np.divide(np.subtract(Outflow_Sum,Outflow_Sum2), Outflow_Sum)
	tmp_percentage[tmp_percentage < -50] = nan
	ax.plot(bins_mid, tmp_percentage, color = colors[i], ls = '-')
	
	label = 'z = %.2f' %(AllVars.SnapZ[SnapList[i]])
	ax.plot(nan, nan, color = colors[i], label = label)
 
#    ax.plot(nan, nan, 'k-', label = Model_Tags[0])
#    ax.plot(nan, nan, 'k--', label = Model_Tags[1])
    #ax.set_xlabel(r'$\log_{10}\ M_{\mathrm{Vir}} \:[M_{\odot}]$', fontsize = talk_fontsize)
    ax.set_xlabel(r'$\log_{10}\ M_{\mathrm{*}} \:[M_{\odot}]$', fontsize = talk_fontsize)
    #ax.set_ylabel(r'$\mathrm{Outflow \: \: Rate} \: [\mathrm{M_\odot}\mathrm{yr}^{-1}]$')
    ax.set_ylabel(r'$\mathrm{Delta Outflow \: \: Rate, IRA - Delayed} \: [\mathrm{M_\odot}\mathrm{yr}^{-1}]$')
    #ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

    #ax.set_yscale('log', nonposy='clip')
    ax.set_xlim([low_mass - binwidth, high_mass + binwidth])
#    ax.set_ylim([0,5])

    leg = plt.legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
    	t.set_fontsize('medium')


    outputFile = Output_Tag + Output_Format 
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

##

#################################

Simulation = 0 # Set 0 for Mini-Millennium, 1 for My_Simulation, 2 for both (kinda).

      
if (Simulation == 0 or Simulation == 2):
    H_Millennium = ReadScripts.ReadHalos('/lustre/projects/p004_swin/jseiler/millennium_trees/trees_063', 0, 7)

    GG_Millennium, Gal_Desc = ReadScripts.ReadGals_SAGE_DelayedSN('/lustre/projects/p004_swin/jseiler/mini_millennium/IRA_z0.000', 0, 7, 64)
    G_Merged_Millennium, Merged_Desc = ReadScripts.ReadGals_SAGE_DelayedSN('/lustre/projects/p004_swin/jseiler/mini_millennium/IRA_MergedGalaxies', 0, 7, 64)
    G_Millennium = ReadScripts.Join_Arrays(GG_Millennium, G_Merged_Millennium, Gal_Desc)

    GG_Millennium2, Gal_Desc = ReadScripts.ReadGals_SAGE_DelayedSN('/lustre/projects/p004_swin/jseiler/mini_millennium/Delayed_SN1Myr_z0.000', 0, 7, 64)
    G_Merged_Millennium2, Gal_Desc = ReadScripts.ReadGals_SAGE_DelayedSN('/lustre/projects/p004_swin/jseiler/mini_millennium/Delayed_SN1Myr_MergedGalaxies', 0, 7, 64)
    G_Millennium2 = ReadScripts.Join_Arrays(GG_Millennium2, G_Merged_Millennium2, Gal_Desc)

    SnapList_Millennium = [15, 25, 30, 45, 58, 63] 
    print "Snapshots analyzing are", SnapList_Millennium
    AllVars.Set_Params_MiniMill()

if (Simulation == 1 or Simulation == 2):

    #H_MySim = ReadScripts.ReadHalos('/lustre/projects/p004_swin/jseiler/Rockstar_output/Halos_noLL/Ltrees/lhalotree.bin', 0, 26)
    #H_MySim2 = ReadScripts.ReadHalos('/lustre/projects/p004_swin/jseiler/Rockstar_output/1024_Halos_noLL/Ltrees/lhalotree.bin', 0, 124)
    #H_MySim2 = ReadScripts.ReadHalos('/lustre/projects/p004_swin/jseiler/Rockstar_output/Halos_noLL/Ltrees/lhalotree.bin', 0, 26)

   
    #print "Minimum halo mass is %.4e and maximum is %.4e (log Msun)" %(np.log10(np.amin(H_MySim2.Mvir) * 1.0e10/AllVars.Hubble_h), np.log10(np.amax(H_MySim2.Mvir) * 1.0e10/AllVars.Hubble_h))

    # 512 
    #GG_MySim, Gal_Desc = ReadScripts.ReadGals_SAGE_Photons('/lustre/projects/p004_swin/jseiler/SAGE_output/512/fiducial_z5.038', 0, 26, 100)
    #G_Merged_MySim, Merged_Desc = ReadScripts.ReadGals_SAGE_Photons('/lustre/projects/p004_swin/jseiler/SAGE_output/512/fiducial_MergedGalaxies', 0, 26, 100)
    #G_MySim = ReadScripts.Join_Arrays(GG_MySim, G_Merged_MySim, Gal_Desc)

    # 1024

    GG_MySim, Gal_Desc = ReadScripts.ReadGals_SAGE_DelayedSN('/lustre/projects/p004_swin/jseiler/SAGE_output/1024/clean/DelayedSN/NewDelayed_Delayed_SN1Myr_SF1Myr_z5.000', 0, 124, 101)
    G_Merged_MySim, Merged_Desc = ReadScripts.ReadGals_SAGE_DelayedSN('/lustre/projects/p004_swin/jseiler/SAGE_output/1024/clean/DelayedSN/NewDelayed_Delayed_SN1Myr_SF1Myr_MergedGalaxies', 0, 124, 101)

    #GG_MySim, Gal_Desc = ReadScripts.ReadGals_SAGE_Photons('/lustre/projects/p004_swin/jseiler/SAGE_output/1024/clean/LenHistory/LenHistory_SF0.01_noreion_z5.000', 0, 124, 101)
    #G_Merged_MySim, Merged_Desc = ReadScripts.ReadGals_SAGE_Photons('/lustre/projects/p004_swin/jseiler/SAGE_output/1024/clean/LenHistory/LenHistory_SF0.01_noreion_MergedGalaxies', 0, 124, 101)

    G_MySim = ReadScripts.Join_Arrays(GG_MySim, G_Merged_MySim, Gal_Desc)


    # 512
    #GG_MySim2, Gal_Desc = ReadScripts.ReadGals_SAGE_Photons('/lustre/projects/p004_swin/jseiler/SAGE_output/512/noreion_z5.038', 0, 26, 100)
    #G_Merged_MySim2, Merged_Desc = ReadScripts.ReadGals_SAGE_Photons('/lustre/projects/p004_swin/jseiler/SAGE_output/512/noreion_z5.038', 0, 26, 100)
    #G_MySim2 = ReadScripts.Join_Arrays(GG_MySim2, G_Merged_MySim2, Gal_Desc)

    #GG_MySim2, Gal_Desc = ReadScripts.ReadGals_SAGE_Photons('/lustre/projects/p004_swin/jseiler/SAGE_output/tmp_z5.038', 0, 26, 100)
    #G_Merged_MySim2, Merged_Desc = ReadScripts.ReadGals_SAGE_Photons('/lustre/projects/p004_swin/jseiler/SAGE_output/512/tmp_z5.038', 0, 26, 100)
    #G_MySim2 = ReadScripts.Join_Arrays(GG_MySim2, G_Merged_MySim2, Gal_Desc)


    # 1024
    GG_MySim2, Gal_Desc = ReadScripts.ReadGals_SAGE_DelayedSN('/lustre/projects/p004_swin/jseiler/SAGE_output/1024/clean/DelayedSN/NewDelayed_Delayed_SN10Myr_SF10Myr_z5.000', 0, 124, 101)
    G_Merged_MySim2, Merged_Desc = ReadScripts.ReadGals_SAGE_DelayedSN('/lustre/projects/p004_swin/jseiler/SAGE_output/1024/clean/DelayedSN/NewDelayed_Delayed_SN10Myr_SF10Myr_MergedGalaxies', 0, 124, 101)

    #GG_MySim2, Gal_Desc = ReadScripts.ReadGals_SAGE_Photons('/lustre/projects/p004_swin/jseiler/SAGE_output/1024/clean/LenHistory/LenHistory_SF0.01_noreion_z5.000', 0, 124, 101)
    #G_Merged_MySim2, Merged_Desc = ReadScripts.ReadGals_SAGE_Photons('/lustre/projects/p004_swin/jseiler/SAGE_output/1024/clean/LenHistory/LenHistory_SF0.01_noreion_MergedGalaxies', 0, 124, 101)

    G_MySim2 = ReadScripts.Join_Arrays(GG_MySim2, G_Merged_MySim2, Gal_Desc) 



#    SnapList_MySim = np.arange(22,78) 


   # SnapList_MySim = [30, 35, 40, 50, 60, 70, 78]
    #SnapList_MySim = [99, 78, 64, 51, 43, 37] 
    # Snapshots = z [5, 6, 7, 8, 9, 10] are [99, 78, 64, 51, 43, 37] (used for UVLF) 

    SnapList_MySim = [78, 64,51]
    # Snapshots for z = [6, 7, 8] are [78, 64, 51]
    # Snapshots for z = [7.23, 6.69] are [60, 67] (used for LyAlpha Luminosity)
 
    print "Snapshots analyzing are", SnapList_MySim
    
HaloPart_Low = 41 # Bounds for where we define the cutoff for a 'Dark Matter Halo'. Below this we can't be sure of the results.
HaloPart_High = 51

calculate_observed_LF = 0



'''
countSnap_low = 30
countSnap_high = 99

idx_low = [i for i in range(len(G_MySim.GridHistory)) if G_MySim.GridHistory[i,countSnap_low] != -1] # Indices for all galaxies that exist at snapshot 'countSnap_low'.
numgals_low = len(idx_low)

numgals_high = len(G_MySim.GridHistory[G_MySim.GridHistory[idx_low, countSnap_high] != -1])
 

print numgals_low
print numgals_high
'''


#HaloPartCount(Simulation, SnapListZ_MySim, HaloCount_MySim, len(SnapList_MySim))
## Halo Initialization ##

w_H_MySim = []
mass_H_MySim = []
Photons_HI_H_MySim = []
Photons_HI_Tot_H_MySim = []
SourceEfficiency_H_MySim = 10

w_H_MySim2 = []
mass_H_MySim2 = []
Photons_HI_H_MySim2 = []
Photons_HI_Tot_H_MySim2 = []
SourceEfficiency_H_MySim2 = 10


## MySim Initialization ## 

w_G_MySim = []
mass_G_MySim = []
Photons_HI_G_MySim = []
Photons_HI_Tot_G_MySim = [] 
SFR_G_MySim = []
sSFR_G_MySim = []
sSFR_min_MySim = 1e100
sSFR_max_MySim = -1e100
Metallicity_Tremonti_G_MySim = []
HaloPart_MySim = []
fesc_MySim = 0.15

mass_Central_MySim = []
Photons_HI_Central_MySim = []
Photons_HI_Tot_Central_MySim = []
SourceEfficiency_MySim = 1 

SnapListZ = []
SnapListZ_MySim = []

w_Ionized_MySim = []
MfiltGnedin_MySim = []
MfiltSobacchi_MySim = []

LymanAlpha_MySim = []
fesc_LymanAlpha_MySim = 0.3

LUV_MySim = []
MUV_MySim = []

MUV_Obs_MySim = []
mean_A_MySim = []

EjectedFraction_MySim = []

fesc_local_MySim = []
fesc_local2_MySim = []
fesc_local3_MySim = []

HaloCount_MySim = []
HaloCut_MySim = 1

## MySim Model 2 ##

w_G_MySim2 = []
mass_G_MySim2 = []
Photons_HI_G_MySim2 = []
Photons_HI_Tot_G_MySim2 = []
SFR_G_MySim2 = []
sSFR_G_MySim2 = []
Metallicity_Tremonti_G_MySim2 = []
HaloPart_MySim2 = []
fesc_MySim2 = 0.15

mass_Central_MySim2 = []
Photons_HI_Central_MySim2 = []
Photons_HI_Tot_Central_MySim2 = []
SourceEfficiency_MySim2 = 1

w_Ionized_MySim2 = []
MfiltGnedin_MySim2 = []
MfiltSobacchi_MySim2 = []

LymanAlpha_MySim2 = []
fesc_LymanAlpha_MySim2 = 0.3 

LUV_MySim2 = []
MUV_MySim2 = []

MUV_Obs_MySim2 = []
mean_A_MySim2 = []

fesc_local_MySim2 = []

HaloCut_MySim2 = 1

## Millennium Initialization ##

w_H_Millennium = []
mass_H_Millennium = []
Photons_HI_H_Millennium = []
Photons_HI_Tot_H_Millennium = []
SourceEfficiency_H_Millennium = 50

w_G_Millennium = []
mass_G_Millennium = []
Photons_HI_G_Millennium = []
Photons_HI_Tot_G_Millennium = []
SFR_G_Millennium = []
sSFR_G_Millennium = []
Metallicity_Tremonti_G_Millennium = []
fesc_Millennium = 0.675

mass_Central_Millennium = []
Photons_HI_Central_Millennium = []
Photons_HI_Tot_Central_Millennium = []
SourceEfficiency_Millennium = 50

SnapListZ_Millennium = []

mass_Central_Millennium2 = []
Outflow_G_Millennium = []

## Millennium Model 2 Initialization ##

w_G_Millennium2 = []
mass_G_Millennium2 = []
mass_Central_Millennium2 = []

## Model 2 Initialization ##


w_H_Millennium2 = []
mass_H_Millennium2 = []
Photons_HI_H_Millennium2 = []
Photons_HI_Tot_H_Millennium2 = []
SourceEfficiency_H_Millennium2 = 10

w_G_Millennium2 = []
mass_G_Millennium2 = []
Photons_HI_G_Millennium2 = []
Photons_HI_Tot_G_Millennium2 = []
SFR_G_Millennium2 = []
sSFR_G_Millennium2 = []
Metallicity_Tremonti_G_Millennium2 = []
fesc_Millennium2 = 0.675

mass_Central_Millennium2 = []
Photons_HI_Central_Millennium2 = []
Photons_HI_Tot_Central_Millennium2 = []
SourceEfficiency_Millennium2 = 10


Outflow_G_Millennium2 = []

##

if (Simulation == 1 or Simulation == 2):
    for i in xrange(0, len(SnapList_MySim)):
 
    ###### MYSIM CALCULATIONS #####

      AllVars.Set_Params_Mysim()

      SnapListZ_MySim.append(AllVars.SnapZ[SnapList_MySim[i]])
      SnapListZ.append(AllVars.SnapZ[SnapList_MySim[i]])

      ## Halo Calculations ##

      #w_H_MySim.append(np.where((H_MySim.SnapNum == SnapList_MySim[i]) & (H_MySim.Mvir > 0.0))[0])  
      #mass_H_MySim.append(np.log10(H_MySim.Mvir[w_H_MySim[i]] * 1.0e10 / AllVars.Hubble_h))
      #Photon_Factor = AllVars.Solar_Mass * SourceEfficiency_H_MySim * AllVars.BaryonFrac / AllVars.Ionized_Mass_H / AllVars.Proton_Mass / (AllVars.Lookback_Time[SnapList_MySim[i]] - AllVars.Lookback_Time[SnapList_MySim[i]+1]) / AllVars.Sec_Per_Megayear / 1e3
      #Photons_HI_H_MySim.append(np.log10(H_MySim.Mvir[w_H_MySim[i]]* 1.0e10 / AllVars.Hubble_h * Photon_Factor))
      #Photons_HI_Tot_H_MySim.append(np.log10(sum(10**Photons_HI_H_MySim[i])))

      ## Model 1 Calculations ##

      #w_G_MySim.append(np.where((G_MySim.GridHistory[:, SnapList_MySim[i]] != -1) & (G_MySim.GridStellarMass[:, SnapList_MySim[i]] > 0.0) & (G_MySim.GridSFR[:, SnapList_MySim[i]] > 0.0) & (G_MySim.GridCentralGalaxyMass[:,SnapList_MySim[i]] > 0.0) & (G_MySim.GridZ[:, SnapList_MySim[i]] > 0.0))[0])
      w_G_MySim.append(np.where((G_MySim.GridHistory[:, SnapList_MySim[i]] != -1) & (G_MySim.GridStellarMass[:, SnapList_MySim[i]] > 0.0) & (G_MySim.GridSFR[:, SnapList_MySim[i]] > 0.0) & (G_MySim.LenHistory[:, SnapList_MySim[i]] > HaloCut_MySim))[0])

      HaloCount_MySim.append(G_MySim.LenHistory[w_G_MySim[i], SnapList_MySim[i]])
      mass_G_MySim.append(np.log10(G_MySim.GridStellarMass[w_G_MySim[i], SnapList_MySim[i]] * 1.0e10 / AllVars.Hubble_h))
#      Photons_HI_G_MySim.append(G_MySim.Photons_HI[w_G_MySim[i], SnapList_MySim[i]] + np.log10(fesc_MySim))
    
#      Photons_HI_Tot_G_MySim.append(np.log10(Sum_Log(Photons_HI_G_MySim[i])))
 
      SFR_G_MySim.append(np.log10(G_MySim.GridSFR[w_G_MySim[i], SnapList_MySim[i]])) # Msun yr^-1.  Log Units.

      HaloPart_MySim.append(G_MySim.Len[w_G_MySim[i]])
      
      Metallicity_Tremonti_G_MySim.append(np.log10(G_MySim.GridZ[w_G_MySim[i], SnapList_MySim[i]] / 0.02) + 9.0)

      print "There were %d galaxies for snapshot %d (Redshift %.4f) for model1." %(len(w_G_MySim[i]), SnapList_MySim[i], AllVars.SnapZ[SnapList_MySim[i]])

      mass_Central_MySim.append(np.log10(G_MySim.GridCentralGalaxyMass[w_G_MySim[i], SnapList_MySim[i]] * 1.0e10 / AllVars.Hubble_h))
      Photon_Factor = AllVars.Solar_Mass * SourceEfficiency_MySim * AllVars.BaryonFrac / AllVars.Ionized_Mass_H / AllVars.Proton_Mass / (AllVars.Lookback_Time[SnapList_MySim[i]] - AllVars.Lookback_Time[SnapList_MySim[i]+1]) / AllVars.Sec_Per_Megayear / 1e3
      Photons_HI_Central_MySim.append(np.log10(10**(mass_Central_MySim[i]) * Photon_Factor))
      Photons_HI_Tot_Central_MySim.append(np.log10(sum(10**Photons_HI_Central_MySim[i])))

      #w_Ionized_MySim.append(np.where((G_MySim.GridHistory[:, SnapList_MySim[i]] != -1) & (G_MySim.GridStellarMass[:, SnapList_MySim[i]] > 0.0) & (G_MySim.GridSFR[:, SnapList_MySim[i]] > 0.0) & (G_MySim.Photons_HI[:, SnapList_MySim[i]] > 0.0) & (G_MySim.MfiltSobacchi[:, SnapList_MySim[i]] > 1.0))[0])

      #MfiltGnedin_MySim.append(G_MySim.MfiltGnedin[w_G_MySim[i], SnapList_MySim[i]])
      #MfiltSobacchi_MySim.append(G_MySim.MfiltSobacchi[w_Ionized_MySim[i], SnapList_MySim[i]])
      #MfiltSobacchi_MySim.append(G_MySim.MfiltSobacchi[w_G_MySim[i], SnapList_MySim[i]])

#      LymanAlpha_MySim.append(np.log10(0.68*(1.0-fesc_MySim)*fesc_LymanAlpha_MySim * (AllVars.LymanAlpha_Energy* AllVars.eV_to_erg)) + Photons_HI_G_MySim[i]) 

      print "Calculating intrinsic Mag"
      LUV_MySim.append((SFR_G_MySim[i] + 39.927)) # Using relationship from STARBURST99, units of erg s^-1 A^-1. Log Units.
      MUV_MySim.append(AllVars.Luminosity_to_ABMag(LUV_MySim[i], 1600))	

      if (calculate_observed_LF == 1):
	      print "Calculating dust"
	      MUV_bins = np.arange(-24, -16, 0.1)
	      A_Mean = np.zeros((len(MUV_bins)))

	      for j in xrange(0, len(MUV_bins)):
		beta = calculate_beta(MUV_bins[j], SnapListZ[i]) 
		dist = np.random.normal(beta, 0.34, 10000)
		A = 4.43 + 1.99*dist
		A[A < 0] = 0
		
		A_Mean[j] = np.mean(A)
		
	      print "Binning the intrinsic into Observed."
	 
	      indices = np.digitize(MUV_MySim[i], MUV_bins) - 1 # Bins the simulation magnitude into the MUV bins. Note that digitize defines an index i if bin[i-1] <= x < bin[i] whereas I prefer bin[i] <= x < bin[i+1] 

	      dust = A_Mean[indices]

	      Flux = AllVars.Luminosity_to_Flux(LUV_MySim[i], 10.0) # Calculate the flux from a distance of 10 parsec, units of erg s^-1 A^-1 cm^-2.  Log Units. 
	      Flux_Observed = Flux - 0.4*dust	
		
	      f_nu = AllVars.spectralflux_wavelength_to_frequency(10**Flux_Observed, 1600) # Spectral flux density in Janksy.
	      MUV_Obs_MySim.append(-2.5 * np.log10(f_nu) + 8.90) # AB Magnitude from http://www.astro.ljmu.ac.uk/~ikb/convert-units/node2.html

      EjectedFraction_MySim.append(G_MySim.EjectedFraction[w_G_MySim[i], SnapList_MySim[i]])

      fesc_local_MySim.append(0.25)
      fesc_local2_MySim.append(pow(10,3.62) * pow(pow(10,mass_Central_MySim[i]),-0.43)) 
      fesc_local3_MySim.append(EjectedFraction_MySim[i]*0.999 + 0.001)
	
      ## Model 2 Calculations ##

      ## Halos ##

#      w_H_MySim2.append(np.where((H_MySim2.SnapNum == SnapList_MySim[i]) & (H_MySim2.Mvir > 0.0))[0]) 
#      mass_H_MySim2.append(np.log10(H_MySim2.Mvir[w_H_MySim2[i]] * 1.0e10 / AllVars.Hubble_h))
#      Photon_Factor = AllVars.Solar_Mass * SourceEfficiency_H_MySim * AllVars.BaryonFrac / AllVars.Ionized_Mass_H / AllVars.Proton_Mass / (AllVars.Lookback_Time[SnapList_MySim[i]] - AllVars.Lookback_Time[SnapList_MySim[i]+1]) / AllVars.Sec_Per_Megayear / 1e3
#      Photons_HI_H_MySim2.append(np.log10(H_MySim2.Mvir[w_H_MySim2[i]]* 1.0e10 / AllVars.Hubble_h * Photon_Factor))
#      Photons_HI_Tot_H_MySim2.append(np.log10(sum(10**Photons_HI_H_MySim2[i])))

      ## Galaxies ##
 
      w_G_MySim2.append(np.where((G_MySim2.GridHistory[:, SnapList_MySim[i]] != -1) & (G_MySim2.GridStellarMass[:, SnapList_MySim[i]] > 0.0) & (G_MySim2.GridSFR[:, SnapList_MySim[i]] > 0.0) & (G_MySim2.LenHistory[:, SnapList_MySim[i]] > HaloCut_MySim))[0])

      #w_G_MySim2.append(np.where((G_MySim2.GridHistory[:, SnapList_MySim[i]] != -1) & (G_MySim2.GridStellarMass[:, SnapList_MySim[i]] > 0.0) & (G_MySim2.GridSFR[:, SnapList_MySim[i]] > 0.0) & (G_MySim2.GridCentralGalaxyMass[:,SnapList_MySim[i]] > 0.0) & (G_MySim2.GridZ[:, SnapList_MySim[i]] > 0.0))[0])
      mass_G_MySim2.append(np.log10(G_MySim2.GridStellarMass[w_G_MySim2[i], SnapList_MySim[i]] * 1.0e10 / AllVars.Hubble_h))
      #Photons_HI_G_MySim2.append(G_MySim2.Photons_HI[w_G_MySim2[i], SnapList_MySim[i]] + np.log10(fesc_MySim2)) 
      #Photons_HI_Tot_G_MySim2.append(np.log10(Sum_Log(Photons_HI_G_MySim2[i])))

      SFR_G_MySim2.append(np.log10(G_MySim2.GridSFR[w_G_MySim2[i], SnapList_MySim[i]])) 
      
      sSFR_G_MySim2.append(SFR_G_MySim2[i] - mass_G_MySim2[i])
      HaloPart_MySim2.append(G_MySim2.Len[w_G_MySim2[i]])
      Metallicity_Tremonti_G_MySim2.append(np.log10(G_MySim2.GridZ[w_G_MySim2[i], SnapList_MySim[i]] / 0.02) + 9.0)

      mass_Central_MySim2.append(np.log10(G_MySim2.GridCentralGalaxyMass[w_G_MySim2[i], SnapList_MySim[i]] * 1.0e10 / AllVars.Hubble_h))
      Photon_Factor = AllVars.Solar_Mass * SourceEfficiency_MySim2 * AllVars.BaryonFrac / AllVars.Ionized_Mass_H / AllVars.Proton_Mass / (AllVars.Lookback_Time[SnapList_MySim[i]] - AllVars.Lookback_Time[SnapList_MySim[i]+1]) / AllVars.Sec_Per_Megayear / 1e3
#      Photons_HI_Central_MySim2.append(np.log10(10**(mass_Central_MySim2[i]) * Photon_Factor))
#      Photons_HI_Tot_Central_MySim2.append(np.log10(sum(10**Photons_HI_Central_MySim2[i]))) 

      print "There were %d galaxies for snapshot %d (Redshift %.4f) for model2." %(len(w_G_MySim2[i]), SnapList_MySim[i], AllVars.SnapZ[SnapList_MySim[i]])

#      w_Ionized_MySim2.append(np.where((G_MySim2.GridHistory[:, SnapList_MySim[i]] != -1) & (G_MySim2.GridStellarMass[:, SnapList_MySim[i]] > 0.0) & (G_MySim2.GridSFR[:, SnapList_MySim[i]] > 0.0) & (G_MySim2.Photons_HI[:, SnapList_MySim[i]] > 0.0) & (G_MySim2.MfiltSobacchi[:, SnapList_MySim[i]] < 1.0))[0])

#      MfiltGnedin_MySim2.append(G_MySim2.MfiltGnedin[w_G_MySim2[i], SnapList_MySim[i]])
#      MfiltSobacchi_MySim2.append(G_MySim2.MfiltSobacchi[w_G_MySim2[i], SnapList_MySim[i]])

      #LymanAlpha_MySim2.append(np.log10(0.68*(1.0-fesc_MySim2)*fesc_LymanAlpha_MySim2 * (AllVars.LymanAlpha_Energy* AllVars.eV_to_erg)) + Photons_HI_G_MySim2[i]) 

      #fesc_local_MySim2.append(10**(1.00 - 0.2 * mass_Central_MySim2[i]))

      print "Calculating intrinsic Mag"
      LUV_MySim2.append((SFR_G_MySim2[i] + 39.927)) # Using relationship from STARBURST99, units of erg s^-1 A^-1. Log Units.
      MUV_MySim2.append(AllVars.Luminosity_to_ABMag(LUV_MySim2[i], 1600))	

      if (calculate_observed_LF == 1):
	      print "Calculating dust"
	      MUV_bins = np.arange(-24, -16, 0.1)
	      A_Mean = np.zeros((len(MUV_bins)))

	      for j in xrange(0, len(MUV_bins)):
		beta = calculate_beta(MUV_bins[j], SnapListZ[i]) 
		dist = np.random.normal(beta, 0.34, 10000)
		A = 4.43 + 1.99*dist
		A[A < 0] = 0
		
		A_Mean[j] = np.mean(A)
		
	      print "Binning the intrinsic into Observed."
	 
	      indices = np.digitize(MUV_MySim2[i], MUV_bins) - 1 # Bins the simulation magnitude into the MUV bins. Note that digitize defines an index i if bin[i-1] <= x < bin[i] whereas I prefer bin[i] <= x < bin[i+1] 

	      dust = A_Mean[indices]

	      Flux = AllVars.Luminosity_to_Flux(LUV_MySim2[i], 10.0) # Calculate the flux from a distance of 10 parsec, units of erg s^-1 A^-1 cm^-2.  Log Units. 
	      Flux_Observed = Flux - 0.4*dust	
		
	      f_nu = AllVars.spectralflux_wavelength_to_frequency(10**Flux_Observed, 1600) # Spectral flux density in Janksy.
	      MUV_Obs_MySim2.append(-2.5 * np.log10(f_nu) + 8.90) # AB Magnitude from http://www.astro.ljmu.ac.uk/~ikb/convert-units/node2.html


if (Simulation == 0 or Simulation == 2):
    for i in xrange(0, len(SnapList_Millennium)):

    ##### MILLENNIUM CALCULATIONS ####
    
   
      AllVars.Set_Params_MiniMill()
      ## Halos ##

      w_H_Millennium.append(np.where((H_Millennium.SnapNum == SnapList_Millennium[i]) & (H_Millennium.Mvir > 0.0))[0]) 
      mass_H_Millennium.append(np.log10(H_Millennium.Mvir[w_H_Millennium[i]] * 1.0e10 / AllVars.Hubble_h))
      #Photon_Factor = AllVars.Solar_Mass * SourceEfficiency_H_Millennium * AllVars.BaryonFrac / AllVars.Ionized_Mass_H / AllVars.Proton_Mass / (AllVars.Lookback_Time[SnapList_Millennium[i]] - AllVars.Lookback_Time[SnapList_Millennium[i]+1]) / AllVars.Sec_Per_Megayear / 1e3
      #Photons_H_Millennium.append(np.log10(H_Millennium.Mvir[w_H_Millennium[i]]* 1.0e10 / AllVars.Hubble_h * Photon_Factor))
      #Photons_Tot_H_Millennium.append(np.log10(sum(10**Photons_H_Millennium[i])))

      SnapListZ.append(AllVars.SnapZ[SnapList_Millennium[i]])
      SnapListZ_Millennium.append(AllVars.SnapZ[SnapList_Millennium[i]])

      ## Model 1 Calculations ##

      w_G_Millennium.append(np.where((G_Millennium.GridHistory[:, SnapList_Millennium[i]] != -1) & (G_Millennium.GridStellarMass[:, SnapList_Millennium[i]] > 0.0) & (G_Millennium.GridSFR[:, SnapList_Millennium[i]] > 0.0) & (G_Millennium.GridCentralGalaxyMass[:,SnapList_Millennium[i]] > 0.0))[0])
      mass_G_Millennium.append(np.log10(G_Millennium.GridStellarMass[w_G_Millennium[i], SnapList_Millennium[i]] * 1.0e10 / AllVars.Hubble_h))
      mass_Central_Millennium.append(np.log10(G_Millennium.GridCentralGalaxyMass[w_G_Millennium[i], SnapList_Millennium[i]] * 1.0e10 / AllVars.Hubble_h))
      Outflow_G_Millennium.append(G_Millennium.GridOutflowRate[w_G_Millennium[i], SnapList_Millennium[i]])

      print Outflow_G_Millennium
      ## Model 2 Calculations ##

      w_G_Millennium2.append(np.where((G_Millennium2.GridHistory[:, SnapList_Millennium[i]] != -1) & (G_Millennium2.GridStellarMass[:, SnapList_Millennium[i]] > 0.0) & (G_Millennium2.GridSFR[:, SnapList_Millennium[i]] > 0.0) & (G_Millennium2.GridCentralGalaxyMass[:,SnapList_Millennium[i]] > 0.0))[0])
      mass_G_Millennium2.append(np.log10(G_Millennium2.GridStellarMass[w_G_Millennium2[i], SnapList_Millennium[i]] * 1.0e10 / AllVars.Hubble_h))
      mass_Central_Millennium2.append(np.log10(G_Millennium2.GridCentralGalaxyMass[w_G_Millennium2[i], SnapList_Millennium[i]] * 1.0e10 / AllVars.Hubble_h))
      Outflow_G_Millennium2.append(G_Millennium2.GridOutflowRate[w_G_Millennium2[i], SnapList_Millennium[i]])

      '''

  ## Model 3 Calculations ##

  w_G_model3.append(np.where((G_model3.GridHistory[:, SnapList[i]] != -1) & (G_model3.GridStellarMass[:, SnapList[i]] > 0.0) & (G_model3.GridSFR[:, SnapList[i]] > 0.0) & (G_model3.GridCentralGalaxyMass[:,SnapList[i]] > 0.0))[0])
  mass_G_model3.append(np.log10(G_model3.GridStellarMass[w_G_model3[i], SnapList[i]] * 1.0e10 / AllVars.Hubble_h))
  Photons_G_model3.append(np.log10(10**(G_model3.Photons[w_G_model3[i], SnapList[i]])*fesc_model3))
  Photons_Tot_G_model3.append(np.log10(sum(10**Photons_G_model3[i])))
  SFR_G_model3.append(np.log10(G_model3.GridSFR[w_G_model3[i], SnapList[i]]))
  sSFR_G_model3.append(SFR_G_model3[i] - mass_G_model3[i])

  mass_Central_model3.append(np.log10(G_model3.GridCentralGalaxyMass[w_G_model3[i], SnapList[i]] * 1.0e10 / AllVars.Hubble_h))
  Photon_Factor = AllVars.Solar_Mass * SourceEfficiency_model3 * AllVars.BaryonFrac / AllVars.Ionized_Mass_H / AllVars.Proton_Mass / (AllVars.Lookback_Time[SnapList[i]] - AllVars.Lookback_Time[SnapList[i]+1]) / AllVars.Sec_Per_Megayear / 1e3
  Photons_Central_model3.append(np.log10(10**(mass_Central_model3[i]) * Photon_Factor))
  Photons_Tot_Central_model3.append(np.log10(sum(10**Photons_Central_model3[i])))
      '''

def Calculate_HaloPartStellarMass(ZZ, HaloPart, StellarMass, BoundLow, BoundHigh):

	for i in xrange(0, len(ZZ)):	
		w = np.where((HaloPart[i] > BoundLow) & (HaloPart[i] < BoundHigh))[0]

		Mass = np.mean(10**(StellarMass[i][w]))
		Mass_std = np.std(10**(StellarMass[i][w]))
		count = len(StellarMass[i][w])
		print "Mass %.4e \t Mass_std = %.4e \t Count = %d" %(Mass, Mass_std, count) 

	return Mass

#print "For 512 model:"
#HaloPartStellarMass_MySim = Calculate_HaloPartStellarMass(SnapListZ, HaloPart_MySim, mass_G_MySim, HaloPart_Low, HaloPart_High)

low = 45
high = 55
#print "For 1024 model:"
#HaloPartStellarMass_MySim = Calculate_HaloPartStellarMass(SnapListZ, HaloPart_MySim, mass_G_MySim, low, high)
#HaloPartStellarMass_MySim = Calculate_HaloPartStellarMass(SnapListZ, HaloPart_MySim, mass_G_Millennium, low, high)
HaloPartStellarMass_MySim = 0

print "Mass_G_MySim", mass_G_MySim
print "Mass_G_MySim2", mass_G_MySim2
#Metallicity(Simulation, SnapListZ, mass_G_MySim, Metallicity_Tremonti_G_model1)
#Photon_Totals(Simulation, [SnapListZ_MySim, SnapListZ_MySim, SnapListZ_MySim, SnapListZ_MySim], [Photons_Tot_Central_MySim, Photons_Tot_G_MySim, Photons_Tot_Central_MySim2, Photons_Tot_G_MySim2], len(SnapList_MySim))
#StellarMassFunction(Simulation, SnapListZ, (mass_G_MySim + mass_G_Millennium + mass_G_MySim2 + mass_G_Millennium2), HaloPartStellarMass_MySim, len(SnapListZ), ["Instant", "Delayed - 1Myr"], "MiniMillenniumSN1_LocalIRA")
#HaloMassFunction(Simulation, SnapListZ, (mass_H_MySim + mass_H_MySim2 + mass_H_Millennium), len(SnapList_MySim)) 
#CentralGalaxy_Comparison(Simulation, SnapListZ_MySim, (mass_Central_MySim2 + mass_Central_MySim2), (Photons_Central_MySim2 + Photons_G_MySim2))
#CentralGalaxy_Comparison_Difference(Simulation, SnapListZ, (mass_Central_MySim + mass_Central_model1), (Photons_Central_model1 + Photons_G_model1))
#CentralGalaxy_Projection(Simulation, SnapListZ, (mass_Central_MySim + mass_Central_MySim2), (Photons_G_MySim + Photons_G_MySim2))  
#FilteringMass(Simulation, SnapListZ_MySim, [MfiltGnedin_MySim, MfiltSobacchi_MySim, 10**mass_Central_MySim], len(SnapList_MySim), [r"\mathrm{Gnedin}", r'Sobacchi Recursive $\texttt{SAGE}, f_\mathrm{esc} \: \propto \: M_H^{\beta}$', r"$M_\mathrm{Vir} \: \mathrm{Base} \: \texttt{SAGE}$"])
#PlotScripts.Plot_Scatter_SixPanel((mass_Central_model1 + mass_Central_model1 + mass_Central_model3), (Photons_Central_model1 + Photons_G_model1 + Photons_G_model3), 1, 2, ['Halos', 'Fiducial Gal', 'NoSN Gal'], [min_mass_H, max_mass_H, -2, 1], [r'Log Halo Mass [$M_{\odot}$]', r'Log (Halo/Galaxy Ionizing Photons) [s$^{-1}$]'], 2, SnapListZ, 'CentralStellar_Photons_Difference_Fiducial_NoSN', '.png')
#StarFormationRate_Redshift(Simulation, SnapListZ_MySim, SFR_G_MySim, SFR_G_MySim2)
#PhotonsStellarMass(Simulation, SnapListZ_MySim, mass_G_MySim2, Photons_HI_G_MySim2)
#LymanAlphaLF(Simulation, SnapListZ_MySim, LymanAlpha_MySim2, len(SnapList_MySim))
#PhotonsVsStellarMass(Simulation, SnapListZ_MySim, mass_G_MySim2, Photons_HI_G_MySim2)
#fesc(Simulation, SnapListZ_MySim, mass_Central_MySim2, fesc_Kimm_MySim2)
#UVLF(Simulation, SnapListZ, (MUV_MySim), (MUV_Obs_MySim), len(SnapList_MySim), r"Recursive $SAGE, f_\mathrm{esc} \: \propto \: M_H^{-\beta}$", "MH_pos_UVLF")
#UVLF(Simulation, SnapListZ, (MUV_MySim2), (MUV_Obs_MySim2), len(SnapList_MySim), r"Recursive $SAGE, f_\mathrm{esc} = 0.10$", "fesc0.10_UVLF")
#SFR_Hist(Simulation, SnapListZ, (SFR_G_MySim + SFR_G_MySim2), len(SnapList_MySim)) 
#SFRVsStellarMass(Simulation, SnapListZ, mass_G_MySim, SFR_G_MySim)
#sSFR_Hist(Simulation, SnapListZ, (sSFR_G_MySim), len(SnapList_MySim)) 
#fesc(Simulation, SnapListZ_MySim, [fesc_local_MySim, fesc_local2_MySim, fesc_local3_MySim], [r'$f_\mathrm{esc} = 0.25$', r'$f_\mathrm{esc} \: \propto \: M_H^{-\beta}$', r'$f_\mathrm{esc} \: \propto \: m_\mathrm{Ejected}$'], "fesc_HaloPartCut100")
#EjectedFracVsStellarMass(Simulation, SnapListZ_MySim, mass_Central_MySim, EjectedFraction_MySim, len(SnapList_MySim), "EjectedMass_HaloMass_HaloPartCut100")
#HaloPartCount(Simulation, SnapListZ_MySim, HaloCount_MySim, len(SnapList_MySim))
#SFR(Simulation, SnapListZ, mass_G_MySim, SFR_G_MySim, len(SnapList_MySim), "SFR")
OutflowRate(SnapList_Millennium, mass_G_Millennium, mass_G_Millennium2, Outflow_G_Millennium, Outflow_G_Millennium2, ["Instant", "Delayed - 1Myr"], "MiniMillennium_Outflow_Stellar")
