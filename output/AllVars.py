#!/usr/bin/env python

from __future__ import print_function
import numpy as np
from astropy import units as u
from astropy import cosmology
from scipy import stats
import os

def set_cosmology(Hubble_h, Omega_m):
	
    cosmo = cosmology.FlatLambdaCDM(H0 = Hubble_h*100, Om0 = Omega_m) 
    t_BigBang = cosmo.lookback_time(100000).value # Lookback time to the Big Bang in Gyr.

    return cosmo, t_BigBang

def Set_Params_MiniMill():
    
    print("Setting parameters to Mini Millennium.")
    
    global Hubble_h
    global Omega_m
    global Omega_L
    global BoxSize
    global Volume
    global SnapZ
    global BaryonFrac
    
    global SnapZ
    global Lookback_Time

    global cosmo
    global t_BigBang
    
    Hubble_h = 0.73
    Omega_m = 0.25
    Omega_L = 0.75
    BoxSize = 62.5 # Mpc/h
    Volume = BoxSize**3
    BaryonFrac = 0.17
    
    SnapZ = [ 1.27000000e+02, 7.99978940e+01, 4.99995900e+01, 3.00000630e+01,
             1.99156900e+01,  1.82437230e+01, 1.67245250e+01, 1.53430740e+01,
             1.40859140e+01,  1.29407800e+01, 1.18965700e+01, 1.09438640e+01,
             1.00734610e+01,  9.27791500e+00, 8.54991200e+00, 7.88320350e+00,
             7.27218800e+00,  6.71158650e+00, 6.19683360e+00, 5.72386400e+00,
             5.28883360e+00,  4.88844900e+00, 4.51955560e+00, 4.17946860e+00,
             3.86568280e+00,  3.57590500e+00, 3.30809780e+00, 3.06041900e+00,
             2.83118270e+00,  2.61886140e+00, 2.42204400e+00, 2.23948550e+00,
             2.07002740e+00,  1.91263270e+00, 1.76633580e+00, 1.63027070e+00,
             1.50363650e+00,  1.38571810e+00, 1.27584620e+00, 1.17341690e+00,
             1.07787450e+00,  9.88708140e-01, 9.05462400e-01, 8.27699100e-01,
             7.55035640e-01,  6.87108800e-01, 6.23590100e-01, 5.64176600e-01,
             5.08591400e-01,  4.56577240e-01, 4.07899440e-01, 3.62340270e-01,
             3.19703430e-01,  2.79801800e-01, 2.42469090e-01, 2.07548630e-01,
             1.74897610e-01,  1.44383420e-01, 1.15883370e-01, 8.92878300e-02,
             6.44933950e-02,  4.14030630e-02, 1.99325420e-02, 0.00000000e+00]

    Lookback_Time = [13.5672,    13.5551,    13.5305,    13.4761,    13.3929,    13.3680,    13.3403,
                 13.3093,    13.2748,    13.2365,    13.1940,    13.1470,    13.0951,    13.0378,
                 12.9748,    12.9055,    12.8296,    12.7465,    12.6558,    12.5569,    12.4494,
                 12.3327,    12.2064,    12.0699,    11.9227,    11.7644,    11.5945,    11.4127,
                 11.2186,    11.0119,    10.7924,    10.5598,    10.3142,    10.0557,    9.78421,
                 9.50011,    9.20376,    8.89563,    8.57635,    8.24662,    7.90726,    7.55927,
                 7.20365,    6.84147,    6.47396,    6.10229,    5.72772,    5.35152,    4.97492,
                 4.59916,    4.22544,    3.85492,    3.48873,    3.12788,    2.77339,    2.42617,
                 2.08709,    1.75692,    1.43640,    1.12622,    0.82696,    0.53917,
                 0.263375,   0.      ] # In Gyr. 

    cosmo, t_BigBang = set_cosmology(Hubble_h, Omega_m)

    print("######################")
    print("BoxSize = %.3f (Mpc/h)" %(BoxSize))
    print("Hubble_h = %.3f" %(Hubble_h))
    print("Omega_m = %.3f" %(Omega_m))
    print("Omega_L = %.3f" %(Omega_L))
    print("BaryonFrac = %.3f" %(BaryonFrac))
    print("t_BigBang = %.3f Gyr" %(t_BigBang))
    print("######################")

def Set_Params_Mysim():

    print("Setting parameters to my simulation.")
    
    global Hubble_h
    global Omega_m
    global Omega_L
    global Omega_b
    global BoxSize
    global Volume
    global SnapZ
    global BaryonFrac
    global Y
    
    global SnapZ
    global Lookback_Time

    global cosmo
    global t_BigBang
    
    Hubble_h = 0.678
    Omega_m = 0.308
    Omega_L = 0.692
    Omega_b = 0.0484
    BoxSize = 100 # Mpc/h
    Volume = BoxSize**3
    BaryonFrac = 0.17
    Y = 0.24

 
    SnapZ = [49.000000,     35.000001,      32.025452,      29.597752,      27.571988,      25.851283,      24.368212,
             23.074276,     21.933618,      20.919111,      20.009818,      19.189299,      18.444455,      17.764701,
             17.141388,     16.567374,      16.036700,      15.544359,      15.086102,      14.658304,      14.257849,
             13.882045,     13.528548,      13.195313,      12.880543,      12.582654,      12.300241,      12.032056,
             11.776987,     11.534035,      11.302305,      11.080989,      10.869359,      10.666753,      10.472570,
             10.286264,     10.107335,      9.935327,       9.769820,       9.610431,       9.456804,       9.308614,
             9.165559,      9.027360,       8.893758,       8.764514,       8.639403,       8.518218,       8.400767,
             8.286867,      8.176350,       8.069058,       7.964843,       7.863565,       7.765095,       7.669310,
             7.576095,      7.485340,       7.396944,       7.310809,       7.226845,       7.144965,       7.065088,
             6.987135,      6.911035,       6.836717,       6.764116,       6.693168,       6.623815,       6.555999,
             6.489667,      6.424768,       6.361251,       6.299072,       6.238185,       6.178547,       6.120119,
             6.062861,      6.006736,       5.951709,       5.897746,       5.844813,       5.792881,       5.741918,
             5.691896,      5.642788,       5.594568,       5.547208,       5.500687,       5.454979,       5.410062,
             5.365914,      5.322515,       5.279845,       5.237883,       5.196612,       5.156013,       5.116069,
             5.076763,      5.038078,       5.000000]

    Lookback_Time = [13.7561, 13.7249, 13.7138, 13.7028, 13.6917, 13.6806, 13.6695,
                 13.6584, 13.6474, 13.6363, 13.6252, 13.6141, 13.6030, 13.5920,
                 13.5809, 13.5698, 13.5587, 13.5477, 13.5366, 13.5255, 13.5144,
                 13.5033, 13.4923, 13.4812, 13.4701, 13.4590, 13.4479, 13.4369,
                 13.4258, 13.4147, 13.4036, 13.3926, 13.3815, 13.3704, 13.3593,
                 13.3482, 13.3372, 13.3261, 13.3150, 13.3039, 13.2928, 13.2818,
                 13.2707, 13.2596, 13.2485, 13.2375, 13.2264, 13.2153, 13.2042,
                 13.1931, 13.1821, 13.1710, 13.1599, 13.1488, 13.1377, 13.1267,
                 13.1156, 13.1045, 13.0934, 13.0824, 13.0713, 13.0602, 13.0491,
                 13.0380, 13.0270, 13.0159, 13.0048, 12.9937, 12.9826, 12.9716,
                 12.9605, 12.9494, 12.9383, 12.9273, 12.9162, 12.9051, 12.8940,
                 12.8829, 12.8719, 12.8608, 12.8497, 12.8386, 12.8275, 12.8165,
                 12.8054, 12.7943, 12.7832, 12.7722, 12.7611, 12.7500, 12.7389,
                 12.7278, 12.7168, 12.7057, 12.6946, 12.6835, 12.6724, 12.6614,
                 12.6503, 12.6392, 12.6281] # In Gyr.


    cosmo, t_BigBang = set_cosmology(Hubble_h, Omega_m)

    print("######################")
    print("BoxSize = %.3f (Mpc/h)" %(BoxSize))
    print("Hubble_h = %.3f" %(Hubble_h))
    print("Omega_m = %.3f" %(Omega_m))
    print("Omega_L = %.3f" %(Omega_L))
    print("BaryonFrac = %.3f" %(BaryonFrac))
    print("t_BigBang = %.3f Gyr" %(t_BigBang))
    print("######################")

    return cosmo

def Set_Constants():

    print("Setting constants (in cgs units).")

    global Proton_Mass
    global Solar_Mass
    global Sec_Per_Year
    global Sec_Per_Megayear
    global Ionized_Mass_H
    global LymanAlpha_Energy
    global eV_to_erg
    global M_Bol_Sun
    global L_Sun
    global W_to_ergs
    global A_to_m
    global c_in_ms
    global pc_to_m
    global m_to_cm
    global Sigmat
    global H0
    global G
    global mpc_to_m 
    global cubic_cm_per_cubic_mpc
    global km_to_m
    global density_baryons  
 
    Proton_Mass = 1.6726219e-24 # Grams.
    Solar_Mass = 1.98855e33 # Grams.
    Ssec_Per_Year = 3.155e7
    Sec_Per_Megayear = 3.155e13
    Ionized_Mass_H = 0.53 # Molecular mass of ionized H.
    LymanAlpha_energy = 10.2 # Energy of Lyman Alpha photon in eV.
    eV_to_erg = 1.6202e-12 # Conversion from eV to erg. 
    M_bol_Sun = 4.74 # Bolometric magnitude of the sun.
    L_sun = 3.828e26 # Luminosity of the Sun in W.
    W_to_ergs = 1.0e7 # Watts to erg s^-1.
    A_to_m = 1.0e-10 # Angstroms to meters.
    c_in_ms = 3.0e8 # Speed of light in m s^-1 
    pc_to_m = 3.086e16 # Parsec to metres. 
    m_to_cm = 1.0e2 # Meters to centimetres.
    Sigmat = 6.652e-29 # Thomson cross-section for an electron in m^2.
    H0 = 3.24078e-18 # In h/sec.
    G = 6.674e-11 # Gravitational constant in m^3 kg^-1 s^-2.
    mpc_to_m = 3.0857e22 # How many metres in a Mpc 
    cubic_cm_per_cubic_mpc = 2.93799e73 # cm^-3 to Mpc_3
    km_to_m = 1.0e3 # Kilometres to metres.
    density_baryons = 2.5e-7 # Number density of baryons in cm^-3 (flat Universe).


def spectralflux_wavelength_to_frequency(Flux, Wavelength):

    # For a given spectral flux at a specific wavelength (f_lamba), this function converts it to the spectral flux density at a specific frequency (f_nu).

    ### INPUT ###
    # Flux: The spectral flux density in units of erg s^-1 A^-1 cm^-2
    # Wavelength: The wavelength we are determining the spectral flux density at in units of Angstroms (A).

    ### OUTPUT ###
    # f_nu: The spectral flux density at a specific frequency in units of Janksy (W Hz^-1 m^-2).

    f_nu = 3.34e4 * pow(Wavelength, 2) * Flux
  
    return f_nu

def Luminosity_to_Flux(Luminosity, Distance):

    # Converts a luminosity to an observed flux at some distance.
    ## NOTE THAT INPUT AND OUTPUT ARE IN LOG UNITS.

    ### INPUT ###
    # Luminosity: The intrinisic luminosity being converted in units of erg s^-1 A^-1. <<<NOTE MUST BE IN LOG UNITS>>>.
    # Distance: The distance at which the flux is being observed in units of parsec.

    ### OUTPUT ###
    # Flux: The observed flux in units of erg s^-1 A^-1 cm^-2. 

    F = Luminosity - np.log10(4*np.pi*pow(Distance * pc_to_m * m_to_cm, 2.0))

    return F

def Luminosity_to_ABMag(Luminosity, Wavelength):
     
    # Converts an intrinsic luminosity into absolute AB magnitude at a specified wavelength.
    ## NOTE THAT INPUT IS IN LOG.

    ### INPUT ###
    # Luminosity: The intrinsic luminosity of the star in units of erg s^-1 A^-1.  <<NOTE MUST BE IN LOG UNITS>>>.
    # Wavelength: The wavelength we want the magnitude at.

    ### OUTPUT ###
    # M: The absolute magnitude in the AB system.

    Flux = Luminosity_to_Flux(Luminosity, 10.0) # Calculate the flux from a distance of 10 parsec, units of erg s^-1 A^-1 cm^-2.  Log Units. 
    f_nu = spectralflux_wavelength_to_frequency(10**Flux, 1600) # Spectral flux density in Janksy.
    M = -2.5 * np.log10(f_nu) + 8.90 # AB Magnitude from http://www.astro.ljmu.ac.uk/~ikb/convert-units/node2.html

    #print "Flux from AllVars.py", Flux
    #print "M from AllVars.py", M
    return M


def Set_Params_Tiamat():

    print("Setting parameters to Tiamat")
    
    global Hubble_h
    global Omega_m
    global Omega_L
    global Omega_b
    global BoxSize
    global Volume
    global SnapZ
    global BaryonFrac
    global Y

    global SnapZ
    global Lookback_Time

    global cosmo
    global t_BigBang
    
    Hubble_h = 0.678
    Omega_m = 0.308
    Omega_L = 0.692
    Omega_b = 0.0484
    BoxSize = 100*Hubble_h # Mpc/h
    Volume = BoxSize**3
    BaryonFrac = 0.17
    Y = 0.24   
 
    SnapZ = [49.000000,     35.000001,      32.025452,      29.597752,      27.571988,      25.851283,      24.368212,
             23.074276,     21.933618,      20.919111,      20.009818,      19.189299,      18.444455,      17.764701,
             17.141388,     16.567374,      16.036700,      15.544359,      15.086102,      14.658304,      14.257849,
             13.882045,     13.528548,      13.195313,      12.880543,      12.582654,      12.300241,      12.032056,
             11.776987,     11.534035,      11.302305,      11.080989,      10.869359,      10.666753,      10.472570,
             10.286264,     10.107335,      9.935327,       9.769820,       9.610431,       9.456804,       9.308614,
             9.165559,      9.027360,       8.893758,       8.764514,       8.639403,       8.518218,       8.400767,
             8.286867,      8.176350,       8.069058,       7.964843,       7.863565,       7.765095,       7.669310,
             7.576095,      7.485340,       7.396944,       7.310809,       7.226845,       7.144965,       7.065088,
             6.987135,      6.911035,       6.836717,       6.764116,       6.693168,       6.623815,       6.555999,
             6.489667,      6.424768,       6.361251,       6.299072,       6.238185,       6.178547,       6.120119,
             6.062861,      6.006736,       5.951709,       5.897746,       5.844813,       5.792881,       5.741918,
             5.691896,      5.642788,       5.594568,       5.547208,       5.500687,       5.454979,       5.410062,
             5.365914,      5.322515,       5.279845,       5.237883,       5.196612,       5.156013,       5.116069,
             5.076763,      5.038078,       5.000000,       5.000000]

    Lookback_Time = [13.7561, 13.7249, 13.7138, 13.7028, 13.6917, 13.6806, 13.6695,
                 13.6584, 13.6474, 13.6363, 13.6252, 13.6141, 13.6030, 13.5920,
                 13.5809, 13.5698, 13.5587, 13.5477, 13.5366, 13.5255, 13.5144,
                 13.5033, 13.4923, 13.4812, 13.4701, 13.4590, 13.4479, 13.4369,
                 13.4258, 13.4147, 13.4036, 13.3926, 13.3815, 13.3704, 13.3593,
                 13.3482, 13.3372, 13.3261, 13.3150, 13.3039, 13.2928, 13.2818,
                 13.2707, 13.2596, 13.2485, 13.2375, 13.2264, 13.2153, 13.2042,
                 13.1931, 13.1821, 13.1710, 13.1599, 13.1488, 13.1377, 13.1267,
                 13.1156, 13.1045, 13.0934, 13.0824, 13.0713, 13.0602, 13.0491,
                 13.0380, 13.0270, 13.0159, 13.0048, 12.9937, 12.9826, 12.9716,
                 12.9605, 12.9494, 12.9383, 12.9273, 12.9162, 12.9051, 12.8940,
                 12.8829, 12.8719, 12.8608, 12.8497, 12.8386, 12.8275, 12.8165,
                 12.8054, 12.7943, 12.7832, 12.7722, 12.7611, 12.7500, 12.7389,
                 12.7278, 12.7168, 12.7057, 12.6946, 12.6835, 12.6724, 12.6614,
                 12.6503, 12.6392, 12.6281, 12.6281] # In Gyr.


    cosmo, t_BigBang = set_cosmology(Hubble_h, Omega_m)

    print("######################")
    print("BoxSize = %.3f (Mpc/h)" %(BoxSize))
    print("Hubble_h = %.3f" %(Hubble_h))
    print("Omega_m = %.3f" %(Omega_m))
    print("Omega_L = %.3f" %(Omega_L))
    print("BaronFrac = %.3f" %(BaryonFrac))
    print("t_BigBang = %.3f Gyr" %(t_BigBang))
    print("######################")

    return cosmo


def Set_Params_Tiamat_extended():

    print("Setting parameters to the extended Tiamat simulation (that ran down to z = 1.8ish).") 
    
    global Hubble_h
    global Omega_m
    global Omega_L
    global Omega_b
    global BoxSize
    global Volume
    global SnapZ
    global BaryonFrac
    global Y
    global PartMass

    global SnapZ
    global Lookback_Time

    global cosmo
    global t_BigBang
    
    Hubble_h = 0.678
    Omega_m = 0.308
    Omega_L = 0.692
    Omega_b = 0.0484
    BoxSize = 100*Hubble_h # Mpc/h
    Volume = BoxSize**3
    BaryonFrac = 0.17
    Y = 0.24   
 
    PartMass = 2.6438e6 # Msun/h 

    SnapZ = [49.000000,     35.000001,      32.025452,      29.597752,      27.571988,      25.851283,      24.368212,
             23.074276,     21.933618,      20.919111,      20.009818,      19.189299,      18.444455,      17.764701,
             17.141388,     16.567374,      16.036700,      15.544359,      15.086102,      14.658304,      14.257849,
             13.882045,     13.528548,      13.195313,      12.880543,      12.582654,      12.300241,      12.032056,
             11.776987,     11.534035,      11.302305,      11.080989,      10.869359,      10.666753,      10.472570,
             10.286264,     10.107335,      9.935327,       9.769820,       9.610431,       9.456804,       9.308614,
             9.165559,      9.027360,       8.893758,       8.764514,       8.639403,       8.518218,       8.400767,
             8.286867,      8.176350,       8.069058,       7.964843,       7.863565,       7.765095,       7.669310,
             7.576095,      7.485340,       7.396944,       7.310809,       7.226845,       7.144965,       7.065088,
             6.987135,      6.911035,       6.836717,       6.764116,       6.693168,       6.623815,       6.555999,
             6.489667,      6.424768,       6.361251,       6.299072,       6.238185,       6.178547,       6.120119,
             6.062861,      6.006736,       5.951709,       5.897746,       5.844813,       5.792881,       5.741918,
             5.691896,      5.642788,       5.594568,       5.547208,       5.500687,       5.454979,       5.410062,
             5.365914,      5.322515,       5.279845,       5.237883,       5.196612,       5.156013,       5.116069,
             5.076763,      5.038078,       5.000000,       4.928760,       4.858360,	    4.788800,
             4.720060,      4.652140,       4.585030,       4.518710,       4.453180,       4.388430,       4.324450,
             4.261230,      4.198750,       4.137020,       4.076030,       4.015750,       3.956190,       3.897350,
             3.839190,      3.781730,       3.724960,       3.668850,       3.613410,       3.558630,       3.504500,
	     3.451020,      3.398170,       3.345940,       3.294340,       3.243350,       3.192960,       3.143170,
	     3.093980,      3.045370,       2.997330,       2.949870,       2.902970,       2.856620,       2.810830,
	     2.765580,      2.720870,       2.676690,       2.633030,       2.589890,       2.547260,       2.505140,
	     2.463520,      2.422400,       2.381760,       2.341610,       2.301930,       2.262720,       2.223980,
	     2.185700,      2.147870,       2.110490,       2.073560,       2.037060,       2.001000,       1.965370,
	     1.930160,      1.895360,       1.860980,       1.82701]

    cosmo, t_BigBang = set_cosmology(Hubble_h, Omega_m)

    Lookback_Time = cosmo.lookback_time(SnapZ).value # In Gyr 
                             
    print("######################")
    print("BoxSize = %.3f (Mpc/h)" %(BoxSize))
    print("Hubble_h = %.3f" %(Hubble_h))
    print("Omega_m = %.3f" %(Omega_m))
    print("Omega_L = %.3f" %(Omega_L))
    print("BaryonFrac = %.3f" %(BaryonFrac))
    print("t_BigBang = %.3f Gyr" %(t_BigBang))
    print("PartMass = %.3f Msun/h" %(PartMass))
    print("######################")

    return cosmo


def Set_Params_Simfast21():

    print("Setting parameters to what I use for Simfast()") 
    
    global Hubble_h
    global Omega_m
    global Omega_L
    global Omega_b
    global BoxSize
    global Volume
    global SnapZ
    global BaryonFrac
    global Y

    global SnapZ
    global Lookback_Time

    global cosmo
    global t_BigBang
    
    Hubble_h = 0.678
    Omega_m = 0.308
    Omega_L = 0.692
    Omega_b = 0.0484
    BoxSize = 100*Hubble_h # Mpc/h
    Volume = BoxSize**3
    BaryonFrac = 0.17
    Y = 0.24   
 
    SnapZ = np.arange(15.0, 5.75, -0.25) 
    cosmo, t_BigBang = set_cosmology(Hubble_h, Omega_m)

    Lookback_Time = cosmo.lookback_time(SnapZ).value # In Gyr 


    print("######################")
    print("BoxSize = %.3f (Mpc/h)" %(BoxSize))
    print("Hubble_h = %.3f" %(Hubble_h))
    print("Omega_m = %.3f" %(Omega_m))
    print("Omega_L = %.3f" %(Omega_L))
    print("BaryonFrac = %.3f" %(BaryonFrac))
    print("t_BigBang = %.3f Gyr" %(t_BigBang))
    print("######################")
                                 
    return cosmo

def Set_Params_Britton():

    print("Setting parameters to Britton's Simulation") 
    
    global Hubble_h
    global Omega_m
    global Omega_L
    global Omega_b
    global BoxSize
    global Volume
    global SnapZ
    global BaryonFrac
    global Y
    global PartMass

    global SnapZ
    global Lookback_Time

    global cosmo
    global t_BigBang
    
    Hubble_h = 0.695
    Omega_m = 0.285
    Omega_L = 0.715
    Omega_b = 0.0
    BoxSize = 50.0 # Mpc/h
    Volume = BoxSize**3
    BaryonFrac = 0.17
    Y = 0.24   
    PartMass = 1.151e6 # Msun/h 

    a = np.loadtxt("/lustre/projects/p134_swin/jseiler/subfind_britton/trees/britton_shifted/a_list.txt")
    SnapZ = 1.0/a - 1
       
    cosmo, t_BigBang = set_cosmology(Hubble_h, Omega_m)

    Lookback_Time = cosmo.lookback_time(SnapZ).value # In Gyr 


    print("######################")
    print("BoxSize = %.3f (Mpc/h)" %(BoxSize))
    print("Hubble_h = %.3f" %(Hubble_h))
    print("Omega_m = %.3f" %(Omega_m))
    print("Omega_L = %.3f" %(Omega_L))
    print("BaryonFrac = %.3f" %(BaryonFrac))
    print("t_BigBang = %.3f Gyr" %(t_BigBang))
    print("PartMass = %.3f Msun/h" %(PartMass))
    print("######################")
                                 
    return cosmo

def Set_Params_Kali():

    print("Setting parameters to Kali") 
    
    global Hubble_h
    global Omega_m
    global Omega_L
    global Omega_b
    global BoxSize
    global Volume
    global SnapZ
    global BaryonFrac
    global Y
    global PartMass

    global SnapZ
    global Lookback_Time

    global cosmo
    global t_BigBang
    
    Hubble_h = 0.681
    Omega_m = 0.302
    Omega_L = 0.698
    Omega_b = 0.0452
    BoxSize = 108.96 # Mpc/h
    Volume = BoxSize**3
    BaryonFrac = 0.17
    Y = 0.24   
    PartMass = 7.8436e6# Msun/h 

    a = np.loadtxt("/lustre/projects/p134_swin/jseiler/kali/a_list.txt")
    a = a[:99]
    SnapZ = 1.0/a - 1
       
    cosmo, t_BigBang = set_cosmology(Hubble_h, Omega_m)

    Lookback_Time = cosmo.lookback_time(SnapZ).value # In Gyr 

    print("######################")
    print("BoxSize = %.3f (Mpc/h)" %(BoxSize))
    print("Hubble_h = %.3f" %(Hubble_h))
    print("Omega_m = %.3f" %(Omega_m))
    print("Omega_L = %.3f" %(Omega_L))
    print("BaryonFrac = %.3f" %(BaryonFrac))
    print("t_BigBang = %.3f Gyr" %(t_BigBang))
    print("PartMass = %.3f Msun/h" %(PartMass))
    print("######################")
                                 
    return cosmo


def depth(l):
    '''
    Determines the nested level of a list.

    Parameters
    ----------
    l : `np.darray'
	List that we are determining the nested level of.

    Returns
    ------
    int
	Nested level of the list.
    '''    


    if isinstance(l, list):
        return 1 + max(depth(item) for item in l)
    else:
        return 0


## What snapshot's in Britton's simulation most closely match the snapshot's in Tiamat.
## Input should be two arrays
## Output should be an array with length equal to Tiamat Snapshot where the array should be [BrittonSnapNum Nearest TiamatSnapshot0, BrittonSnapNum Nearest TiamatSnapshot1, ..., BrittonSnapNum Nearest TiamatSnapshotN]

def find_nearest_redshifts(SnapZ_Sim1, SnapZ_Sim2):
    '''
    Given the redshift for each snapshot in two simulations, this function calculates the snapshot in simulation 2 that has the nearest redshift to that in simulation 1.

    Parameters
    ----------
    SnapZ_Sim1, SnapZ_Sim2 : 'array-like' of floats.
        The redshift at each snapshot for the two simulations.

    Returns
    -------

    map_sim1_sim2 : 'array-like' of floats with same length as SnapZ_sim1.
        The snapshots in simulation 2 that have the closest redshift to those in simulation 1. 

    '''


    map_sim1_sim2 = np.empty((len(SnapZ_Sim1)), dtype = np.int32) 

    # Let's do it the stupid way first.
    # For each redshift in Sim1 we will find the redshift difference for each Sim2 redshift. 
    # Then will find the lowest difference and record it.

    for snapshot_idx_sim1 in range(len(SnapZ_Sim1)):
        difference = []
        for snapshot_idx_sim2 in range(len(SnapZ_Sim2)):
            difference.append(np.abs(SnapZ_Sim1[snapshot_idx_sim1] - SnapZ_Sim2[snapshot_idx_sim2]))
        map_sim1_sim2[snapshot_idx_sim1] = np.argmin(difference)
 
    return map_sim1_sim2

def Calculate_Histogram(data, bin_width, weights, min_hist=None, max_hist=None):
    '''
    Calculates a 1D histogram for the input data.  Allows the calculation of the count or probability within each bin.

    Parameters
    ---------
    data : `array-like'
        Data used in the histogram.
    bin_width : `array-like'
        Width of the bins.
    weights : either 0 or 1
        Selects the binning mode.
        0 : Histogram will be a frequency (count) histogram.
        1 : Histogram will be a probability histogram.  
    min_hist, max_hist : float (optional)
        Defines the bounds that we will be binning over.
        If no values defined, range will be the minimum/maximum data point +/- 10 times the bin_width.

    Returns
    -------
    counts : `array-like'
    Array that contains the count of probability in each bin.
    bin_edges : `array-like'
    Array containing the location of the bin edges.
    bin_middle : `array-like'
    Array containing the location of the bin middles.

    Units
    -----
    All units are kept the same as the inputs.
    '''

    if (min_hist == None): 
        range_low = np.floor(min(data)) - 10*bin_width
        range_high = np.floor(max(data)) + 10*bin_width
    else:
        range_low = min_hist 
        range_high = max_hist 

    if not np.isfinite(range_low):
        raise ValueError("The low range should be finite (it's not).")

    if not np.isfinite(range_high):
        raise ValueError("The high range should be finite (it's not).")

    if range_high <= range_low:
        print("Upper bin range = %.4f, lower bing range = %.4f" %(range_high, range_low)) 
        raise ValueError("The upper bin range should be less than the lower bin range")
 
    NB = int((range_high - range_low) / bin_width) 

    if NB < 1:
        print("Number of bins = %d" %(NB))
        raise ValueError("The number of bins should be greater than one.")
    
    if (weights == 0): # Running in frequency mode.
        (counts, bin_edges) = np.histogram(data, range=(range_low, range_high), bins=NB)
    else: # Running in probability mode.
        weights = np.ones_like(data)/len(data)
        (counts, bin_edges) = np.histogram(data, range=(range_low, range_high), bins=NB, weights = weights)

    bin_middle = bin_edges[:-1] + 0.5 * bin_width

    return (counts, bin_edges, bin_middle)

##

def Calculate_2D_Mean(data_x, data_y, bin_width, min_hist_x = None, max_hist_x = None):     

    '''
    Calculates the mean of the y-data that lies within binned x-data.  
 
    Note: Scipy.binned_statistic returns NaN when there are no data points within the bin.  
    I have adjusted this so it returns 0.0 instead so I can collate the data properly.    

    Parameters
    ----------
    data_x : array-like 
        Data that will be binned and provide the bins for the y-mean.
    data_y : array-like 
        Data that will be averaged in each of the bins defined by the x-data.
    bin_width : float
        Width of each x-bin.
    min_hist_x, max_hist_x: float (optional)
        Defines the x-bounds that we will be binning over.
        If no values defined, range will be the minimum/maximum data point +/- 10 times the bin_width.

    Returns
    -------
    mean_data_y, std_data_y, sum_data_y : array-like 
        Arrays that contain the mean, standard deviation and sum for the y-data as binned by the x-data.
    N_data_y : array-like 
        Array that contains the number of data points in each x-bin.    
    bins_mid : array-like 
        The mid-point coordinate for the x-bins. 

    Units
    -----
    All units are kept the same as the inputs.
    '''

    if not np.isfinite(min_hist_x):
        raise ValueError("xmin should be finite")

    if not np.isfinite(max_hist_x):
        raise ValueError("xmax should be finite")

    if (min_hist_x == None): 
        range_low = np.floor(min(data_x)) - 10*bin_width
        range_high = np.floor(max(data_x)) + 10*bin_width
    else:
        range_low = min_hist_x 
        range_high = max_hist_x 

    if range_high <= range_low: 
        raise ValueError("The upper bin range should be less than the lower bin range")
 
    NB = round((range_high - range_low) / bin_width) 

    bins = np.arange(range_low, range_high + bin_width, bin_width)
    bins_mid = bins + bin_width/2.0
    bins_mid = bins_mid[:-1] # len(bins_mid) should be 1 less than len(bins) as the last bin doesn't have a midpoint.   

    mean_data_y, bin_edges, bin_number = stats.binned_statistic(data_x, data_y, statistic='mean', bins = bins)
    sum_data_y, bin_edges, bin_number = stats.binned_statistic(data_x, data_y, statistic='sum', bins = bins)
    N_data_y, bin_edges, bin_number = stats.binned_statistic(data_x, data_y, statistic='count', bins = bins)
    std_data_y, bin_edges, bin_number = stats.binned_statistic(data_x, data_y, statistic=np.std, bins = bins)
    
    mean_data_y[N_data_y == 0] = 0.0
    std_data_y[N_data_y == 0] = 0.0   
 
    return mean_data_y, std_data_y, N_data_y, sum_data_y, bins_mid
 
##

def ensure_dir(file_path):
    print("Checking to see if directory {0} exists".format(file_path))
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        print("It does not exist, creating.")
        os.makedirs(directory)
    else:
        print("It does exist")

## Function to find the value (and index of said value) that is nearest to a given value within an array.
## Taken from https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx], idx


## Following two functions taken from Anne's grid-model Analysis code.

def modes_to_pspec(modes, boxsize):
    """
    From a given set of fourier modes, compute the (binned) power spectrum
    with errors.

    modes   - (n,n,n) array (from ifftn of overdensity values)
    boxsize - size of box (e.g. in Mpc/h)

    returns

    kmid    - k values of power spec
    pspec   - estimate of power spec at k
    perr    - error on estimate
    """

    ngrid = modes.shape[0]
    assert(modes.shape==(ngrid, ngrid,ngrid))
    kmin, kmax, kbins, kvol = powerspec_bins(ngrid, boxsize)

    wts = np.square(modes.ravel().real) + np.square(modes.ravel().imag)

    v1 = np.bincount(kbins, weights=wts)
    powerspec = v1 * (1.0 / kvol)

    # work out error on power spectrum
    v2 = np.bincount(kbins, weights=np.square(wts))
    v0 = np.bincount(kbins)
    p_err = np.sqrt((v2*v0 - v1*v1)/(v0-1)) / kvol

    kmid_bins = 0.5 * (kmin+kmax)

    return kmid_bins, powerspec, p_err

def powerspec_bins(ngrid, boxsize):
    """
    find power spectrum bins to for a cubic grid of size ngrid^3 of fourier modes.
    Assumes the FFT convention of 0, ..., n/2, -n/2+1, ..., -1
    ngrid   - num cells on side of cube
    boxsize - size of the box in real space

    returns kmin, kmax, kbins, kvol
    kmin  - the lower bound of the bin
    kmax  - the upper bound of the bin
    kbins - index (0, ..., m) of the bin of each cell
    kvol  - the volume in k space (number of cells of that type divided by k vol of each cell)
    """

    mid = int(ngrid/2)
    # find the magnitude of the indices (i.e. ix**2+iy**2+iz**2 in the FFT convention)
    n1 = np.arange(ngrid)
    n1[1+mid:] -= ngrid
    n2 = np.square(n1)
    nmag = np.sqrt(np.add.outer(np.add.outer(n2, n2), n2)).ravel()

    nbins = (-1,) + tuple(np.arange(mid-1)+1.5) + (ngrid*2,)
    #print 'nbins', nbins
    kbins = np.digitize(nmag, nbins) - 1
    assert(kbins.min()==0)
    assert(kbins.max()==len(nbins)-2)

    # multiplier to go to k-space
    kmult = 2.0 * np.pi / boxsize

    kmin = (np.array(nbins) * kmult)[:-1]
    kmin[0] = 0

    kmax = (np.array(nbins) * kmult)[1:]
    kmax[-1] = mid * kmult * np.sqrt(3.0)

    kvol = np.bincount(kbins) * (kmult * kmult * kmult)
    return kmin, kmax, kbins, kvol

