#!/usr/bin/env python
from __future__ import division
import matplotlib
matplotlib.use('Agg')

import numpy as np
import random
import timeit
import pylab as plt
from tqdm import tqdm 

import sys
sys.path.append('/home/jseiler/SAGE-stuff/output/')

import AllVars
import PlotScripts
import ReadScripts

PlotScripts.Set_Params_Plot()
N_loops = 10
verbose = 0

def match(N_key, N_match, method):
    
    times = []
    match_fraction = []
    match_possible_fraction = []

    start_time = timeit.default_timer()
    for i in range(0, N_loops + 1):
        
        if N_key > N_match:
            largest_value = N_key
        else:
            largest_value = N_match
        key_list = random.sample(range(largest_value*2), N_key)
        if method == 0:
            key_list = sorted(key_list) 
        to_match = random.sample(range(largest_value*2), N_match)

        if method == 0:
            matched_numbers = []
            for match_idx in range(0, len(to_match)):
                is_found = 0
                search_idx = int(np.ceil(len(key_list) / 2))
                count = 0
                while is_found == 0:
                #    print("Index is {0} and the value here is {1}".format(search_idx, key_list[search_idx]))
                    count += 1
        #            print("Match_idx = {0} search_idx = {1}".format(match_idx, search_idx))
                    if (to_match[match_idx] == key_list[search_idx]):
                        is_found = 1
                    elif (len(key_list) / pow(2,count) < 1.0):
                        break
                    elif (to_match[match_idx] > key_list[search_idx]):
                        search_idx = int(search_idx + np.ceil(len(key_list) / pow(2,count + 1)))
                    else:
                        search_idx = int(search_idx - np.ceil(len(key_list) / pow(2,count + 1)))

                    if search_idx >= len(key_list):
                        search_idx = len(key_list) - 1        

                if is_found == 1:
                    matched_numbers.append(to_match[match_idx])

            if verbose == 1:
                print("We matched these numbers: {0}".format(matched_numbers))
                print("Which means we matched {1} numbers, or {0:.4f} fraction of numbers".format(len(matched_numbers) / len(to_match), len(matched_numbers)))
                print("The maximum amount of matches we could do is {0} meaning that we got {1:4f} fraction of the total amount matched.".format( len(key_list), len(matched_numbers) / len(key_list)))

        else:
  
            matched_numbers = np.intersect1d(to_match, key_list)

        match_fraction.append(len(matched_numbers) / len(to_match))
        match_possible_fraction.append(len(matched_numbers) / len(key_list))


    elapsed = timeit.default_timer() - start_time

    time = elapsed / (N_loops)
    match_fraction = np.mean(match_fraction)
    match_possible_fraction = np.mean(match_possible_fraction)

    if method == 0:
        print("Mine took {0} ms".format(time*1.0e3)) 
    else:
        print("Numpy took {0} ms".format(time*1.0e3))
 
    print("Mine took {0} ms".format(time*1.0e3)) 
    return time, match_fraction, match_possible_fraction

ax1 = plt.subplot(121)
ax2 = plt.subplot(122)

#for N_match in tqdm(np.logspace(1, 6, num = 10)):
#for N_match in tqdm(range(100, 1000, 200)):
    #N_match = int(N_match)
N_key = 3326 # One snapshot has 3326 FoF Particles.
N_match =  33499774 # One snapshot has 33,499,774 particles. 

my_time, my_match_fraction, my_match_possible_fraction = match(N_key, N_match, 0)
numpy_time, numpy_match_fraction, numpy_match_possible_fraction = match(N_key, N_match, 1)

ax1.scatter(N_match, my_time * 1.0e3, color = 'r')
ax2.scatter(N_match, my_match_fraction, color = 'r')

ax1.scatter(N_match, numpy_time * 1.0e3, color = 'b')
ax2.scatter(N_match, numpy_match_fraction, color = 'b')

ax1.scatter(np.nan, np.nan, label = 'Mine', color = 'r')
ax1.scatter(np.nan, np.nan, label = 'Numpy', color = 'b')

ax1.set_xlabel("Number to Match")
ax1.set_ylabel("Time (ms)")
ax1.set_xscale('log', nonposy='clip')

ax2.set_xlabel("Number to Match")
ax2.set_ylabel("Fraction of Numbers Matched")
ax2.set_xscale('log', nonposy='clip')

plt.tight_layout()
leg = ax1.legend(markerscale = 5.0, loc='upper left', numpoints=1, labelspacing=0.1)
leg.draw_frame(False)  # Don't want a box frame
for t in leg.get_texts():  # Reduce the size of the text
    t.set_fontsize(PlotScripts.global_legendsize)


outputFile = './MatchingComp.png'
plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
print('Saved file to {0}'.format(outputFile))
plt.close()

