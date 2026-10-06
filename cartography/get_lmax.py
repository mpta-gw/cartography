# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 15:27:19 2025

@author: kgrun
"""

import argparse
import glob

import numpy as np
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
import astropy.units as u



plt.rcParams['font.family'] = "serif"
plt.rcParams['font.sans-serif'] = "Times"

plt.rcParams['text.usetex']= False
plt.rcParams['xtick.labelsize'] = 10.0
plt.rcParams['ytick.labelsize'] = 10.0
plt.rcParams['axes.labelsize'] = 14.0







def distance_pulsars(psr_SkyCoord):
    npsr = len(psr_SkyCoord)
    separations = []
    for i in range(npsr-1):
        temp = []
        for j in range(i+1, npsr):
            temp.append(psr_SkyCoord[i].separation(psr_SkyCoord[j]).radian)
        
        separations.append(min(temp))
    
    return np.array(separations), np.array(separations)*180/np.pi
        
        
    


def dnode(l,m):
    f1 = (l-m+2)/(l-m)
    cosd = np.sin(f1*np.pi/2)
    return np.arccos(cosd)



def set_up_global_options():
    parser = argparse.ArgumentParser(description='Calculate lmax from par/tim files or pulsar pickle')
    
    parser.add_argument('--parpath', type=str, default=None,
                        help='Path to the directory with par and tim files')
    
    parser.add_argument('--outdir', type=str, default=None,
                        help='Directory for the output')
    
    return parser.parse_args()



args = set_up_global_options()


#psrlist = glob.glob(args.parpath + '/*.par')[0]
psrlist = glob.glob("../../datasets/pta_2/pars/*.par")

print('LOAD DATA AND CALCULATE MINIMAL DISTANCES', '\n')
psr_SkyCoord = []
for parname in psrlist:
    with open(parname) as pf:
        ra = None
        dec = None
        for line in pf:
            if "RA" in line:
                ra = line.split()[1]
            elif "DEC" in line:
                dec = line.split()[1]
                
            if ra is not None and dec is not None:
                break
        
        psr_SkyCoord.append(SkyCoord(ra=ra, dec=dec, unit = (u.hourangle, u.deg), frame='icrs'))
    

min_distances = distance_pulsars(psr_SkyCoord)[0]
mean_md = np.mean(min_distances)
median_md = np.median(min_distances)

lpta_mean = int(2*np.pi/mean_md)
lpta_median = int(2*np.pi/median_md)



### PLOT SPH NODE DISTANCES ###################################################
###############################################################################
print('CREATE PLOT')

ll = np.min([lpta_mean-20, lpta_median-20])
ul = np.max([lpta_mean+20, lpta_median+20])
lrange = np.arange(ll, ul)
node_distances = np.zeros(len(lrange))

for ii, l in enumerate(lrange):
    node_distances[ii] = dnode(int(l), m=int(l/2))


delta_sph = np.abs(node_distances - mean_md)
l_suitable = lrange[np.argmin(delta_sph)]



fig, axs = plt.subplots(1,2, sharey=True, figsize=(7,5), gridspec_kw={'width_ratios':[2,1]})
plt.subplots_adjust(wspace=0.1)

axs[0].errorbar(lrange, node_distances, label=r'$d_\mathrm{SpH nodes}$')
axs[0].axhline(mean_md, ls=':', color='k', label=r'mean $d_\mathrm{min,PTA}$')
axs[0].axhline(median_md, ls='--', color='k', label=r'median $d_\mathrm{min,PTA}$')

axs[0].axvline(lpta_mean, ls=':', color='navy')
axs[0].axvline(lpta_median, ls='--', color='navy')

axs[1].hist(min_distances, orientation = 'horizontal', color='gray', bins=40, density=True)
axs[1].axhline(mean_md, ls=':', color='k')
axs[1].axhline(median_md, ls='--', color='k')


axs[0].set_xlim(ll,ul-1)
axs[0].set_ylim(0)
axs[0].legend(loc='best')

axs[0].set_xlabel('Sph. harmonics degree $l$')
axs[0].set_ylabel('angular distance / rad')

#plt.savefig('lmax_from_psr_separation.png', bbox_inches='tight', dpi=400)
plt.show()

print("lmax from mean of pulsar distances:\t {}".format(lpta_mean))
print("lmax from median of pulsar distances:\t {}".format(lpta_median))





