#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 16:19:40 2024

@author: kgrunthal
"""

import os
import numpy as np
import matplotlib.pyplot as plt





###############################################################################
#   font and label settings 

plt.rcParams['font.family'] = "serif"
plt.rcParams['font.sans-serif'] = "Times"

plt.rcParams['text.usetex']= False
plt.rcParams['xtick.labelsize'] = 10.0
plt.rcParams['ytick.labelsize'] = 10.0
plt.rcParams['axes.labelsize'] = 20.0


def get_normalisation(map_values, resolution=1, angle_mode='deg'):
    if angle_mode == 'deg':
        sqdeg2sqrad = (np.pi/180.)**2
        res = resolution * sqdeg2sqrad
    elif angle_mode =='rad':
        res = resolution
    else:
        print('Unsupported angle format')
        return None
    pixsum = np.sum(map_values)
    normalisation = res/pixsum
    print('map sum: \t', pixsum)
    print('resolution: \t', res)
    print('normalisation: \t', normalisation) 
    return res
        

def convert2healpix(map_values, pxs, NPIX):
    '''
    Convert a skymap with arbitrary (theta, phi) distribution to a healpix 
    pixel distribution.
    
    Calculates mean value in each pixel
    
    param map_values:   skymap
    param pxs:          healpix pixel indices for the map
    param NPIX:         number of healpix pixels
    '''
    
    hpmap = np.zeros(NPIX)
    added = np.zeros(NPIX)
    
    for ii, px in enumerate(pxs):
        hpmap[px] += map_values[ii]
        added[px] += 1
    
    return hpmap/added
    
    
def transform_phi(phi):
    #if phi < np.pi:
    #    return -1.*phi
    #else:
    #    return 2*np.pi-phi
    return np.pi-phi

def plot_thy_pulsars(ax, PSRs, rms_dict, color='black', pulsar_names = True):

    # sort pulsars by RMS
    s_val = []
    for psr in PSRs:
        s = 28 * (1 / rms_dict[psr.name])
        #print('{}\t{}\t{}\t{}'.format(psr.name, psr.phi, psr.theta, rms_dict[psr.name]))
        s_val.append(s)
        ax.scatter(transform_phi(psr.phi), (np.pi/2-psr.theta),
                   color = color, marker='*', s=s,
                   edgecolors='white', linewidths=0.3)
    
    if pulsar_names is True:
        zipped = zip(s_val, PSRs)
        pulsar_importance = sorted(zipped, key = lambda x: x[0], reverse=True)
    
        # assert names to the seven best pulsars
        for jj in range(7):
            s_j, psr_j = pulsar_importance[jj]
            # plot the pulsars
            ax.text(transform_phi(psr_j.phi-0.1), (np.pi/2-psr_j.theta+0.05), psr_j.name,
                    color = 'white', fontsize=5)
    
    else:
        None
        
    return None


def crosshair(ax, phi, theta, size=6):
    color='dimgray'
    ax.plot(phi, theta, ls='', marker='o', mfc='none', color=color, ms=size, mew=1.5)
    ax.plot(phi, theta, ls='', marker='+', mfc='none', color=color, ms=2*size, mew=1.5)
    
    return None

def plotmap(map_values, PSRs, rms_dict, fbin, cutoff, lmax, nside, res,
            figuredir = '', normalisation = None,
            maptype = '', plot_style = 'detailed',
            mark_max = False, injected = None,
            sn_range=None, sv_value=False, cmap = 'plasma'):

    
    label = {'clean': r"$\tilde{\mathcal{P}}'(\hat{\Omega})$ / sr$^{-1}$",
             'clean_norm' : r"$\tilde{\mathcal{P}}'(\hat{\Omega})$",
             'sigma' : r"$\left[\text{Var}(\tilde{\mathcal{P}}'(\hat{\Omega}) )\right]^{-1/2}$ / sr",
             'sigma_norm': r"$\text{Var}(\tilde{\mathcal{P}}'(\hat{\Omega}))^{-1/2}$",
             'sn': '$S/N_\mathrm{clean}$',
             'sn_radio': '$S/N_\mathrm{radiometer}$',
             'sensitivity': 'sensitivity'}
    
    title = {'clean': 'Clean map',
             'clean_norm': 'Clean map',
             'sigma': 'Sensitivity map',
             'sigma_norm': 'Sensitivity map',
             'sn': 'S/N map',
             'sn_radio': 'Radiometer map',
             'sensitivity': 'Sensitivity map'}
    
    name = {'clean': 'cleanmap',
            'clean_norm': 'cleanmap_norm',
            'sigma': 'sensitivitymap',
            'sigma_norm': 'sensitivitymap_norm',
            'sn': 'snrmap',
            'sn_radio': 'radiometer',
            'sensitivity': 'sensitivity'}
    
    subtitle = '$f = $'+ '{}'.format(fbin)+'$/T_\mathrm{obs}$ \n' + \
                'resolution: $A_\mathrm{px} = $' +'{} '.format(res) + 'deg$^2$, ' + \
                '$l_\mathrm{max} = $' + '{}\n'.format(lmax) + \
                'singular value cutoff: $n_{s_v} = $' + '{}'.format(cutoff)
    subtitle_radio = '$f = $'+ '{}'.format(fbin)+'$/T_\mathrm{obs}$ \n' + \
                'resolution: $A_\mathrm{px} = $' +'{} '.format(res) + 'deg$^2$, ' + \
                '$l_\mathrm{max} = $' + '{}\n'.format(lmax)
            
    
    folder = {'clean': 'clean',
              'clean_norm': 'clean_norm',
              'sigma': 'sigma',
              'sigma_norm': 'sigma_norm',
              'sn': 'snr',
              'sn_radio': 'radiometer',
              'sensitivity': 'sensitivity'
             }
  
    
    ax = plt.axes([0.05, 0.05, 0.9, 0.9], projection='mollweide')
    ax.grid(linestyle='dotted', color='black')
   

    # arrange map amd mesh --> flip map according to PTA convention  
    #  - x-axis is reversed (flip the array along 0 axis)
    
  
    lon = np.linspace(-np.pi, np.pi, num=360)
    lat = np.linspace(-np.pi/2., np.pi/2., num=181, endpoint=True)
    LAT, LON = np.meshgrid(lat,lon)

    map_racorr = np.flip(map_values, axis=0)

    
    # set up the grid
    ax.set_latitude_grid(30)
    ax.set_longitude_grid(45)
    ax.set_longitude_grid_ends(90)
    ax.set_xticks(np.linspace(-np.pi, np.pi, num=9),
                  ['', '21$^\mathregular{h}$', '18$^\mathregular{h}$',
                   '15$^\mathregular{h}$', '12$^\mathregular{h}$', '9$^\mathregular{h}$',
                   '6$^\mathregular{h}$', '3$^\mathregular{h}$', '0$^\mathregular{h}$'])

    ax.grid(linestyle='dotted', color='black')

    ## plot the skymap ########################################################
    # sn range fixing
    if (maptype == 'sn' or maptype == 'sn_radio') and sn_range is not None:
        cax = ax.pcolormesh(LON, LAT, map_racorr, cmap=cmap, vmin=-1*sn_range, vmax=sn_range)
        
    else:
        cax = ax.pcolormesh(LON, LAT, map_racorr, cmap=cmap)
        #cax = ax.pcolormesh(LON, LAT, map_radeccorr, cmap='plasma')
        
    # print sv value
    if sv_value == True:
        plt.text(4,1,'$i = ${}'.format(cutoff))
    
    #mark hot spot on S/N map
    if maptype == 'sn' and mark_max == True:
        lon_idx, lat_idx = np.where(map_values == map_values.max())
        lon_val = lon[lon_idx]
        lat_val = lat[lat_idx]
        ax.errorbar(lon_val, lat_val, ls ='', fmt='Xb', ms=5)
    

    ## mark injected value on S/N map ########################################

    #crosshair(ax, np.pi/2, np.pi/4., size=11)		# RA06h DEC45deg
    #crosshair(ax, np.pi/2, 0., size=11)                # RA06h DEC00deg
    #crosshair(ax, np.pi/2, -np.pi/4., size=11)       	# RA06h DEC-45deg
    
    #crosshair(ax, 0., np.pi/4., size=11)		# RA12h DEC45deg
    #crosshair(ax, 0., 0., size=11)               # RA12h DEC00deg
    #crosshair(ax, 0., -np.pi/4., size=11)                # RA12h DEC-45deg
 
    #crosshair(ax, -np.pi/2, np.pi/4., size=11)         # RA18h DEC45deg
    #crosshair(ax, -np.pi/2, 0., size=11)               # RA18h DEC00deg
    #crosshair(ax, -np.pi/2, -np.pi/4., size=11)        # RA18h DEC-45deg
    
    
    ## output style dependent plot finish #####################################
    if plot_style == 'detailed':
        # colorbar
        cbar = plt.colorbar(cax, orientation='horizontal',pad=0.05, label=label[maptype])
        cbar.ax.tick_params(labelsize=18)

        # plot the pulsars
        plot_thy_pulsars(ax, PSRs, rms_dict, color='white', pulsar_names=True)    
        
        plt.suptitle(title[maptype], fontsize=20,y=1.02)
        if maptype != 'sn_radio':
            plt.title(subtitle, fontsize=10)
        else:
            plt.title(subtitle_radio, fontsize=10)

        plt.savefig(figuredir + '/detailed/{}/{}_lmax{}_sv{}_nside{}_fbin{}.png'.format(folder[maptype], name[maptype], lmax, cutoff, nside, fbin),
                    bbox_inches='tight', dpi=400)
        
    
    elif plot_style == 'paper_large':
        os.makedirs(figuredir + '/paper_large/', exist_ok=True)
        os.makedirs(figuredir + '/paper_large/{}'.format(folder[maptype]), exist_ok=True)
        # colorbar
        cbar = plt.colorbar(cax, orientation='horizontal',pad=0.05, label=label[maptype])
        cbar.ax.tick_params(labelsize=18)

        # plot the pulsars
        plot_thy_pulsars(ax, PSRs, rms_dict, color='white', pulsar_names=True)
        
        if sv_value == True and folder[maptype] == 'snr':
            plt.savefig(figuredir + '/movie/{}/{}_paperlarge_lmax{}_sv{}_nside{}_fbin{}.png'.format(folder[maptype], name[maptype], lmax, cutoff, nside, fbin), bbox_inches='tight', dpi=400)
        else:
            plt.savefig(figuredir + '/paper_large/{}/{}_paperlarge_lmax{}_sv{}_nside{}_fbin{}.png'.format(folder[maptype], name[maptype], lmax, cutoff, nside, fbin), bbox_inches='tight', dpi=400)

    elif plot_style == 'paper_small':
        # colorbar
        cbar = plt.colorbar(cax, orientation='horizontal',pad=0.05, label=label[maptype])
        cbar.ax.tick_params(labelsize=18)

        # plot the pulsars
        plot_thy_pulsars(ax, PSRs, rms_dict, color='white', pulsar_names=False)
        
        plt.savefig(figuredir + '/paper_small/{}/{}_papersmall_lmax{}_sv{}_nside{}_fbin{}.png'.format(folder[maptype], name[maptype], lmax, cutoff, nside, fbin),
                    bbox_inches='tight', dpi=400)
    

    elif plot_style == 'paper_tiny':
        plt.rcParams['xtick.labelsize'] = 14.0
        plt.rcParams['ytick.labelsize'] = 14.0
        plt.rcParams['axes.labelsize'] = 22.0

        # colorbar
        cbar = plt.colorbar(cax, orientation='horizontal',pad=0.05, label=label[maptype])
        cbar.ax.tick_params(labelsize=20)

        # plot the pulsars
        plot_thy_pulsars(ax, PSRs, rms_dict, color='white', pulsar_names=False)

        plt.savefig(figuredir + '/paper_tiny/{}/{}_papertiny_lmax{}_sv{}_nside{}_fbin{}.png'.format(folder[maptype], name[maptype], lmax, cutoff, nside, fbin),
                    bbox_inches='tight', dpi=400)

    plt.clf()
    return None
    



def distribution_processing(distribution, skymap, fbin, cutoff, lmax, nside, res, figuredir='', maptype=''):
    """
    old function, for individual distribution files
    """

    os.makedirs(figuredir + '/distribution_plots/', exist_ok=True)

    max_snr0, min_snr0 = np.max(skymap), np.min(skymap)
    max_snr_dist, min_snr_dist = distribution[0], distribution[1]
    
    
    pvalue_max = len(np.where(max_snr_dist > max_snr0)[0]) / len(max_snr_dist)
    pvalue_min = len(np.where(min_snr_dist < min_snr0)[0]) / len(min_snr_dist)
    print('\t p(S/N,max) = {}'.format(pvalue_max))
    print('\t p(S/N,min) = {}'.format(pvalue_min))
    
    print('\t plot histograms')
    hist_bins=50

    plt.hist(min_snr_dist, bins=hist_bins, density=True,
             color='midnightblue', alpha=0.5,
             label='distribution $S/N_\mathrm{min}$')
    plt.hist(min_snr_dist, bins=hist_bins, density=True, histtype='step',
             color='midnightblue', lw=1)
    
    plt.hist(max_snr_dist, bins=hist_bins,
             color='goldenrod', alpha=0.5, density=True,
             label='distribution $S/N_\mathrm{max}$')
    plt.hist(max_snr_dist, bins=hist_bins, density=True, histtype='step',
             color='goldenrod', lw=1)

    plt.axvline(min_snr0, color='navy', ls=':', label='$S/N_\mathrm{min}^\mathrm{meas.}$')
    plt.axvline(max_snr0, color='darkorange', ls=':', label='$S/N_\mathrm{max}^\mathrm{meas.}$')
    
    plt.xlabel('$S/N$', fontsize = 18)
    plt.legend(bbox_to_anchor=(1.05,1.0))
    plt.savefig(figuredir + '/distribution_plots/distribution_{}_lmax{}_sv{}_nside{}_fbin{}.png'.format(maptype, lmax, cutoff, nside, fbin), bbox_inches='tight', dpi=400)
    plt.clf()
    



def obtain_pvalues(distr_array, setup_dict, figuredir,
                   sn_radio_map, sn_map, sn_map_nmp):
    
    # data structure distribution
    # 0 - max radiometer snr
    # 1 - min radiometer snr
    # 2 - max clean snr
    # 3 - min clean snr
    # 4 - max clean snr w/o monopole
    # 5 - min clean snr w/o monopole
    
    print('CALCULATE P-VALUES')
    
    print('... radiometer S/N')
    distribution_processing(distr_array[:,0:2], sn_radio_map,
                            setup_dict["fbin"], setup_dict["cutoff"], setup_dict["lmax"], setup_dict["NSIDE"], setup_dict["res"],
                            figuredir=figuredir, maptype='radiometer')
    print()
    print('... clean S/N')
    distribution_processing(distr_array[:,2:4], sn_map,
                            setup_dict["fbin"], setup_dict["cutoff"], setup_dict["lmax"], setup_dict["NSIDE"], setup_dict["res"],
                            figuredir=figuredir, maptype='clean')
    print()
    print('... clean S/N w/o monopole')
    distribution_processing(distr_array[:,4:], sn_map_nmp,
                            setup_dict["fbin"], setup_dict["cutoff"], setup_dict["lmax"], setup_dict["NSIDE"], setup_dict["res"],
                            figuredir=figuredir, maptype='clean')
    
    return None
    



    

    



def plot_skymaps(data_map, PSRs, setup_dict, rms_dict, figuredir,
                 plot_style = "paper_large", cmap='plasma',
                 sn_range = None, mark_max=False):   
    """
    SUMMARY.

    Parameters
    ----------
    data_map : TYPE
        DESCRIPTION.
    PSRs : TYPE
        DESCRIPTION.
    setup_dict : TYPE
        DESCRIPTION.
    rms_dict : TYPE
        DESCRIPTION.
    plot_style : TYPE, optional
        DESCRIPTION. The default is "paper_large".
        choices=['detailed', 'paper_large', 'paper_small', 'paper_tiny'],
        

    Returns
    -------
    int
        DESCRIPTION.
    """
    
    RA = data_map[0].reshape(360,181)
    DEC = data_map[1].reshape(360,181)

    ## reshape the sky maps for plotting ######################################
    print('RESHAPE MAPS')
    print('... clean map')
    clean_map = data_map[2].reshape(360,181)
    print(clean_map)
    print('... sigma map')
    sigma_map = data_map[3].reshape(360,181)
    print('... S/N map')
    sn_map = data_map[4].reshape(360,181)
    print('... S/N radiometer map')
    sn_radio_map = data_map[5].reshape(360,181)
    print('... sensitivity map')
    sensitivity_map = data_map[6].reshape(360,181)
    print('... clean map w/o monopole')
    sn_map_nmp = data_map[7].reshape(360,181)
    print()


    
    idx = np.where(data_map[4] == data_map[4].max())
    RA_val = data_map[0][idx]
    DEC_val = data_map[1][idx]
    print('\t Maximum of ', np.round(np.max(sn_map), 2), 'is at RA  = ', RA_val, 'DEC = ', DEC_val)
   
 
    ## clean map normalisation ################################################
    #normalisation = get_normalisation(clean_map)
    
    ## plot maps ##############################################################
    print('PLOT MAPS')
    
    print('... clean map')
    plotmap(clean_map, PSRs, rms_dict,
            setup_dict["fbin"], setup_dict["cutoff"], setup_dict["lmax"], setup_dict["NSIDE"], setup_dict["res"],
            maptype='clean', figuredir = figuredir, plot_style = plot_style, cmap=cmap)
    #plotmap(clean_map*normalisation, PSRs,
    #        args.fbins, args.cutoff, args.lmax, args.nside, args.res,
    #        maptype='clean_norm', figuredir=figuredir, plot_style=args.plot_style)
    
    print('... sigma map')
    plotmap(1/sigma_map, PSRs, rms_dict,
            setup_dict["fbin"], setup_dict["cutoff"], setup_dict["lmax"], setup_dict["NSIDE"], setup_dict["res"],
            maptype='sigma', figuredir=figuredir, plot_style= plot_style, cmap=cmap)
    #plotmap((1./sigma_map)/normalisation, PSRs,
            #args.fbins, args.cutoff, args.lmax, args.nside, args.res,
            #maptype='sigma_norm', figuredir=figuredir, plot_style=args.plot_style)
    
    print('... S/N clean map')
    plotmap(sn_map, PSRs, rms_dict,
            setup_dict["fbin"], setup_dict["cutoff"], setup_dict["lmax"], setup_dict["NSIDE"], setup_dict["res"],
            maptype='sn', figuredir=figuredir, plot_style = plot_style,
            mark_max = mark_max, sn_range = sn_range, cmap=cmap)
    
    print('... S/N radiometer map')
    plotmap(sn_radio_map, PSRs, rms_dict,
            setup_dict["fbin"], setup_dict["cutoff"], setup_dict["lmax"], setup_dict["NSIDE"], setup_dict["res"],
            maptype='sn_radio', figuredir=figuredir, plot_style=plot_style,
            sn_range = sn_range, cmap=cmap)

    print('... sensitivity map')
    plotmap(sensitivity_map, PSRs, rms_dict,
            setup_dict["fbin"], setup_dict["cutoff"], setup_dict["lmax"], setup_dict["NSIDE"], setup_dict["res"],
            maptype='sensitivity', figuredir=figuredir, plot_style = plot_style, cmap=cmap)

    print()
    

    
    return sn_map, sn_radio_map, sn_map_nmp



    
