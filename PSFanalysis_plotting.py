# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 16:24:39 2024

@author: kgrun
"""


import numpy as np
import matplotlib.pyplot as plt
import colormaps as cmaps
import argparse
import pickle,json
import astropy

###############################################################################
#   font and label settings 

plt.rcParams['font.family'] = "serif"
plt.rcParams['font.sans-serif'] = "Times"

plt.rcParams['text.usetex']= False
plt.rcParams['xtick.labelsize'] = 10.0
plt.rcParams['ytick.labelsize'] = 10.0
plt.rcParams['axes.labelsize'] = 20.0




 


def plot_thy_pulsars(ax, PSRs, rms_dict, color='black', pulsar_names = True):

    # sort pulsars by RMS
    s_val = []
    for psr in PSRs:
        s = 10 * (1 / rms_dict[psr.name])
        s_val.append(s)
        ax.scatter(np.pi - psr.phi, (np.pi/2-psr.theta),
                   color = color, marker='*', s=s,
                   edgecolors='white', linewidths=0.3)
    
    if pulsar_names is True:
        zipped = zip(s_val, PSRs)
        pulsar_importance = sorted(zipped, key = lambda x: x[0], reverse=True)
    
        # assert names to the seven best pulsars
        for jj in range(7):
            s_j, psr_j = pulsar_importance[jj]
            # plot the pulsars
            ax.text(np.pi - (psr_j.phi-0.1), (np.pi/2-psr_j.theta+0.05), psr_j.name,
                    color = 'white', fontsize=5)
    
    else:
        None
        
    return None



def plt_mw(LON, LAT, mapvalues, labeltext='', outfile='', pplot = None):
    ax = plt.axes([0.05, 0.05, 0.9, 0.9], projection='mollweide')

    # map
    cax = ax.pcolormesh(LON, LAT, np.flip(mapvalues, axis = 0) ,cmap=cmaps.ice_r.shift(0.15)) 
    cbar = plt.colorbar(cax, orientation='horizontal', label = labeltext)
    cbar.ax.tick_params(labelsize=18)

    # pulsars
    if pplot is not None:
        plot_thy_pulsars(ax, pplot['PSRs'], pplot['rms_dict'], color='white', pulsar_names=pplot['PN'])

    # grid
    ax.set_latitude_grid(30)
    ax.set_longitude_grid(45)
    ax.set_longitude_grid_ends(90)
    ax.set_xticks(np.linspace(-np.pi, np.pi, num=9),
                  ['', '21$^\mathregular{h}$', '18$^\mathregular{h}$',
                   '15$^\mathregular{h}$', '12$^\mathregular{h}$', '9$^\mathregular{h}$',
                   '6$^\mathregular{h}$', '3$^\mathregular{h}$', '0$^\mathregular{h}$'])

    ax.grid(linestyle='dotted', color='black')

    plt.savefig(outfile, bbox_inches='tight', dpi = 400)
    plt.clf()
    
    return None





def plot_mollweide_withpos(LON, LAT, mapvalues, pix, labeltext='', outfile='', pplot = None):
    ax = plt.axes([0.05, 0.05, 0.9, 0.9], projection='mollweide')

    # map
    cax = ax.pcolormesh(LON, LAT, np.flip(mapvalues, axis = 0) ,cmap=cmaps.ice_r.shift(0.15)) 
    cbar = plt.colorbar(cax, orientation='horizontal', label = labeltext)
    cbar.ax.tick_params(labelsize=18)

    # pulsars
    if pplot is not None:
        plot_thy_pulsars(ax, pplot['PSRs'], pplot['rms_dict'], color='white', pulsar_names=pplot['PN'])

    #plot position
    idx = np.unravel_index(pix, mapvalues.shape)
    phi = np.linspace(-np.pi, np.pi, num = mapvalues.shape[0])[idx[0]]
    theta = np.linspace(-np.pi/2., np.pi/2., num=mapvalues.shape[1], endpoint=True)[idx[1]]

    color='k'
    size=7
    ax.plot(phi, theta, ls='', marker='o', mfc='none', color=color, ms=size, mew=1)
    ax.plot(phi, theta, ls='', marker='+', mfc='none', color=color, ms=2*size, mew=1)
    
    
    
    # grid
    ax.set_latitude_grid(30)
    ax.set_longitude_grid(45)
    ax.set_longitude_grid_ends(90)
    ax.set_xticks(np.linspace(-np.pi, np.pi, num=9),
                  ['', '21$^\mathregular{h}$', '18$^\mathregular{h}$',
                   '15$^\mathregular{h}$', '12$^\mathregular{h}$', '9$^\mathregular{h}$',
                   '6$^\mathregular{h}$', '3$^\mathregular{h}$', '0$^\mathregular{h}$'])

    ax.grid(linestyle='dotted', color='black')

    plt.savefig(outfile, bbox_inches='tight', dpi = 400)
    plt.clf()
    
    return None




def set_up_global_options():
    parser = argparse.ArgumentParser(description='perform PSF analysis')

    parser.add_argument('--head_dir', type=str, default=None,
                        help='Path to the directory with analysis output')
    
    parser.add_argument('--psrpickle', type=str, default=None,
                        help='Path to the directory with analysis output')

    parser.add_argument('--psrlist', type=str, default=None,
                        help='txt-file with psr names in the array')

    parser.add_argument('--tres', type=str, default=None,
                        help='pulsar RMS')

    parser.add_argument('--lmax', type=int, default=8,
                        help='maximum spherical harmonics')

    parser.add_argument('--nside', type=int, default=2,
                        help='healpix NSIDE value')

    parser.add_argument('--cutoff', type=int, default=None,
                        help='cutoff index for regularisation')

    parser.add_argument('--fbin', type=str,
                        help='denote which frequency bins have been used to calculate the map')


    parser.add_argument('--res', type=int, default=1,
                        help='map resolution')


    parser.add_argument('--level', type=float, default=0.5,
                        help='contour value to determine the PSF')

    parser.add_argument('--calculate', action='store_true',
                        help='perform matrix evaluation')

    parser.add_argument('--incCrossCorr', action='store_true',
                        help='data with or without crosscorrelations')

    return parser.parse_args()





args = set_up_global_options()


psr_list = list(np.loadtxt(args.psrlist, dtype='str'))

if args.incCrossCorr == True:
    which = 'cc'
else:
    which = 'nocc'

if args.calculate == True:

    print('loading pre-calculated matrices', flush=True)

    print('... MpinvM', flush=True)
    MpinvM_re = np.loadtxt(args.head_dir + '/maps/MpinvM-re_' + which + '_lmax{}_res{}_sv{}_fbin{}.txt'.format(args.lmax, args.res, args.cutoff, args.fbin))
    MpinvM_im = np.loadtxt(args.head_dir + '/maps/MpinvM-im_' + which + '_lmax{}_res{}_sv{}_fbin{}.txt'.format(args.lmax, args.res, args.cutoff, args.fbin))

    MpinvM = MpinvM_re + 1j*MpinvM_im


    print('... conversion matrix U', flush=True)
    U_re = np.loadtxt(args.head_dir + '/conversion_matrix/U-re_lmax{}_res{}.txt'.format(args.lmax, args.res))
    U_im = np.loadtxt(args.head_dir + '/conversion_matrix/U-im_lmax{}_res{}.txt'.format(args.lmax, args.res))

    U = U_re + 1j * U_im
    Uh = U.conj().transpose()


    print('... sensitivity map', flush=True)
    sigmaPix = np.loadtxt(args.head_dir + '/maps/maps_{}_{}_{}_{}_fbin{}.txt'.format(which, args.lmax, args.cutoff, args.nside, args.fbin), usecols = 3).transpose()

    print()
    ##########


    print('Setting up map grid', flush=True)
    # map specification
    res = args.res
    lons = int(360/res)
    lats = int(180/res) + 1

    #create grid
    lon = np.linspace(-np.pi, np.pi, num=lons)
    lat = np.linspace(-np.pi/2., np.pi/2., num=lats, endpoint=True)
    LAT, LON = np.meshgrid(lat,lon)


    print()
    ################

    print('preparing columnwise calculation of Lambda', flush=True)
    # number of pixel slices
    slices = len(U)

    UMpinvM = np.matmul(U, MpinvM)
    sigmaPix = sigmaPix.reshape((lons,lats))

    print()
    ####################################################


    print('calculate tesselation areas', flush=True)
    ## Calculate solid angles of the grid ##

    # longitude stripes
    avg_dlon = np.mean(np.diff(lon))
    dlon = np.append(np.diff(lon), 0.5*avg_dlon)
    dlon[0] = 0.5*avg_dlon
    
    # latitude stripes
    d_lat = np.zeros(len(lat))
                
    for i in range(len(lat)):
        if i == 0:
            lat_up = lat[i] + (lat[i+1] - lat[i])/2
            lat_low = lat[i]
        elif i == len(lat)-1:
            lat_up = lat[i]
            lat_low = lat[i] - (lat[i] - lat[i-1])/2
        else:
            lat_up = lat[i] + (lat[i+1] - lat[i])/2
            lat_low = lat[i] - (lat[i] - lat[i-1])/2
  
        d_lat[i] = np.sin(lat_up) - np.sin(lat_low)
      

    # Pixel area
    A_mesh = np.zeros((len(lon), len(lat)))
    for j in range(len(lon)):
        A_mesh[j] = dlon[j]*d_lat
    
    # Check that the area adds up to 4pi
    check_diff = np.sum(A_mesh) - 4*np.pi

    if np.abs(check_diff) < 1e-6:
        print('... area is calculated properly', flush=True)
    
    else:
        if check_diff < 0:
            print('... area is underestimated', flush=True)
        else:
            print('... area is overestimated', flush=True)
     

    print()
    #################################################


    # level for area calculation
    #levels = [np.exp(-1)]
    levels = [args.level]
    

    

    #iterate over all column vectors
    CLEANareas = np.zeros(slices)
    SNRareas = np.zeros(slices)

    for slc in range(slices):
        print('Slice {}/{}'.format(slc+1, slices), flush=True)
    
        # calculate, reshape and normalize map
        real_Lambda = np.real( np.dot(UMpinvM, Uh[:,slc]).reshape((lons,lats)) )
 
        CLEANmap_norm = real_Lambda/np.max(real_Lambda)

        snr = real_Lambda / sigmaPix
        SNRmap_norm = snr/np.max(snr)

    
        ## AREA OF THE REC HOTSPOT ##
        SNRindx = np.where(SNRmap_norm > args.level)
        CLEANindx = np.where(CLEANmap_norm > args.level)


        CLEANarea_temp = 0
        for i in range(len(CLEANindx[0])):
            CLEANidx_lo = CLEANindx[0][i]
            CLEANidx_la = CLEANindx[1][i]
            CLEANarea_temp += A_mesh[CLEANidx_lo,CLEANidx_la]
        CLEANareas[slc]=CLEANarea_temp
        del CLEANarea_temp    
    

    CLEAN_A_MAP = CLEANareas.reshape(lons,lats)
    
    np.savetxt(args.head_dir + '/maps/clean_PSFarea_{}_level{}_res{}_lmax{}_cut{}_nside{}_fbin{}.txt'.format(which, args.level, res, args.lmax, args.cutoff, args.nside, args.fbin), CLEAN_A_MAP)
    
 
    
else:
    print('Loading maps from files')
    CLEAN_A_MAP = np.loadtxt(args.head_dir + '/maps/clean_PSFarea_{}_level{}_res{}_lmax{}_cut{}_nside{}_fbin{}.txt'.format(which, args.level, args.res, args.lmax, args.cutoff, args.nside, args.fbin) )
    
    print('Setting up map grid', flush=True)
    # map specification
    res = args.res
    lons = int(360/res)
    lats = int(180/res) + 1

    #create grid
    lon = np.linspace(-np.pi, np.pi, num=lons)
    lat = np.linspace(-np.pi/2., np.pi/2., num=lats, endpoint=True)
    LAT, LON = np.meshgrid(lat,lon)
    print()


    


#-- PLOTTING ------------------------------------------------------------------

# Pulsars
if args.psrpickle is not None:
    print('Loading pulsars to plot')

    with open(args.psrpickle, 'rb') as psrpickle:
        allPSRs = pickle.load(psrpickle)
    psrpickle.close()
    PSRs = [p for p in allPSRs if p.name in psr_list]   #include this for dropout test

    with open(args.tres, 'r') as tresjson:
        rms_dict_temp = json.load(tresjson)
        if type(rms_dict_temp[list(rms_dict_temp.keys())[0]]) is list:
            rms_dict = {}
            for k in rms_dict_temp.keys():
                rms_dict[k] = rms_dict_temp[k][2]
        else:
            rms_dict = rms_dict_temp
    tresjson.close()
    
    pulsarplot_dict = {'PSRs': PSRs, 'rms_dict': rms_dict, 'PN': False}

else:
    pulsarplot_dict = None
    print('Plotting without pulsars')
        
    


    
### distortion size map ###

# CLEAN 
plt_mw(LON, LAT, CLEAN_A_MAP,
       labeltext = '$A$ / sr',
       outfile = args.head_dir + '/figures/PSFmap_{}_lmax{}_cutoff{}_res{}_level{}.png'.format(which, args.lmax, args.cutoff, res, args.level),
       pplot = pulsarplot_dict)






