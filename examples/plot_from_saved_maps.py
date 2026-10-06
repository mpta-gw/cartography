# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 17:03:38 2026

@author: kgrun
"""


import sys, argparse
import pickle, json
import numpy as np

from PTA_AnalysisUtils import PTA_utils
from cartography import plot_maps as pm



def set_up_global_options():
    parser = argparse.ArgumentParser(description='Anisotropy analysis')
    
    #-  data arguments
    parser.add_argument('--datapath', type=str, default=None,
                        help='Path to the directory with par and tim files')
    
    parser.add_argument('--psrlist', type=str, default=None,
                        help='file with list of pulsars to use')
    
    parser.add_argument('--psrpickle', type=str, default=None,
                        help='pickle file with PSRs object')
    
    
    #-  anisotropy arguments  -------------------------------------------------
    parser.add_argument('--fbin', type=int, default=1,
                        help='Frequncy bin')
    
    parser.add_argument('--lmax', type=int, default=6,
                        help='maximum spherical harmonics')
    
    parser.add_argument('--nside', type=int, default=16,
                        help='healpix NSIDE value')
    
    parser.add_argument('--cutoff', type=int, default=None,
                        help='cutoff index for regularisation')
    
    #-  performance arguments  ------------------------------------------------
    parser.add_argument('--basepath', type=str, default=None,
                        help='path to the directory, in which the analysis files are stored')
    
    #-  output  ---------------------------------------------------------------
    parser.add_argument('--outdir', type=str, default=None,
                        help='Directory for the output')
    
    return parser.parse_args()




def main():
    args = set_up_global_options()
    sys.path.append(args.basepath + '/PTA_AnalysisUtils/')

    
    print("\n\n\n")
    #==========================================================================
    print('LOAD DATA INTO PTA MODEL', '\n')
    
    
    # --- load pulsars --------------------------------------------------------
    
    psr_list = list(np.loadtxt(args.psrlist, dtype='str'))
    
    
    if args.psrpickle == None:
        print('... loading pulsars from files')
        PSRs = PTA_utils.load_pulsars(args.datapath + '/pars/',
                                     args.datapath + '/tims/',
                                     psr_list)
         
        with open(args.datapath+'/psrs.pkl', 'wb') as psrpickle:
            pickle.dump(PSRs, psrpickle)
        psrpickle.close()
    
    else:
        print('... loading pulsars from pickle')
        with open(args.psrpickle, 'rb') as psrpickle:
            allPSRs = pickle.load(psrpickle)
        psrpickle.close()
        PSRs = [p for p in allPSRs if p.name in psr_list]   #include this for dropout test
        
   
    # get positions and rms of pulsars
    rms_dict = {}
    
    for pp, psr in enumerate(PSRs):
        rms_dict[psr.name] =  np.sqrt(np.mean(psr.residuals**2))
        
        
    setup_dict = {"NSIDE": args.nside,
                  "lmax": args.lmax,
                  "cutoff": args.cutoff,
                  "fbin": args.fbin,
                  "res": 1
                  }
    
    # =========================================================================    
    print('LOAD FILES FROM MAPS TO {}\n'.format(args.outdir))
    
    
    figuredir = args.outdir + "/figures/"
    
    
    out_maps = np.loadtxt(args.outdir + '/maps_{}_{}_{}_fbin{}.txt'.format(setup_dict["lmax"], setup_dict["cutoff"], setup_dict["NSIDE"], setup_dict["fbin"]))
    maps_array = out_maps.transpose() 
    sn_radio_map, sn_map, sn_map_nmp = pm.plot_skymaps(maps_array, PSRs, setup_dict, rms_dict, figuredir,
                                                       plot_style = "paper_large", cmap='plasma')
    
    
    calculate_pvalue = True
    if calculate_pvalue is True:
        distr_array = np.loadtxt(args.outdir + '/distribution_maxmin_{}_{}_{}_fbin{}.txt'.format(setup_dict["lmax"], setup_dict["cutoff"], setup_dict["NSIDE"], setup_dict["fbin"]) )
        pm.obtain_pvalues(distr_array, setup_dict, figuredir, sn_radio_map, sn_map, sn_map_nmp)
    # =========================================================================





if __name__ == '__main__':
    main()

        
