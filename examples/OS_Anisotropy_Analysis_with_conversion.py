# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 15:02:36 2026

@author: kgrun
"""



import sys, argparse
import pickle, json

import numpy as np


from defiant import OptimalStatistic

from PTA_AnalysisUtils import PTA_utils
from PTA_AnalysisUtils.IPTA_model import PTA
from cartography import anisotropy_analysis as aa
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
    
    
    parser.add_argument('--noisefile', type=str, default=None,
                        help='json file with noise parameters')
    
    parser.add_argument('--noisemodels', type=str, default=None,
                        help='json file with SPNA model specifications')
    
    #-  anisotropy arguments  -------------------------------------------------
    parser.add_argument('--fbin', type=int, default=1,
                        help='Frequncy bin')
    
    parser.add_argument('--lmax', type=int, default=6,
                        help='maximum spherical harmonics')
    
    parser.add_argument('--nside', type=int, default=16,
                        help='healpix NSIDE value')
    
    parser.add_argument('--cutoff', type=int, default=None,
                        help='cutoff index for regularisation')

    parser.add_argument('--incCrossCorr', action='store_true',
                        help='include cross correlations in OS sigma')
    
    #-  performance arguments  ------------------------------------------------
    parser.add_argument('--basepath', type=str, default=None,
                        help='path to the directory, in which the analysis files are stored')
    
    parser.add_argument('--create_plots', action='store_true',
                        help='create output plots')
    
    parser.add_argument('--check_iso', action='store_true',
                        help='compare c_00 to isotropic OS analysis')
    
    

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
    psrs_theta, psrs_phi = np.zeros(len(PSRs)), np.zeros(len(PSRs))
    rms_dict = {}
    
    for pp, psr in enumerate(PSRs):
        psrs_theta[pp] = psr.theta
        psrs_phi[pp] = psr.phi
        rms_dict[psr.name] =  np.sqrt(np.mean(psr.residuals**2))
        

    # --- load noise dictionary -----------------------------------------------
    print('... loading noise dictionary')
    with open(args.noisefile, 'r') as nf:
        noisedict = json.load(nf)
    nf.close()   
    
    

    # --- create PTA object ---------------------------------------------------
    print('... setting up the PTA model (loading SPNA models and common signals)')
    # load the individual pulsar models
    with open(args.noisemodels, 'r') as nm:
        noisemodel_dict = json.load(nm)
    nm.close()
        
        
    print('... creating the PTA')
    print('... pulsars and noise components: \n')
    pta = PTA(PSRs, noisemodel_dict, noisedict)
    pta.set_default_params(noisedict)

    
    # --- check noisefile -----------------------------------------------------
    checkpar = []
    for parameter in pta.param_names:
        if parameter not in noisedict.keys():
            checkpar.append(parameter)
    
    print('... check on the following parameters: ')
    print(checkpar, '\n')
    #==========================================================================
    
    

    #==========================================================================    
    print('\n' + 50*'-' + '\n' + 50*'-' + '\n',
          'OS ANALYSIS')

    print('... build OS object')
    os_obj = OptimalStatistic(PSRs, pta, gwb_name='gw')
    os_obj.set_orf(['hd'])

    if "gw_log10_A" not in noisedict.keys():
        print('... log10A not given in noisedict, obtain from OS')
        A2, _, _ = os_obj.compute_OS(params = noisedict, return_pair_vals = False)
        log10_Agw = np.log10(np.sqrt(A2))
        print('    log10A = {}'.format(log10_Agw)) 
        noisedict['gw_log10_A'] = log10_Agw
        pta.set_default_params(noisedict)

    # helpful output
    with open(args.outdir+'/summary_fbin{}.txt'.format(args.fbin), 'w') as sumfile:
        sumfile.write(pta.summary())
    sumfile.close()
    #==========================================================================
        
     
     
    
        
    #==========================================================================
    print('\n' + 50*'-' + '\n' + 50*'-' + '\n',
          'ANISOTROPY ANALYSIS')
    
    setup_dict = {"NSIDE": args.nside,
                  "lmax": args.lmax,
                  "cutoff": args.cutoff,
                  "fbin": args.fbin,
                  "res": 1
                  }

    clms, clms_err, X, M, cov_M, Mprime_inv = aa.MLsolution_with_reg(PSRs, os_obj, noisedict, setup_dict, args.outdir,
                                                                     narrowband = True, incCrossCorr = args.incCrossCorr,
                                                                     check_iso = False, plot_sv = args.create_plots)
                                                                     

    calculate_pvalue = True
    if calculate_pvalue is True:
        out_maps_array, out_fisher, distr_array = aa.convert_maps_fast(clms, M, X, cov_M, Mprime_inv, setup_dict, args.outdir,
                                                                       calculate_pvalue = True)
    else:
        out_maps_array, out_fisher = aa.convert_maps_fast(clms, M, X, cov_M, Mprime_inv, setup_dict, args.outdir,
                                                          calculate_pvalue = False)

    #==========================================================================
    

    # == Plot =================================================================
    print("PLOTTING")
    # data structure maps
    # 0 - RA in hours
    # 1 - DEC in degree [-90, 90]
    # 2 - clean map
    # 3 - sigma map
    # 4 - S/N map
    # 5 - radiometer map
    # 6 - sensitivity map
    
    # data structure distribution
    # 0 - max radiometer snr
    # 1 - min radiometer snr
    # 2 - max clean snr
    # 3 - min clean snr
    # 4 - max clean snr w/o monopole
    # 5 - min clean snr w/o monopole
    
    
    maps_array = out_maps_array.transpose() 
    figuredir = args.outdir + "/figures/"
    sn_radio_map, sn_map, sn_map_nmp = pm.plot_skymaps(maps_array, PSRs, setup_dict, rms_dict, figuredir,
                                                       plot_style = "paper_large", cmap='plasma')

    if calculate_pvalue is True:
        pm.obtain_pvalues(distr_array, setup_dict, figuredir, sn_radio_map, sn_map, sn_map_nmp)
                
        
    return None





if __name__ == '__main__':
    main()


