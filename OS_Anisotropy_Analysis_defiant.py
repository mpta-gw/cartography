# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 17:03:59 2024

@author: kgrunthal
"""


import sys, argparse
import pickle, json

import numpy as np
import matplotlib.pyplot as plt


from defiant import OptimalStatistic

from PTA_AnalysisUtils import PTA_utils
from PTA_AnalysisUtils.PTA_model import PTA
from PTA_AnalysisUtils.maps_mod import anis_pta as anis_pta






def set_up_global_options():
    parser = argparse.ArgumentParser(description='Anisotropy analysis')
    
    #-  data arguments
    parser.add_argument('--datapath', type=str, default=None,
                        help='Path to the directory with  par and tim files')
    
    parser.add_argument('--psrlist', type=str, default=None,
                        help='file with list of pulsars to use')
    
    parser.add_argument('--psrpickle', type=str, default=None,
                        help='pickle file with PSRs object')
    
    parser.add_argument('--ptaobject', type=str, default=None,
                        help='pickle file with PTA object')
    
    parser.add_argument('--noisefile', type=str, default=None,
                        help='json file with noise parameters')
    
    parser.add_argument('--SPNAmodels', type=str, default=None,
                        help='json file with SPNA model specifications')
    
    parser.add_argument('--ptamodel', type=str, default=None,
                        help='Which signals to include in the PTA model')
       
    
    #-  anisotropy arguments  -------------------------------------------------
    parser.add_argument('--fbin', type=int, default=1,
                        help='Frequncy bin')
    
    parser.add_argument('--lmax', type=int, default=6,
                        help='maximum spherical harmonics')
    
    parser.add_argument('--nside', type=int, default=2,
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
    
    figuredir = args.outdir + '/figures/OS_Anisotropy_diagnostics/'
    basis = 'pbi'
    NSIDE = args.nside
    
    
    #==========================================================================
    print('LOAD DATA INTO PTA MODEL', '\n')
    
    
    # --- load pulsars --------------------------------------------------------
    
    psr_list = list(np.loadtxt(args.psrlist, dtype='str'))
    
    
    if args.psrpickle == None:
        print('... loading pulsars from files')
        PSRs = PTA_utils.load_pulsars(args.datapath + '/par/',
                                     args.datapath + '/tim/',
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
    
   
    # get positions of pulsars
    psrs_theta, psrs_phi = np.zeros(len(PSRs)), np.zeros(len(PSRs))
    for pp, psr in enumerate(PSRs):
        psrs_theta[pp] = psr.theta
        psrs_phi[pp] = psr.phi
        

    # --- load noise dictionary -----------------------------------------------
    print('... loading noise dictionary')
    with open(args.noisefile, 'r') as nf:
        noisedict = json.load(nf)
    nf.close()   
    
    

    # --- create PTA object ---------------------------------------------------
    if args.ptaobject is None:
        print('... setting up the PTA model')
        
        print('       - loading SPNA models')
        # load the individual pulsar models
        with open(args.SPNAmodels, 'r') as nm:
            spnamodel_dict = json.load(nm)
        nm.close()
        
        print('       - loading common signals')
        with open(args.ptamodel, 'r') as ptam:
            common_signal = json.load(ptam)
        ptam.close()
        
        print('... creating the PTA')
        pta = PTA(PSRs, noisedict, spnamodel_dict, common_signal)
      
    
    else:
        with open(args.ptaobject, 'rb') as ptaobj:
            pta = pickle.load(ptaobj)
        ptaobj.close()
        

    pta.set_default_params(noisedict)

    # helpful output
    with open(args.outdir+'/summary_fbin{}.txt'.format(args.fbin), 'w') as sumfile:
        sumfile.write(pta.summary())
    sumfile.close()
    
    
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
    os_obj = OptimalStatistic(PSRs, pta, gwb_name='gw')
    os_obj.set_orf(['hd'])
    
    xi, rhok, sigk, covk, Sk, _, _ = os_obj.compute_PFOS(select_freq = args.fbin, 
                                                         return_pair_vals = True,
                                                         params = noisedict,
                                                         pair_covariance = args.incCrossCorr,
                                                         narrowband=True)
    
    print("... OS calculation done for frequency bin {}".format(args.fbin),
          '\n')
    #==========================================================================
        
     
        
     
        
    #==========================================================================
    print('\n' + 50*'-' + '\n' + 50*'-' + '\n',
          'ANISOTROPY ANALYSIS')
    

    if args.incCrossCorr == True:
        os_cov = covk
    else:
        os_cov = np.diag(covk)
 
    
     
    # creating anis PTA with complex spherical harmonics
    apta = anis_pta.anis_pta(psrs_theta, psrs_phi, xi, rhok, sigk, os_cov = os_cov,
                             nside = NSIDE, l_max = args.lmax,
                             os = Sk, mode = 'sph_basis',
                             include_pta_monopole = args.incMonopole)

    print('... built APTA object', flush=True)
    
 
    
    # crosscheck with S/N ratios
    if args.check_iso is True:
        print('... check with isotropic OS analysis')
        _, _, _, _, os, os_sigma, _ = os_obj.compute_OS(params=noisedict, pair_covariance=False, return_pair_vals = True)
        os_snr = os/os_sigma
        print(f"       Traditional OS A^2 estimate: {os} +/- {os_sigma} with SNR: {os_snr}",
              '\n')
        
        M_raw = apta.fisher_matrix_sph(cross_corr = False)
        _, _, _, _, _, X_raw, _, _ = apta.max_lkl_clm(use_svd_reg=True, cross_corr = False)
 
        P_00 = X_raw[0]/(M_raw[0][0])
        sigma_00 = (M_raw[0][0])**(-1./2.)
        sn_00 = P_00/sigma_00
        
        print("      P_00 = [M_(00)(00)]^-1 X_(00) = ", P_00)
        print("      sigma_00 = [M_(00)(00)]^-1/2 = ", sigma_00)
        print("      (S/N)_00 := P_00 / sigma_00 = ", sn_00)
        print("      (S/N)_OS = ", os_snr)
        if np.abs(sn_00 - os_snr) < 0.019:
            print('      Analyses agree.\n\n', flush=True)
        else:
            print('      Analyses disagree.\n\n', flush=True)

    
    
    if args.cutoff == 0:
        print('... deriving solution for recovered power without regularisation')

        print('      calculate c_lm values')
        clms, clm_err, cn, sv, Minv, X, U, Vh = apta.max_lkl_clm(use_svd_reg=True,
                                                                 cross_corr = args.incCrossCorr)
            
        # refer to Minv as Mprimeinv, for storage reasons
        Mprime_inv = Minv
        cov_M = Minv
        M = apta.fisher_matrix_sph(cross_corr=args.incCrossCorr)

    else:
        print('... deriving solution for recovered power with regularisation')

        print('      calculate c_lm values')
        clms, clms_err, cn, sv, Mprime_inv, X, U, Vh = apta.max_lkl_clm(use_svd_reg=True,
                                                                        cutoff=args.cutoff-1,
                                                                        cross_corr=args.incCrossCorr)
                    

        M = apta.fisher_matrix_sph(cross_corr=args.incCrossCorr)
        sv_prime = np.copy(sv)
        print('      calculate covariance matrix')
        cov_M = Mprime_inv @ M @ Mprime_inv


        if args.create_plots is True:
            print('      plot singular values')
            plt.errorbar(range(len(sv)), sv, ls='', fmt='ko', ms=3)
            plt.xlabel('number of eigenvalue')
            plt.ylabel('eigenvalue')
            plt.yscale('log')
            plt.grid(visible=True, which='major')
            plt.savefig(figuredir + '/singularvalues_{}_lmax{}_nside{}_fbins{}.png'.format(basis, args.lmax, NSIDE, args.fbin),
                        bbox_inches='tight', dpi=400)
            plt.clf()
    #==========================================================================




    #==========================================================================
    print('\n')
    
    if args.outdir is not None:
        if args.incCrossCorr == False :
            save_complex_matrix(args.outdir + '/M_sph/M_sph_nocc_{}_{}_{}_fbin{}.dat'.format(args.lmax, args.cutoff, args.nside, args.fbin), M)
            save_complex_matrix(args.outdir + '/Mprimeinv_sph/Mprimeinv_sph_nocc_{}_{}_{}_fbin{}.dat'.format(args.lmax, args.cutoff, args.nside, args.fbin), Mprime_inv)
            save_complex_matrix(args.outdir + '/X_sph/X_sph_nocc_{}_{}_{}_fbin{}.dat'.format(args.lmax, args.cutoff, args.nside, args.fbin), X)
            save_complex_matrix(args.outdir + '/Cov_sph/Cov_sph_nocc_{}_{}_{}_fbin{}.dat'.format(args.lmax, args.cutoff, args.nside, args.fbin), cov_M)
            save_complex_matrix(args.outdir + '/popt_sph/popt_sph_nocc_{}_{}_{}_fbin{}.dat'.format(args.lmax, args.cutoff, args.nside, args.fbin), clms)
        else:
            save_complex_matrix(args.outdir + '/M_sph/M_sph_cc_{}_{}_{}_fbin{}.dat'.format(args.lmax, args.cutoff, args.nside, args.fbin), M)
            save_complex_matrix(args.outdir + '/Mprimeinv_sph/Mprimeinv_sph_cc_{}_{}_{}_fbin{}.dat'.format(args.lmax, args.cutoff, args.nside, args.fbin), Mprime_inv)
            save_complex_matrix(args.outdir + '/X_sph/X_sph_cc_{}_{}_{}_fbin{}.dat'.format(args.lmax, args.cutoff, args.nside, args.fbin), X)
            save_complex_matrix(args.outdir + '/Cov_sph/Cov_sph_cc_{}_{}_{}_fbin{}.dat'.format(args.lmax, args.cutoff, args.nside, args.fbin), cov_M)
            save_complex_matrix(args.outdir + '/popt_sph/popt_sph_cc_{}_{}_{}_fbin{}.dat'.format(args.lmax, args.cutoff, args.nside, args.fbin), clms)

        print('SAVED RESULTS TO {}'.format(args.outdir))
    else:
        print('I have not saved anything :O ... da dumm...')
    #==========================================================================
    
    
    
    return None





if __name__ == '__main__':
    main()


