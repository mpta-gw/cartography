# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 15:02:36 2026

@author: kgrun
"""


import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm

from .maps_mod import anis_pta
from .map_conversion_utils import makemap, diag_pixel, get_sigma_map
from .map_conversion_utils import build_pixel_conversion_fast as conv_fast





def MLsolution_with_reg(PSRs, os_obj, noisedict, setup_dict, outdir,
                        narrowband = True, incCrossCorr = True, check_iso = False,
                        plot_sv = False):
    
    
    xi, rhok, sigk, covk, Sk, _, _ = os_obj.compute_PFOS(select_freq = setup_dict["fbin"], 
                                                         return_pair_vals = True,
                                                         params = noisedict,
                                                         pair_covariance = incCrossCorr,
                                                         narrowband=narrowband)
    
    print("... PFOS calculation done for frequency bin {}".format(setup_dict["fbin"]),
          '\n')

    if incCrossCorr == True:
        os_cov = covk
        cc_tag = 'cc'
    else:
        os_cov = np.diag(covk)
        cc_tag = 'nocc'
 
    
     
    # creating anis PTA with complex spherical harmonics
    # get positions of pulsars
    psrs_theta, psrs_phi = np.zeros(len(PSRs)), np.zeros(len(PSRs))
    for pp, psr in enumerate(PSRs):
        psrs_theta[pp] = psr.theta
        psrs_phi[pp] = psr.phi

    print('... builing APTA object', flush=True)
        
    apta = anis_pta.anis_pta(psrs_theta, psrs_phi, xi, rhok, sigk, os_cov = os_cov,
                             nside = setup_dict["NSIDE"], l_max = setup_dict["lmax"],
                             os = Sk, mode = 'sph_basis_fast',
                             include_pta_monopole = True)

    print('\t done', flush=True)
    print()
 
    
    # crosscheck with S/N ratios
    if check_iso is True:
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

    


    #==========================================================================
    if setup_dict["cutoff"] == 0:
        print('... deriving solution for recovered power without regularisation')

        print('      calculate c_lm values')
        clms, clm_err, cn, sv, Minv, X, U, Vh = apta.max_lkl_clm(use_svd_reg=True,
                                                                 cross_corr = incCrossCorr)
            
        # refer to Minv as Mprimeinv, for storage reasons
        Mprime_inv = Minv
        cov_M = Minv
        M = apta.fisher_matrix_sph(cross_corr=incCrossCorr)



    else:
        print('... deriving solution for recovered power with regularisation')

        print('      calculate c_lm values')
        clms, clms_err, cn, sv, Mprime_inv, X, _, _ = apta.max_lkl_clm(use_svd_reg = True,
                                                                        cutoff = setup_dict["cutoff"]-1,
                                                                        cross_corr = incCrossCorr)
                    

        M = apta.fisher_matrix_sph(cross_corr = incCrossCorr)
        print('      calculate covariance matrix')
        cov_M = Mprime_inv @ M @ Mprime_inv


        if plot_sv is True:
            print('      plot singular values')
            plt.errorbar(range(len(sv)), sv, ls='', fmt='ko', ms=3)
            plt.xlabel('number of eigenvalue')
            plt.ylabel('eigenvalue')
            plt.yscale('log')
            plt.grid(visible=True, which='major')
            plt.savefig(outdir + '/figures/singularvalues_lmax{}_nside{}_fbins{}.png'.format(setup_dict["lmax"], setup_dict["cutoff"], setup_dict["fbin"]),
                        bbox_inches='tight', dpi=400)
            plt.clf()
    #==========================================================================


    return clms, clms_err, X, M, cov_M, Mprime_inv


    
    
    
    

def convert_maps(clms, M, X, cov_M, Mprime_inv, setup_dict, outdir,
                 calculate_pvalue = False):
    print('\n' + 50*'-' + '\n')
    print('MAP CONVERSION')
    
    
    res = setup_dict["res"]    # pixel resolution in degrees
    
    
    # ---- 1. radiometer map --------------------------------------------------
    print('... computing radiometer map')
    diag_fisher, RA, DEC, _, U = diag_pixel(M, res)
    X_pix         = U @ X
    p_radio       = np.real(X_pix / diag_fisher)
    s_radio       = np.real(diag_fisher ** -0.5)
    sn_radiometer = p_radio / s_radio
    
    
    # ---- 2. clean map -------------------------------------------------------
    print('... computing clean map')
    clean_map, _, _, _ = makemap(clms, res)
    
    clms_nmp = np.copy(clms)
    clms_nmp[0] = 0
    clean_map_nmp, _, _, _ = makemap(clms_nmp, res)
    
    # ---- 3. S/N map ---------------------------------------------------------
    print('... computing S/N map')
    sigmaPix, _, _, _, _ = get_sigma_map(cov_M, res)
    sn = clean_map / sigmaPix
    
    
    # ---- 4. sensitivity map -------------------------------------------------
    print('... computing sensitivity map')
    Lambda = Mprime_inv @ M
    diag_LambdaPix = np.einsum('ij,jk,ik->i', U, Lambda, U.conj())
    sn_sensitivity = np.real(diag_LambdaPix) / sigmaPix
    
    
    # =========================================================================    
    out_maps = np.column_stack([RA, DEC, clean_map, sigmaPix,
                                sn, sn_radiometer, sn_sensitivity, clean_map_nmp])
    np.savetxt(outdir + '/maps_{}_{}_{}_fbin{}.txt'.format(setup_dict["lmax"], setup_dict["cutoff"], setup_dict["NSIDE"], setup_dict["fbin"]),
               out_maps, delimiter='\t')
    
    out_fisher = np.column_stack([RA, DEC, np.real(diag_fisher)])
    np.savetxt(outdir + '/map_M_{}_{}_{}_fbin{}.txt'.format(setup_dict["lmax"], setup_dict["cutoff"], setup_dict["NSIDE"], setup_dict["fbin"]),
               out_fisher, delimiter='\t')
    
    print('SAVED MAPS TO {}\n'.format(outdir))
    # =========================================================================
    
    
    
    
    
    # ---- 5. Monte Carlo significance ----------------------------------------
    if calculate_pvalue is True:
        # 
        print('... computing Monte Carlo significance')
        # draw random dirty maps
        # Fisher matrix is the covariance matrix of the dirty map, X
        n_trials = 10000
        rng      = np.random.default_rng()
        
        ev_M, V   = np.linalg.eigh(M)
        sqrt_ev   = np.sqrt(np.real(ev_M))
    
        distr_snr = np.empty((n_trials, 6))
        
        UVsqrtev = U @ V * sqrt_ev[None, :]
        print(np.shape(U @ V), np.shape(UVsqrtev))
        MinvVsqrtev = Mprime_inv @ V * sqrt_ev[None, :]
    
        for kk in tqdm(range(n_trials)):
            rv                  = rng.standard_normal(len(clms))
            X_pix_rand          = UVsqrtev @ rv
            p_rand              = np.real(X_pix_rand / diag_fisher)
            radio_snr           = p_rand / s_radio
            distr_snr[kk, 0:2]  = np.array([radio_snr.max(), radio_snr.min()])
            
            clm_rand = MinvVsqrtev * rv
            clean_map_rand, _, _, _ = makemap(clm_rand, res)
            clean_snr               = clean_map_rand / sigmaPix
            distr_snr[kk, 2:4]     = np.array([clean_snr.max(), clean_snr.min()])
            
            clm_rand[0] = 0
            clean_map_rand_nmp, _, _, _ = makemap(clm_rand, res)
            clean_nmp_snr               = clean_map_rand_nmp / sigmaPix
            distr_snr[kk, 4:]           = np.array([clean_nmp_snr.max, clean_nmp_snr.min])
            

        np.savetxt(outdir + '/distribution_maxmin_{}_{}_{}_fbin{}.txt'.format(
                   setup_dict["lmax"], setup_dict["cutoff"], setup_dict["NSIDE"], setup_dict["fbin"]),
                   distr_snr, delimiter='\t')
        
        return out_maps, out_fisher, distr_snr
    
    
    else:
        return out_maps, out_fisher
    
    
    

def convert_maps_fast(clms, M, X, cov_M, Mprime_inv, setup_dict, outdir,
                      calculate_pvalue = False):
    print('\n' + 50*'-' + '\n')
    print('MAP CONVERSION ( FAST :D )')
    
    
    res = setup_dict["res"]    # pixel resolution in degrees
    
    pc = conv_fast(setup_dict["lmax"], res)
    U  = pc['U']
    
    # ---- 1. radiometer map with S/N -----------------------------------------
    print('... computing radiometer map')
    diag_fisher   = np.einsum('ij,jk,ik->i', U, M, np.conj(U) )
    X_pix         = U @ X
    p_radio       = np.real(X_pix / diag_fisher)
    s_radio       = np.real(diag_fisher ** -0.5)
    sn_radiometer = p_radio / s_radio
    
    
    # ---- 2. clean map with S/N ----------------------------------------------
    print('... computing clean map')
    clean_map     = np.real( U @ clms)
    sigmaPix      = np.sqrt( np.real( np.einsum('ij,jk,ik->i', U, cov_M, np.conj(U)) ) )
    sn            = clean_map / sigmaPix
    
    clms_nmp      = np.copy(clms)
    clms_nmp[0]   = 0
    clean_map_nmp = np.real( U @ clms_nmp)
    
     
    # ---- 4. sensitivity map -------------------------------------------------
    print('... computing sensitivity map')
    Lambda = Mprime_inv @ M
    diag_LambdaPix = np.einsum('ij,jk,ik->i', U, Lambda, U.conj())
    sn_sensitivity = np.real(diag_LambdaPix) / sigmaPix
    
    
    # =========================================================================    
    out_maps = np.column_stack([pc["RA"], pc["DEC"], clean_map, sigmaPix,
                                sn, sn_radiometer, sn_sensitivity, clean_map_nmp])
    np.savetxt(outdir + '/maps_{}_{}_{}_fbin{}.txt'.format(setup_dict["lmax"], setup_dict["cutoff"], setup_dict["NSIDE"], setup_dict["fbin"]),
               out_maps, delimiter='\t')
    
    out_fisher = np.column_stack([pc["RA"], pc["DEC"], np.real(diag_fisher)])
    np.savetxt(outdir + '/map_M_{}_{}_{}_fbin{}.txt'.format(setup_dict["lmax"], setup_dict["cutoff"], setup_dict["NSIDE"], setup_dict["fbin"]),
               out_fisher, delimiter='\t')
    
    print('SAVED MAPS TO {}\n'.format(outdir))
    # =========================================================================
    
    
    
    
    
    # ---- 5. Monte Carlo significance ----------------------------------------
    if calculate_pvalue is True:
        # 
        print('... computing Monte Carlo significance')
        # draw random dirty maps
        # Fisher matrix is the covariance matrix of the dirty map, X
        n_trials = 1000
        rng      = np.random.default_rng()
        
        ev_M, V   = np.linalg.eigh(M)
        sqrt_ev   = np.sqrt(np.real(ev_M))
    
        distr_snr = np.empty((n_trials, 6))
        
        UVsqrtev = U @ V * sqrt_ev[None, :]
        print(np.shape(U @ V), np.shape(UVsqrtev))
        MinvVsqrtev = Mprime_inv @ V * sqrt_ev[None, :]
    
        for kk in tqdm(range(n_trials)):
            rv                 = rng.standard_normal(len(clms))
            X_pix_rand         = UVsqrtev @ rv
            p_rand             = np.real(X_pix_rand / diag_fisher)
            radio_snr          = p_rand / s_radio
            distr_snr[kk, 0:2] = np.array([radio_snr.max(), radio_snr.min()])
            
            clm_rand           = MinvVsqrtev @ rv
            clean_map_rand     = np.real( U @ clm_rand )
            clean_snr          = clean_map_rand / sigmaPix
            distr_snr[kk, 2:4] = np.array([clean_snr.max(), clean_snr.min()])
            
            clm_rand[0] = 0
            clean_map_rand_nmp = np.real( U @ clm_rand )
            clean_nmp_snr      = clean_map_rand_nmp / sigmaPix
            distr_snr[kk, 4:]  = np.array([clean_nmp_snr.max(), clean_nmp_snr.min()])
            

        np.savetxt(outdir + '/distribution_maxmin_{}_{}_{}_fbin{}.txt'.format(
                   setup_dict["lmax"], setup_dict["cutoff"], setup_dict["NSIDE"], setup_dict["fbin"]),
                   distr_snr, delimiter='\t')
        
        return out_maps, out_fisher, distr_snr
    
    
    else:
        return out_maps, out_fisher





