#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 12:55:29 2023

@author: kgrunthal
"""


import numpy as np
import glob

#matplotlib.use('TKAgg')
import time, datetime

from enterprise.pulsar import Pulsar
from enterprise_extensions import sampler

from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc



#------------------------------------------------------------------------------
def get_psrname(file,name_sep='.'):
    return file.split('/')[-1].split(name_sep)[0]
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
def load_pulsars(pars_dir, tims_dir, psr_list):
    PSRs = []
    for pulsar_name in psr_list:
        parfile = pars_dir + pulsar_name + '.par'
        timfile = tims_dir + pulsar_name + '.tim'
        
        psr = Pulsar(parfile, timfile, ephem='DE440')
        PSRs.append(psr)
    
    return PSRs
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
def add_customised_jump_proposals(jp, sp, pta):
    # Red noise prior draw
    if 'red noise' in jp.snames:
        print('Adding red noise prior draws...')
        sp.addProposalToCycle(jp.draw_from_red_prior, 5)

    # DM GP noise prior draw
    if 'dm_gp' in jp.snames:
        print('Adding DM GP noise prior draws...')
        sp.addProposalToCycle(jp.draw_from_dm_gp_prior, 5)
        
        # DM annual prior draw
    if 'dm_s1yr' in jp.snames:
        print('Adding DM annual prior draws...')
        sp.addProposalToCycle(jp.draw_from_dm1yr_prior, 2)

    # DM dip prior draw
    if 'dmexp' in jp.snames:
        print('Adding DM exponential dip prior draws...')
        sp.addProposalToCycle(jp.draw_from_dmexpdip_prior, 2)

    # DM cusp prior draw
    if 'dm_cusp' in jp.snames:
        print('Adding DM exponential cusp prior draws...')
        sp.addProposalToCycle(jp.draw_from_dmexpcusp_prior, 2)

    # DMX prior draw
    if 'dmx_signal' in jp.snames:
        print('Adding DMX prior draws...')
        sp.addProposalToCycle(jp.draw_from_dmx_prior, 10)

    # Chromatic GP noise prior draw
    if 'chrom_gp' in jp.snames:
        print('Adding chromatic GP noise prior draws...')
        sp.addProposalToCycle(jp.draw_from_chrom_gp_prior, 2)

    # SW prior draw
    if 'gp_sw' in jp.snames:
        print('Adding Solar Wind DM GP prior draws...')
        sp.addProposalToCycle(jp.draw_from_dm_sw_prior, 2)
        #sp.addProposalToCycle(jp.draw_from_par_prior('gp_sw'), 2)
        
    if 'n_earth' in jp.snames:
        print('Adding deterministic Solar Wind prior draws...')
        sp.addProposalToCycle(jp.draw_from_det_sw_prior, 2)
        #sp.addProposalToCycle(jp.draw_from_par_prior('n_earth'), 2)

    # Chromatic GP noise prior draw
    if 'dm_bump' in jp.snames:
        print('Adding chromatic Bump noise prior draws...')
        sp.addProposalToCycle(jp.draw_from_chrombump_prior, 2)
        #sp.addProposalToCycle(jp.draw_from_par_prior('chrom_bump'), 2)
        
    # Chromatic GP noise prior draw
    if 'chrom1yr' in jp.snames:
        print('Adding annual Chromatic Variation prior draws...')
        sp.addProposalToCycle(jp.draw_from_annualchrom_prior, 2)
        #sp.addProposalToCycle(jp.draw_from_par_prior('chrom1yr'), 2)
    
    # Ephemeris prior draw
    if 'd_jupiter_mass' in pta.param_names:
        print('Adding ephemeris model prior draws...')
        sp.addProposalToCycle(jp.draw_from_ephem_prior, 10)

    # GWB uniform distribution draw
    if np.any([('gw' in par and 'log10_A' in par) for par in pta.param_names]):
        print('Adding GWB uniform distribution draws...')
        sp.addProposalToCycle(jp.draw_from_gwb_log_uniform_distribution, 10)

    # Dipole uniform distribution draw
    if 'dipole_log10_A' in pta.param_names:
        print('Adding dipole uniform distribution draws...')
        sp.addProposalToCycle(jp.draw_from_dipole_log_uniform_distribution, 10)

    # Monopole uniform distribution draw
    if 'monopole_log10_A' in pta.param_names:
        print('Adding monopole uniform distribution draws...')
        sp.addProposalToCycle(jp.draw_from_monopole_log_uniform_distribution, 10)

    # Altpol uniform distribution draw
    if 'log10Apol_tt' in pta.param_names:
        print('Adding alternative GW-polarization uniform distribution draws...')
        sp.addProposalToCycle(jp.draw_from_altpol_log_uniform_distribution, 10)

    # BWM prior draw
    if 'bwm_log10_A' in pta.param_names:
        print('Adding BWM prior draws...')
        sp.addProposalToCycle(jp.draw_from_bwm_prior, 10)

    # FDM prior draw
    if 'fdm_log10_A' in pta.param_names:
        print('Adding FDM prior draws...')
        sp.addProposalToCycle(jp.draw_from_par_prior('n_earth'), 2)

    # CW prior draw
    if 'log10_h' in pta.param_names:
        print('Adding CGW strain prior draws...')
        sp.addProposalToCycle(jp.draw_from_cw_log_uniform_distribution, 6)

    if 'log10_Mc' in pta.param_names:
        print('Adding CGW prior draws...')
        sp.addProposalToCycle(jp.draw_from_cw_prior, 6)
        #sp.addProposalToCycle(jp.draw_from_par_log_uniform('cw'), 5)

    # free spectrum prior draw
    if np.any(['log10_rho' in par for par in pta.param_names]):
        print('Adding free spectrum prior draws...')
        sp.addProposalToCycle(jp.draw_from_gw_rho_prior, 25)
#------------------------------------------------------------------------------
       
   
   
   

#------------------------------------------------------------------------------
def setup_sampler_single(pta, outdir='MCMC_out', loglkwargs={}, logpkwargs={}):
       
    ### SETUP SAMPLER ###
    # start from prior draw, update from empirical distr.
    p0 = {p.name: p.sample() for p in pta.params}
    x0 = np.hstack(list(p0.values()))
    ndim = len(x0)

    # set initial cov stdev to (starting order of magnitude)/10
    stdev = np.array([10**np.floor(np.log10(abs(val)))/10 for val in x0])
    cov = np.diag(stdev**2)
    
    # get groups

    groups = sampler.get_parameter_groups(pta)
    groups.extend(sampler.get_psr_groups(pta))

    sp = ptmcmc(ndim, pta.get_lnlikelihood, pta.get_lnprior, cov,
                     groups=groups, outDir=outdir, resume=False,
                     loglkwargs=loglkwargs, logpkwargs=logpkwargs)

    sampler.save_runtime_info(pta, sp.outDir, human=None)

    # additional jump proposals
    jp = sampler.JumpProposal(pta, outdir=outdir)
    sp.jp = jp

    # always add draw from prior
    sp.addProposalToCycle(jp.draw_from_prior, 5)
    add_customised_jump_proposals(jp, sp, pta)
    '''
    sel_sig = ["rn", "red_noise", "dm_gp", "fcn", "chrom-rn", "srn", "dm_srn",
               "freechrom-srn", "chrom-srn", "dm-expd", "freechrom-expd",
               "chrom-expd", "dm-y", "freechrom-y", "chrom-y"]
    for s in sel_sig:
        if any([s in p for p in pta.param_names]):
            #pnames = [p.name for p in pta.params if s in p.name]
            #print('Adding %s prior draws with parameters :'%s, pnames, '\n')
            print('Adding %s prior draws.'%s)
            sp.addProposalToCycle(jp.draw_from_par_prior(s), 10)
    '''

    # try adding empirical proposals
    #if empirical_distr is not None:
    #    print('Attempting to add empirical proposals...\n')
    #    sp.addProposalToCycle(jp.draw_from_empirical_distr, 10)

    return sp, x0
#------------------------------------------------------------------------------





#------------------------------------------------------------------------------
def setup_sampler_hyper(hmodel, outdir='MCMC_out',
                        empirical_distr=None, groups=None,
                        loglkwargs={}, logpkwargs={}):
    
    x0 = hmodel.initial_sample()
    while hmodel.get_lnlikelihood(x0) < 0. or hmodel.get_lnprior(x0) == float(-np.inf):
        x0 = hmodel.initial_sample()
        print('... looking for a suitable initial sample', end='\r')

    ndim = len(x0)

    # set initial cov stdev to (starting order of magnitude)/10
    stdev = np.array([10**np.floor(np.log10(abs(val)))/10 for val in x0])
    cov = np.diag(stdev**2)
    
    # get groups
    groups = hmodel.get_parameter_groups()
    
    sp = ptmcmc(ndim, hmodel.get_lnlikelihood, hmodel.get_lnprior, cov,
                     groups=groups, outDir=outdir, resume=False,
                     loglkwargs=loglkwargs, logpkwargs=logpkwargs)

    sampler.save_runtime_info(hmodel, sp.outDir, human=None)

    # additional jump proposals
    jp = sampler.JumpProposal(hmodel, hmodel.snames, outdir=outdir,
                              empirical_distr=empirical_distr)
    sp.jp = jp

    # always add draw from prior
    sp.addProposalToCycle(jp.draw_from_prior, 5)
    
    add_customised_jump_proposals(jp, sp, hmodel)
    
    '''
    sel_sig = ["rn", "red_noise", "dm_gp", "fcn", "chrom-rn", "srn", "dm_srn",
               "freechrom-srn", "chrom-srn", "dm-expd", "freechrom-expd",
               "chrom-expd", "dm-y", "freechrom-y", "chrom-y"]
    for s in sel_sig:
        if any([s in p for p in hpta.param_names]):
            #pnames = [p.name for p in pta.params if s in p.name]
            #print('Adding %s prior draws with parameters :'%s, pnames, '\n')
            print('Adding %s prior draws.'%s)
            sp.addProposalToCycle(jp.draw_from_par_prior(s), 10)
    '''
    if 'nmodel' in hmodel.param_names:
        print('Adding nmodel uniform distribution draws...\n')
        sp.addProposalToCycle(hmodel.draw_from_nmodel_prior, 25)


    # try adding empirical proposals
    if empirical_distr is not None:
        print('Attempting to add empirical proposals...\n')
        sp.addProposalToCycle(jp.draw_from_empirical_distr, 10)

    return sp, x0

    
