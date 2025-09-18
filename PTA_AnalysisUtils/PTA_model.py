#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 15:50:20 2024

@author: kgrunthal
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 14:34:32 2024

@author: kgrunthal
"""




import numpy as np

import json

import enterprise.signals.parameter as eparameter
from enterprise.signals import utils
from enterprise.signals import signal_base
from enterprise.signals import selections
from enterprise.signals import white_signals
from enterprise.signals import gp_signals
from enterprise.signals import deterministic_signals
from enterprise.signals import gp_bases
import enterprise.constants as const
from enterprise_extensions import deterministic
from enterprise_extensions.chromatic.solar_wind import solar_wind, createfourierdesignmatrix_solar_dm

from enterprise_extensions import model_utils, blocks


def get_cgw_parameter(key, cgw_dict):
    if len(cgw_dict[key]) == 1:
        parameter = eparameter.Constant(cgw_dict[key])(key)
        return parameter
    elif len(cgw_dict[key]) == 2:
        parameter = eparameter.Uniform(cgw_dict[key][0], cgw_dict[key][1])(key)
        return parameter
    else:
        print(key, 'is set wrong')
        return None

def get_cgw_parameter_pi(key, cgw_dict):
    if len(cgw_dict[key]) == 1:
        parameter = eparameter.Constant(cgw_dict[key]*np.pi)(key)
        return parameter
    elif len(cgw_dict[key]) == 2:
        parameter = eparameter.Uniform(cgw_dict[key][0]*np.pi, cgw_dict[key][1]*np.pi)(key)
        return parameter
    else:
        print(key, 'is set wrong')
        return None



def dm_noise(log10_A,gamma,Tspan,components=30,option="powerlaw"):
    """
    A term to account for stochastic variations in DM. It is based on spin
    noise model, with Fourier amplitudes depending on radio frequency nu
    as ~ 1/nu^2.
    """
    nfreqs = 30
    if option=="powerlaw":
      #pl = utils.powerlaw(log10_A=log10_A, gamma=gamma, components=components)
      pl = utils.powerlaw(log10_A=log10_A, gamma=gamma)
      #pl = enterprise.signals.gp_priors.powerlaw_no_components(log10_A=log10_A, gamma=gamma)

    #elif option=="turnover":
    #  fc = parameter.Uniform(self.params.sn_fc[0],self.params.sn_fc[1])
    #  pl = powerlaw_bpl(log10_A=log10_A, gamma=gamma, fc=fc,
    #                    components=components)
    dm_basis = utils.createfourierdesignmatrix_dm(nmodes = components,
                                                  Tspan=Tspan)
    dmn = gp_signals.BasisGP(pl, dm_basis, name='dm_gp')

    return dmn

@signal_base.function
def chrom_yearly_sinusoid(toas, freqs, log10_Amp, phase, idx):
    """
    Chromatic annual sinusoid.
    :param log10_Amp: amplitude of sinusoid
    :param phase: initial phase of sinusoid
    :param idx: index of chromatic dependence
    :return wf: delay time-series [s]
    """

    wf = 10**log10_Amp * np.sin(2 * np.pi * const.fyr * toas + phase)
    return wf * (1400 / freqs) ** idx


@signal_base.function
def chrom_gaussian_bump(toas, freqs, log10_Amp=-2.5, sign_param=1.0,
                    t0=53890, sigma=81, idx=2):
    """
    Chromatic time-domain Gaussian delay term in TOAs.
    Example: J1603-7202 in Lentati et al, MNRAS 458, 2016.
    """
    #t0 *= const.day
    #sigma *= const.day
    wf = 10**log10_Amp * np.exp(-(toas - t0)**2/2/sigma**2)
    return np.sign(sign_param) * wf * (1400 / freqs) ** idx

def dm_gaussian_bump(tmin, tmax, idx=2, sigma_min=600000, sigma_max=1000,
    log10_A_low=-10, log10_A_high=-1, name='dm_bump'):
    """
    Returns chromatic Gaussian bump (i.e. TOA advance):
    :param tmin, tmax:
        search window for exponential cusp time.
    :param idx:
        index of radio frequency dependence (i.e. DM is 2). If this is set
        to 'vary' then the index will vary from 1 - 6
    :param sigma_min, sigma_max:
        standard deviation of a Gaussian in MJD
    :param sign:
        [boolean] allow for positive or negative exponential features.
    :param name: Name of signal
    :return dm_bump:
        chromatic Gaussian bump waveform.
    """
    sign_param = eparameter.Uniform(-1,1)
    t0_dm_bump = eparameter.Uniform(tmin,tmax)
    sigma_dm_bump = eparameter.Uniform(sigma_min,tmax-tmin)
    log10_Amp_dm_bump = eparameter.Uniform(log10_A_low, log10_A_high)
    if idx == 'vary':
        idx = eparameter.Uniform(0, 14)
    wf = chrom_gaussian_bump(log10_Amp=log10_Amp_dm_bump,
                         t0=t0_dm_bump, sigma=sigma_dm_bump,
                         sign_param=sign_param, idx=idx)
    dm_bump = deterministic_signals.Deterministic(wf, name=name)

    return dm_bump



def mk_band_split(freqs, backend_flags):
    """ Selection for splitting the band in 3 """
    d =  dict(zip(['low', 'high', 'all'],
                  [(freqs <= 1100), (freqs > 1100), (freqs >= 0)]))
    delkeys = []
    for key in d.keys():
        if np.sum(d[key]) == 0:  # all False
            delkeys.append(key)
    for key in delkeys:
        del d[key]
    return d








###############################################################################
# actual MPTA model
###############################################################################

def PTA(PSRs, ml_noiseparam_dict, SPNA_models, cs,
         orfmatrix=None, psrspos = None, verbose=True,
         **ptaparams):
    
    no_selection=selections.Selection(selections.no_selection)
    by_backend = selections.Selection(selections.by_backend)
    
    

    Tspan=model_utils.get_tspan(PSRs)

    efac = eparameter.Constant()
    equad = eparameter.Constant()
    ecorr = eparameter.Constant()

    models = []
    s_common = []
    
    #-- common signals --------------------------------------------------------
    string_common = 'PTA common signal:'
    for key in cs.keys():
        if key == 'CRN':
        
            if 'fbin' in ptaparams.keys():
                f_bin = ptaparams['fbin']
                if f_bin==1:
                    CRN = blocks.common_red_noise_block(psd=cs['CRN']['psd'], prior='log-uniform',
                                                        Tspan=Tspan,
                                                        components=1,gamma_val=3.41,
                                                        orf=cs['CRN']['orf'],
                                                        orf_matrix = orfmatrix,
                                                        psrs_pos = psrspos,
                                                        name='gw')
                    s_common.append(CRN)
                    string_common += ' CRN '

                else:
                    DUMMY = blocks.common_red_noise_block(psd=cs['CRN']['psd'], prior='log-uniform',
                                                          Tspan=Tspan,
                                                          components=(f_bin-1),gamma_val=3.41,
                                                          orf=cs['CRN']['orf'],
                                                          orf_matrix = orfmatrix,
                                                          psrs_pos = psrspos,
                                                          name='dummy')
            
                    CRN = blocks.common_red_noise_block(psd=cs['CRN']['psd'], prior='log-uniform',
                                                        Tspan=Tspan/f_bin,
                                                        components=1,gamma_val=3.41,
                                                        orf=cs['CRN']['orf'],
                                                        orf_matrix = orfmatrix,
                                                        psrs_pos = psrspos,
                                                        name='gw')
                    s_common.append(DUMMY)
                    s_common.append(CRN)
                    string_common += ' CRN+DUMMY '
                
            
            else:
                CRN = blocks.common_red_noise_block(psd=cs['CRN']['psd'], prior='log-uniform',
                                                Tspan=Tspan,
                                                components=cs['CRN']['comp'],
                                                gamma_val=cs['CRN']['gamma'],
                                                orf=cs['CRN']['orf'],
                                                orf_matrix = orfmatrix,
                                                psrs_pos = psrspos,
                                                name='gw')
                s_common.append(CRN)
                string_common += ' CRN '
        
        
        
        
        elif key == 'CGW':
            cgw_dict = cs['CGW']
            
            cos_gwtheta = get_cgw_parameter('cos_gwtheta', cgw_dict)
            gwphi = get_cgw_parameter_pi('gwphi', cgw_dict)		# position of source
            log10_mc = get_cgw_parameter('log10_Mc', cgw_dict)		# chirp mass
            log10_h = get_cgw_parameter('log10_h', cgw_dict)		# strain amplitude
            log10_fgw = get_cgw_parameter('log10_fgw', cgw_dict)	# gw frequency
            phase0 = get_cgw_parameter_pi('phase0', cgw_dict)		# gw phase
            psi = get_cgw_parameter_pi('psi', cgw_dict)			# gw polarization 
            cos_inc = get_cgw_parameter('cos_inc', cgw_dict)		# inclination of binary with respect to Earth
            '''
            cos_gwtheta = eparameter.Uniform(-1, 1)                  # position of source
            gwphi = eparameter.Uniform(0, 2*np.pi)                   # position of source
            log10_mc = eparameter.Uniform(6.5, 10.)                  # chirp mass
            log10_h = eparameter.Uniform(-16, -11)                   # strain amplitude
            log10_fgw = eparameter.Uniform(-8., -7.)                 # gw frequency
            phase0 = eparameter.Uniform(0, 2*np.pi)                  # gw phase
            psi = eparameter.Uniform(0, np.pi)                       # gw polarization 
            cos_inc = eparameter.Uniform(-1, 1)                      # inclination of binary with respect to Earth 
            '''
            cw_wf = deterministic.cw_delay(cos_gwtheta=cos_gwtheta, gwphi=gwphi, log10_mc=log10_mc, 
                                           log10_h=log10_h, log10_fgw=log10_fgw, phase0=phase0, 
                                           psi=psi, cos_inc=cos_inc)
            CGW = deterministic.CWSignal(cw_wf, ecc=False, psrTerm=False)
            s_common.append(CGW)
            string_common += ' CGW '
        
        
        elif key == 'EPH':
            EPH = deterministic_signals.PhysicalEphemerisSignal()
            s_common.append(EPH)
            string_common += ' EPH '
                    



    for psr in PSRs:
        outstring = psr.name + ': '

        s = gp_signals.TimingModel(use_svd=True)

        if psr.name+"_log10_tnequad" in ml_noiseparam_dict.keys():
            eq = white_signals.TNEquadNoise(log10_tnequad=equad, selection=no_selection)
            s += eq
            outstring += 'EQUAD, '
       
        elif psr.name+"_KAT_MKBF_log10_tnequad" in ml_noiseparam_dict.keys():
            eq = white_signals.TNEquadNoise(log10_tnequad=equad, selection=by_backend)
            s += eq
            outstring += 'EQUAD, '
 
   
        if psr.name+"_efac" in ml_noiseparam_dict.keys():
            ef = white_signals.MeasurementNoise(efac=efac, selection=no_selection)
            s += ef
            outstring += 'EFAC, '

        elif psr.name+"_KAT_MKBF_efac" in ml_noiseparam_dict.keys():
            ef = white_signals.MeasurementNoise(efac=efac, selection=by_backend)
            s += ef
            outstring += 'EFAC, '

        
        if any('low' in key for key in ml_noiseparam_dict.keys()):
            mk_ecorr_selection = selections.Selection(mk_band_split)
            if psr.name+"_basis_ecorr_all_log10_ecorr" in ml_noiseparam_dict.keys():
                ec = gp_signals.EcorrBasisModel(log10_ecorr=ecorr,selection=mk_ecorr_selection)
                s += ec
                outstring += 'ECORR(gp), '
            elif psr.name+"_all_log10_ecorr" in ml_noiseparam_dict.keys():
                ec = white_signals.EcorrKernelNoise(log10_ecorr=ecorr, selection=mk_ecorr_selection)
                s += ec
                outstring += 'ECORR(kernel), '  
            
        else:
            if psr.name+"_basis_ecorr_KAT_MKBF_log10_ecorr" in ml_noiseparam_dict.keys():
                ec = gp_signals.EcorrBasisModel(log10_ecorr=ecorr,selection=by_backend)
                s += ec
                outstring += 'ECORR(gp), '
            elif psr.name+"_KAT_MKBF_log10_ecorr" in ml_noiseparam_dict.keys():
                ec = white_signals.EcorrKernelNoise(log10_ecorr=ecorr, selection=by_backend)
                s += ec
                outstring += 'ECORR(kernel), '  
            
        tmin = psr.toas.min()
        tmax = psr.toas.max()
        Tspan=(tmax-tmin)

        max_cadence = 60  # days
        components = 120
        high_comps = 120
        
        # Get list of models
        #psrmodels = [ psr_model for psr_model in SPNA_models if psr.name in psr_model ][0].split("_")[1:]
        
        psrmodels = SPNA_models[psr.name].split(",")

        
        for i, pmodel in enumerate(psrmodels):
            comp = pmodel.split('_')

            pm = comp[0]
            if pm =="DM":
                if ( i+1 < len(psrmodels) and psrmodels[i+1] != "WIDE" ):
                    log10_A_dm = eparameter.Uniform(-20, -11)
                    gamma_dm = eparameter.Uniform(0, 7)
                    dm = dm_noise(log10_A = log10_A_dm, gamma = gamma_dm,
                                  Tspan = Tspan,
                                  components = int(comp[1]),option="powerlaw")
                    s += dm
                    outstring += 'DM, '
                    
                elif i+1 == len(psrmodels):
                    log10_A_dm = eparameter.Uniform(-20, -11)
                    gamma_dm = eparameter.Uniform(0, 7)
                    dm = dm_noise(log10_A=log10_A_dm,gamma=gamma_dm,Tspan=Tspan,components=int(comp[1]),option="powerlaw")
                    s += dm 
                    outstring += 'DM, '
            
            if pm == "DMWIDE" or ( pm == "DM" and ( i+1 < len(psrmodels) and psrmodels[i+1] == "WIDE" ) ):
                log10_A_dm = eparameter.Uniform(-20, -11)
                gamma_dm = eparameter.Uniform(0, 14)
                dm = dm_noise(log10_A=log10_A_dm,gamma=gamma_dm,Tspan=Tspan,components=120,option="powerlaw")
                s += dm
                outstring += 'DM, '

            if pm == "RN":
                if int(comp[1]) == 30:
                    if 'RNsampling' in ptaparams.keys() and ptaparams['RNsampling'] == True:
                        log10_A_red = eparameter.Uniform(-21., -11.)
                        gamma_red = eparameter.Uniform(0., 7.)
                    elif 'RNsampling' in ptaparams.keys() and ptaparams['RNsampling'] == False:
                        log10_A_red = eparameter.Constant()
                        gamma_red = eparameter.Constant()
                    else:
                        log10_A_red = eparameter.Constant()
                        gamma_red = eparameter.Constant()

                
                else:
                    log10_A_red = eparameter.Uniform(-20, -11)
                    gamma_red = eparameter.Uniform(0, 7)
                    
                pl = utils.powerlaw(log10_A=log10_A_red, gamma=gamma_red)
                rn = gp_signals.FourierBasisGP(spectrum=pl, components=int(comp[1]), Tspan=Tspan)
                s += rn
                outstring += 'RN, '

            if pm == "CHROM":
                if ( i+1 < len(psrmodels) and psrmodels[i+1] != "WIDE" ):
                    log10_A_chrom_prior = eparameter.Uniform(-20, -11)
                    gamma_chrom_prior = eparameter.Uniform(0, 14)
                    chrom_gp_idx = eparameter.Uniform(0,14)
                    chrom_model = utils.powerlaw(log10_A=log10_A_chrom_prior, gamma=gamma_chrom_prior)
                    idx = chrom_gp_idx
                    components = high_comps
                    chrom_basis = gp_bases.createfourierdesignmatrix_chromatic(nmodes=components,
                                                                            idx=idx)
                    chrom = gp_signals.BasisGP(chrom_model, chrom_basis, name='chrom_gp')
                    s += chrom
                    outstring += 'CHROM, '
                    
                elif i+1 == len(psrmodels):
                    log10_A_chrom_prior = eparameter.Uniform(-20, -11)
                    gamma_chrom_prior = eparameter.Uniform(0, 14)
                    chrom_gp_idx = eparameter.Uniform(0,14)
                    chrom_model = utils.powerlaw(log10_A=log10_A_chrom_prior, gamma=gamma_chrom_prior)
                    idx = chrom_gp_idx
                    components = high_comps
                    chrom_basis = gp_bases.createfourierdesignmatrix_chromatic(nmodes=components,
                                                                            idx=idx)
                    chrom = gp_signals.BasisGP(chrom_model, chrom_basis, name='chrom_gp')
                    s += chrom
                    outstring += 'CHROM, '

            if pm == "CHROMWIDE" or ( pm == "CHROM" and ( i+1 < len(psrmodels) and psrmodels[i+1] == "WIDE" ) ):
                log10_A_chrom_prior = eparameter.Uniform(-20, -11)
                gamma_chrom_prior = eparameter.Uniform(0, 14)
                chrom_gp_idx = eparameter.Uniform(0,14)
                chrom_model = utils.powerlaw(log10_A=log10_A_chrom_prior, gamma=gamma_chrom_prior)
                idx = chrom_gp_idx
                components = 120
                chrom_basis = gp_bases.createfourierdesignmatrix_chromatic(nmodes=components,
                                                                        idx=idx)
                chrom = gp_signals.BasisGP(chrom_model, chrom_basis, name='chrom_gp')
                s += chrom
                outstring += 'CHROM, '

            if pm == "CHROMCIDX":
                if ( i+1 < len(psrmodels) and psrmodels[i+1] != "WIDE" ):
                    log10_A_chrom_prior = eparameter.Uniform(-20, -11)
                    gamma_chrom_prior = eparameter.Uniform(0, 7)
                    chrom_gp_idx = eparameter.Constant(4)
                    chrom_model = utils.powerlaw(log10_A=log10_A_chrom_prior, gamma=gamma_chrom_prior)
                    idx = chrom_gp_idx
                    components = high_comps
                    chrom_basis = gp_bases.createfourierdesignmatrix_chromatic(nmodes=components,
                                                                            idx=idx)
                    chrom = gp_signals.BasisGP(chrom_model, chrom_basis, name='chrom_gp')
                    s += chrom
                    outstring += 'CHROMCIDX, '
                    
                elif i+1 == len(psrmodels):
                    log10_A_chrom_prior = eparameter.Uniform(-20, -11)
                    gamma_chrom_prior = eparameter.Uniform(0, 7)
                    chrom_gp_idx = eparameter.Constant(4)
                    chrom_model = utils.powerlaw(log10_A=log10_A_chrom_prior, gamma=gamma_chrom_prior)
                    idx = chrom_gp_idx
                    components = high_comps
                    chrom_basis = gp_bases.createfourierdesignmatrix_chromatic(nmodes=components,
                                                                            idx=idx)
                    chrom = gp_signals.BasisGP(chrom_model, chrom_basis, name='chrom_gp')
                    s += chrom
                    outstring += 'CHROMCIDX, '

            if pm == "CHROMCIDXWIDE" or ( pm == "CHROMCIDX" and ( i+1 < len(psrmodels) and psrmodels[i+1] == "WIDE" ) ):
                log10_A_chrom_prior = eparameter.Uniform(-20, -11)
                gamma_chrom_prior = eparameter.Uniform(0, 14)
                chrom_gp_idx = eparameter.Constant(4)
                chrom_model = utils.powerlaw(log10_A=log10_A_chrom_prior, gamma=gamma_chrom_prior)
                idx = chrom_gp_idx
                components = 120
                chrom_basis = gp_bases.createfourierdesignmatrix_chromatic(nmodes=components,
                                                                        idx=idx)
                chrom = gp_signals.BasisGP(chrom_model, chrom_basis, name='chrom_gp')
                s += chrom
                outstring += 'CHROMCIDX, '

            if pm == "CHROMANNUAL":
                log10_Amp_chrom1yr = eparameter.Uniform(-20, -5)
                phase_chrom1yr = eparameter.Uniform(0, 2*np.pi)
                idx_chrom1yr = eparameter.Uniform(0, 14)
                wf = chrom_yearly_sinusoid(log10_Amp=log10_Amp_chrom1yr, phase=phase_chrom1yr, idx=idx_chrom1yr)
                chrom1yr = deterministic_signals.Deterministic(wf, name="chrom1yr")
                s += chrom1yr
                outstring += 'CHROM1YR, '

            if pm == "SW":
                n_earth = eparameter.Uniform(0, 20)
                deter_sw = solar_wind(n_earth=n_earth)
                mean_sw = deterministic_signals.Deterministic(deter_sw, name='n_earth')

                max_cadence = 60
                sw_components = 120

                log10_A_sw = eparameter.Uniform(-10, 1)
                gamma_sw = eparameter.Uniform(-4, 4)
                sw_prior = utils.powerlaw(log10_A=log10_A_sw, gamma=gamma_sw)
                sw_basis = createfourierdesignmatrix_solar_dm(nmodes=sw_components, Tspan=Tspan)

                sw = mean_sw + gp_signals.BasisGP(sw_prior, sw_basis, name='gp_sw')
                s += sw
                outstring += 'SWvary, '
            
            if pm == "SWDET":
                n_earth = eparameter.Uniform(0, 20)
                deter_sw = solar_wind(n_earth=n_earth)
                mean_sw = deterministic_signals.Deterministic(deter_sw, name='n_earth')

                sw = mean_sw

                s += sw
                outstring += 'SWdet, '

            if pm == "CHROMBUMP":
                chrom_gauss_bump = dm_gaussian_bump(tmin, tmax, idx="vary", sigma_min=604800, sigma_max=tmax-tmin, log10_A_low=-10, log10_A_high=-1, name='chrom_bump')
                s += chrom_gauss_bump
                outstring += 'CHROMBUMP, '
            
        if "SW" not in psrmodels and "SWDET" not in psrmodels:
            n_earth = eparameter.Constant(4)
            deter_sw = solar_wind(n_earth=n_earth)
            mean_sw = deterministic_signals.Deterministic(deter_sw, name='n_earth')

            sw = mean_sw
            s += sw
            outstring += 'SWfix, '
        

        for common_signal in s_common:
            s += common_signal
        models.append(s(psr))
        if verbose == True:
            print(outstring[:-2])
        

    if verbose == True:
        print(string_common)
    pta = signal_base.PTA(models)
    return pta
