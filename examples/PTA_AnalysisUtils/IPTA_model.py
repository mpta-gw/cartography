#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 09:53:48 2023

@author: kgrunthal
"""

# This script runs enterprise for individual pulsars, or for an entire gravitational wave source

from __future__ import division
import numpy as np



import enterprise.signals.parameter as parameter
from enterprise.signals import utils
from enterprise.signals import signal_base
from enterprise.signals import selections
from enterprise.signals import white_signals
from enterprise.signals import gp_signals
from enterprise.signals import deterministic_signals
from enterprise.signals import gp_bases
import enterprise.constants as const
from enterprise_extensions.chromatic.solar_wind import solar_wind, createfourierdesignmatrix_solar_dm
from enterprise_extensions import deterministic
from enterprise_extensions import chromatic

from enterprise_extensions import model_utils, blocks
from enterprise_extensions import timing





### Auxiliary functions ###

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
    sign_param = parameter.Uniform(-1,1)
    t0_dm_bump = parameter.Uniform(tmin[0],tmax[0])
    sigma_dm_bump = parameter.Uniform(sigma_min,tmax[0]-tmin[0])
    log10_Amp_dm_bump = parameter.Uniform(log10_A_low, log10_A_high)
    if idx == 'vary':
        idx = parameter.Uniform(0, 14)
    wf = chrom_gaussian_bump(log10_Amp=log10_Amp_dm_bump,
                         t0=t0_dm_bump, sigma=sigma_dm_bump,
                         sign_param=sign_param, idx=idx)
    dm_bump = deterministic_signals.Deterministic(wf, name=name)

    return dm_bump




#------------------------------------------------------------------------------
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






#------------------------------------------------------------------------------
@signal_base.function
def cw_delay(toas, pos, pdist,
             cos_gwtheta=0, gwphi=0, cos_inc=0,
             log10_mc=9, log10_fgw=-8, log10_dist=None, log10_h=None,
             phase0=0, psi=0, k_drop=1.,
             psrTerm=False, p_dist=1, p_phase=None,
             evolve=False, phase_approx=False, check=False,
             tref=0):
    """
    Function to create GW incuced residuals from a SMBMB as
    defined in Ellis et. al 2012,2013.
    :param toas:
        Pular toas in seconds
    :param pos:
        Unit vector from the Earth to the pulsar
    :param pdist:
        Pulsar distance (mean and uncertainty) [kpc]
    :param cos_gwtheta:
        Cosine of Polar angle of GW source in celestial coords [radians]
    :param gwphi:
        Azimuthal angle of GW source in celestial coords [radians]
    :param cos_inc:
        cosine of Inclination of GW source [radians]
    :param log10_mc:
        log10 of Chirp mass of SMBMB [solar masses]
    :param log10_fgw:
        log10 of Frequency of GW (twice the orbital frequency) [Hz]
    :param log10_dist:
        log10 of Luminosity distance to SMBMB [Mpc],
        used to compute strain, if not None
    :param log10_h:
        log10 of GW strain,
        used to compute distance, if not None
    :param phase0:
        Initial GW phase of source [radians]
    :param psi:
        Polarization angle of GW source [radians]
    :param psrTerm:
        Option to include pulsar term [boolean]
    :param p_dist:
        Pulsar distance parameter
    :param p_phase:
        Use pulsar phase to determine distance [radian]
    :param evolve:
        Option to include/exclude full evolution [boolean]
    :param phase_approx:
        Option to include/exclude phase evolution across observation time
        [boolean]
    :param check:
        Check if frequency evolves significantly over obs. time [boolean]
    :param tref:
        Reference time for phase and frequency [s]
    :return: Vector of induced residuals
    """

    # convert units to time
    mc = 10**log10_mc * const.Tsun
    fgw = 10**log10_fgw
    gwtheta = np.arccos(cos_gwtheta)
    inc = np.arccos(cos_inc)
    p_dist = (pdist[0] + pdist[1]*p_dist)*const.kpc/const.c

    if log10_h is None and log10_dist is None:
        raise ValueError("one of log10_dist or log10_h must be non-None")
    elif log10_h is not None and log10_dist is not None:
        raise ValueError("only one of log10_dist or log10_h can be non-None")
    elif log10_h is None:
        dist = 10**log10_dist * const.Mpc / const.c
    else:
        dist = 2 * mc**(5/3) * (np.pi*fgw)**(2/3) / 10**log10_h

    if check:
        # check that frequency is not evolving significantly over obs. time
        fstart = fgw * (1 - 256/5 * mc**(5/3) * fgw**(8/3) * toas[0])**(-3/8)
        fend = fgw * (1 - 256/5 * mc**(5/3) * fgw**(8/3) * toas[-1])**(-3/8)
        df = fend - fstart

        # observation time
        Tobs = toas.max()-toas.min()
        fbin = 1/Tobs

        if np.abs(df) > fbin:
            print('WARNING: Frequency is evolving over more than one '
                  'frequency bin.')
            print('f0 = {0}, f1 = {1}, df = {2}, fbin = {3}'.format(fstart, fend, df, fbin))
            return np.ones(len(toas)) * np.nan

    # get antenna pattern funcs and cosMu
    # write function to get pos from theta,phi
    fplus, fcross, cosMu = utils.create_gw_antenna_pattern(pos, gwtheta, gwphi)

    # get pulsar time
    toas -= tref
    if p_dist > 0:
        tp = toas-p_dist*(1-cosMu)
    else:
        tp = toas

    # orbital frequency
    w0 = np.pi * fgw
    phase0 /= 2  # convert GW to orbital phase
    # omegadot = 96/5 * mc**(5/3) * w0**(11/3) # Not currently used in code

    # evolution
    if evolve:
        # calculate time dependent frequency at earth and pulsar
        omega = w0 * (1 - 256/5 * mc**(5/3) * w0**(8/3) * toas)**(-3/8)
        omega_p = w0 * (1 - 256/5 * mc**(5/3) * w0**(8/3) * tp)**(-3/8)

        if p_dist > 0:
            omega_p0 = w0 * (1 + 256/5
                             * mc**(5/3) * w0**(8/3) * p_dist*(1-cosMu))**(-3/8)
        else:
            omega_p0 = w0

        # calculate time dependent phase
        phase = phase0 + 1/32/mc**(5/3) * (w0**(-5/3) - omega**(-5/3))

        if p_phase is None:
            phase_p = phase0 + 1/32/mc**(5/3) * (w0**(-5/3) - omega_p**(-5/3))
        else:
            phase_p = (phase0 + p_phase
                       + 1/32*mc**(-5/3) * (omega_p0**(-5/3) - omega_p**(-5/3)))

    elif phase_approx:
        # monochromatic
        omega = w0
        if p_dist > 0:
            omega_p = w0 * (1 + 256/5
                            * mc**(5/3) * w0**(8/3) * p_dist*(1-cosMu))**(-3/8)
        else:
            omega_p = w0

        # phases
        phase = phase0 + omega * toas
        if p_phase is not None:
            phase_p = phase0 + p_phase + omega_p*toas
        else:
            phase_p = (phase0 + omega_p*toas
                       + 1/32/mc**(5/3) * (w0**(-5/3) - omega_p**(-5/3)))

    # no evolution
    else:
        # monochromatic
        omega = np.pi*fgw
        omega_p = omega

        # phases
        phase = phase0 + omega * toas
        phase_p = phase0 + omega * tp

    # define time dependent coefficients
    At = -0.5*np.sin(2*phase)*(3+np.cos(2*inc))
    Bt = 2*np.cos(2*phase)*np.cos(inc)
    At_p = -0.5*np.sin(2*phase_p)*(3+np.cos(2*inc))
    Bt_p = 2*np.cos(2*phase_p)*np.cos(inc)

    # now define time dependent amplitudes
    alpha = mc**(5./3.)/(dist*omega**(1./3.))
    alpha_p = mc**(5./3.)/(dist*omega_p**(1./3.))

    # define rplus and rcross
    rplus = alpha*(-At*np.cos(2*psi)+Bt*np.sin(2*psi))
    rcross = alpha*(At*np.sin(2*psi)+Bt*np.cos(2*psi))
    rplus_p = alpha_p*(-At_p*np.cos(2*psi)+Bt_p*np.sin(2*psi))
    rcross_p = alpha_p*(At_p*np.sin(2*psi)+Bt_p*np.cos(2*psi))

    # k dropout
    k = np.rint(k_drop)

    # residuals
    if psrTerm:
        res = k*fplus*(rplus-rplus_p) + k*fcross*(rcross-rcross_p)
    else:
        res = k*fplus*rplus + k*fcross*rcross

    return res


@signal_base.function
def compute_eccentric_residuals(toas, theta, phi, cos_gwtheta, gwphi,
                                log10_mc, log10_dist, log10_h, log10_F, cos_inc,
                                psi, gamma0, e0, l0, q, nmax=400, pdist=1.0,
                                pphase=None, pgam=None, psrTerm=False,
                                tref=0, check=False):
    """
    Simulate GW from eccentric SMBHB. Waveform models from
    Taylor et al. (2015) and Barack and Cutler (2004).
    WARNING: This residual waveform is only accurate if the
    GW frequency is not significantly evolving over the
    observation time of the pulsar.
    :param toa: pulsar observation times
    :param theta: polar coordinate of pulsar
    :param phi: azimuthal coordinate of pulsar
    :param gwtheta: Polar angle of GW source in celestial coords [radians]
    :param gwphi: Azimuthal angle of GW source in celestial coords [radians]
    :param log10_mc: Base-10 lof of chirp mass of SMBMB [solar masses]
    :param log10_dist: Base-10 uminosity distance to SMBMB [Mpc]
    :param log10_F: base-10 orbital frequency of SMBHB [Hz]
    :param inc: Inclination of GW source [radians]
    :param psi: Polarization of GW source [radians]
    :param gamma0: Initial angle of periastron [radians]
    :param e0: Initial eccentricity of SMBHB
    :param l0: Initial mean anomoly [radians]
    :param q: Mass ratio of SMBHB
    :param nmax: Number of harmonics to use in waveform decomposition
    :param pdist: Pulsar distance [kpc]
    :param pphase: Pulsar phase [rad]
    :param pgam: Pulsar angle of periastron [rad]
    :param psrTerm: Option to include pulsar term [boolean]
    :param tref: Fidicuial time at which initial parameters are referenced [s]
    :param check: Check if frequency evolves significantly over obs. time
    :returns: Vector of induced residuals
    """

    # convert from sampling
    F = 10.0**log10_F
    mc = 10.0**log10_mc
    dist = 10.0**log10_dist
    if log10_h is not None:
        h0 = 10.0**log10_h
    else:
        h0 = None
    inc = np.arccos(cos_inc)
    gwtheta = np.arccos(cos_gwtheta)

    # define variable for later use
    cosgwtheta, cosgwphi = np.cos(gwtheta), np.cos(gwphi)
    singwtheta, singwphi = np.sin(gwtheta), np.sin(gwphi)
    sin2psi, cos2psi = np.sin(2*psi), np.cos(2*psi)

    # unit vectors to GW source
    m = np.array([singwphi, -cosgwphi, 0.0])
    n = np.array([-cosgwtheta*cosgwphi, -cosgwtheta*singwphi, singwtheta])
    omhat = np.array([-singwtheta*cosgwphi, -singwtheta*singwphi, -cosgwtheta])

    # pulsar position vector
    phat = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi),
                     np.cos(theta)])

    fplus = 0.5 * (np.dot(m, phat)**2 - np.dot(n, phat)**2) / (1+np.dot(omhat, phat))
    fcross = (np.dot(m, phat)*np.dot(n, phat)) / (1 + np.dot(omhat, phat))
    cosMu = -np.dot(omhat, phat)

    # get values from pulsar object
    toas = toas.copy() - tref

    if check:
        # check that frequency is not evolving significantly over obs. time
        y = utils.solve_coupled_ecc_solution(F, e0, gamma0, l0, mc, q,
                                             np.array([0.0, toas.max()]))

        # initial and final values over observation time
        Fc0, ec0, gc0, phic0 = y[0, :]
        Fc1, ec1, gc1, phic1 = y[-1, :]

        # observation time
        Tobs = 1/(toas.max()-toas.min())

        if np.abs(Fc0-Fc1) > 1/Tobs:
            print('WARNING: Frequency is evolving over more than one frequency bin.')
            print('F0 = {0}, F1 = {1}, delta f = {2}'.format(Fc0, Fc1, 1/Tobs))
            return np.ones(len(toas)) * np.nan

    # get gammadot for earth term
    gammadot = utils.get_gammadot(F, mc, q, e0)

    # get number of harmonics to use
    if not isinstance(nmax, int):
        if isinstance(nmax, str):
            f1yr = 1/(3600*24*365.25)
            nharm = int(f1yr/F)
        elif e0 < 0.999 and e0 > 0.001:
            nharm = int(nmax(e0))
        elif e0 < 0.001:
            nharm = 2
        else:
            nharm = int(nmax(0.999))
    else:
        nharm = nmax

    # no more than 100 harmonics
    nharm = min(nharm, 100)

    ##### earth term #####
    splus, scross = utils.calculate_splus_scross(nmax=nharm, mc=mc, dl=dist,
                                                 h0=h0, F=F, e=e0, t=toas.copy(),
                                                 l0=l0, gamma=gamma0,
                                                 gammadot=gammadot, inc=inc)

    ##### pulsar term #####
    if psrTerm:
        # pulsar distance
        pd = pdist

        # convert units
        pd *= const.kpc / const.c

        # get pulsar time
        tp = toas.copy() - pd * (1-cosMu)

        # solve coupled system of equations to get pulsar term values
        y = utils.solve_coupled_ecc_solution(F, e0, gamma0, l0, mc,
                                             q, np.array([0.0, tp.min()]))

        # get pulsar term values
        if np.any(y):
            Fp, ep, gp, phip = y[-1, :]

            # get gammadot at pulsar term
            gammadotp = utils.get_gammadot(Fp, mc, q, ep)

            # get phase at pulsar
            if pphase is None:
                lp = phip
            else:
                lp = pphase

            # get angle of periastron at pulsar
            if pgam is None:
                gp = gp
            else:
                gp = pgam

            # get number of harmonics to use
            if not isinstance(nmax, int):
                if e0 < 0.999 and e0 > 0.001:
                    nharm = int(nmax(e0))
                elif e0 < 0.001:
                    nharm = 2
                else:
                    nharm = int(nmax(0.999))
            else:
                nharm = nmax

            # no more than 1000 harmonics
            nharm = min(nharm, 100)
            splusp, scrossp = utils.calculate_splus_scross(nmax=nharm, mc=mc,
                                                           dl=dist, h0=h0,
                                                           F=Fp, e=ep,
                                                           t=toas.copy(),
                                                           l0=lp, gamma=gp,
                                                           gammadot=gammadotp,
                                                           inc=inc)

            rr = (fplus*cos2psi - fcross*sin2psi) * (splusp - splus) + \
                (fplus*sin2psi + fcross*cos2psi) * (scrossp - scross)

        else:
            rr = np.ones(len(toas)) * np.nan

    else:
        rr = - (fplus*cos2psi - fcross*sin2psi) * splus - \
            (fplus*sin2psi + fcross*cos2psi) * scross

    return rr

# Extra model components not part of base enterprise ####
def cw_block_circ(amp_prior='log-uniform', dist_prior=None, drop=False, drop_psr=False,
                  skyloc=None, log10_fgw=None, patch=None,
                  psrTerm=False, tref=0, name='cw'):
    """
    Returns deterministic, cirular orbit continuous GW model:
    :param amp_prior:
        Prior on log10_h. Default is "log-uniform."
        Use "uniform" for upper limits, or "None" to search over
        log10_dist instead.
    :param dist_prior:
        Prior on log10_dist. Default is "None," meaning that the
        search is over log10_h instead of log10_dist. Use "log-uniform"
        to search over log10_h with a log-uniform prior.
    :param skyloc:
        Fixed sky location of CW signal search as [cos(theta), phi].
        Search over sky location if ``None`` given.
    :param log10_fgw:
        Fixed log10 GW frequency of CW signal search.
        Search over GW frequency if ``None`` given.
    :param patch:
        Bounds for sky position prior [[costh_min, costh_max], [phi_min, phi_max]]
    :param ecc:
        Fixed log10 distance to SMBHB search.
        Search over distance or strain if ``None`` given.
    :param psrTerm:
        Boolean for whether to include the pulsar term. Default is False.
    :param name:
        Name of CW signal.
    """

    if dist_prior == None:
        log10_dist = None

        if amp_prior == 'uniform':
            log10_h = parameter.LinearExp(-16.0, -11.0)('{}_log10_h'.format(name))
        elif amp_prior == 'log-uniform':
            log10_h = parameter.Uniform(-16.0, -11.0)('{}_log10_h'.format(name))

    elif dist_prior == 'log-uniform':
        log10_dist = parameter.Uniform(-2.0, 4.0)('{}_log10_dL'.format(name))
        log10_h = None

    # chirp mass [Msol]
    log10_Mc = parameter.Uniform(6.0, 10.0)('{}_log10_Mc'.format(name))

    # GW frequency [Hz]
    if log10_fgw is None:
        log10_fgw = parameter.Uniform(-9.0, -7.0)('{}_log10_fgw'.format(name))
    elif isinstance(log10_fgw, list):
        log10_fgw = parameter.Uniform(log10_fgw[0], log10_fgw[1])('{}_log10_fgw'.format(name))
    else:
        log10_fgw = parameter.Constant(log10_fgw)('{}_log10_fgw'.format(name))
    # orbital inclination angle [radians]
    cosinc = parameter.Uniform(-1.0, 1.0)('{}_cosinc'.format(name))
    # initial GW phase [radians]
    phase0 = parameter.Uniform(0.0, np.pi)('{}_phase0'.format(name))

    # polarization
    psi_name = '{}_psi'.format(name)
    psi = parameter.Uniform(0, np.pi)(psi_name)

    # sky location
    costh_name = '{}_costheta'.format(name)
    phi_name = '{}_phi'.format(name)
    if skyloc is None:
        if patch is None:
            costh = parameter.Uniform(-1, 1)(costh_name)
            phi = parameter.Uniform(0, 2*np.pi)(phi_name)
        else:
            costh = parameter.Uniform(patch[0, 0], patch[0, 1])(costh_name)
            phi = parameter.Uniform(patch[1, 0], patch[1, 1])(phi_name)
    else:
        costh = parameter.Constant(skyloc[0])(costh_name)
        phi = parameter.Constant(skyloc[1])(phi_name)

    if psrTerm:
        p_phase = parameter.Uniform(0, 2*np.pi)
        p_dist = parameter.Normal(0, 1)
    else:
        p_phase = None
        p_dist = 0

    if drop:
        if drop_psr:
            k_drop = parameter.Uniform(0, 1)
        else:
            kdrop_name = '{}_k_drop'.format(name)
            k_drop = parameter.Uniform(0, 1)(kdrop_name)
    else:
        k_drop = 1.

    # continuous wave signal
    wf = cw_delay(cos_gwtheta=costh, gwphi=phi, cos_inc=cosinc,
                  log10_mc=log10_Mc, log10_fgw=log10_fgw,
                  log10_h=log10_h, log10_dist=log10_dist,
                  phase0=phase0, psi=psi, k_drop=k_drop,
                  psrTerm=psrTerm, p_dist=p_dist, p_phase=p_phase,
                  phase_approx=True, check=False,
                  tref=tref)
    cw = deterministic.CWSignal(wf, ecc=False, psrTerm=psrTerm)

    return cw

def cw_block_ecc(amp_prior='log-uniform', skyloc=None, log10_F=None, patch=None, nmax='1yr',
                 ecc=None, psrTerm=False, tref=0, name='cw'):
    """
    Returns deterministic, eccentric orbit continuous GW model:
    :param amp_prior:
        Prior on log10_h and log10_Mc/log10_dL. Default is "log-uniform" with
        log10_Mc and log10_dL searched over. Use "uniform" for upper limits,
        log10_h searched over.
    :param skyloc:
        Fixed sky location of CW signal search as [cos(theta), phi].
        Search over sky location if ``None`` given.
    :param log10_F:
        Fixed log-10 orbital frequency of CW signal search.
        Search over orbital frequency if ``None`` given.
    :param patch:
        Bounds for sky position prior [[costh_min, costh_max], [phi_min, phi_max]]
    :param ecc:
        Fixed eccentricity of SMBHB search.
        Search over eccentricity if ``None`` given.
    :param psrTerm:
        Boolean for whether to include the pulsar term. Default is False.
    :param name:
        Name of CW signal.
    """

    if amp_prior == 'uniform':
        log10_h = parameter.LinearExp(-18.0, -11.0)('{}_log10_h'.format(name))
    elif amp_prior == 'log-uniform':
        log10_h = parameter.Uniform(-18.0, -11.0)('{}_log10_h'.format(name))
    # chirp mass [Msol]
    log10_Mc = parameter.Uniform(6.0, 10.0)('{}_log10_Mc'.format(name))
    # luminosity distance [Mpc]
    log10_dL = parameter.Uniform(-2.0, 4.0)('{}_log10_dL'.format(name))

    # orbital frequency [Hz]
    if log10_F is None:
        log10_Forb = parameter.Uniform(-9.0, -7.0)('{}_log10_Forb'.format(name))
    elif isinstance(log10_F, list):
        log10_Forb = parameter.Uniform(log10_F[0], log10_F[1])('{}_log10_Forb'.format(name))
    else:
        log10_Forb = parameter.Constant(log10_F)('{}_log10_Forb'.format(name))
    # orbital inclination angle [radians]
    cosinc = parameter.Uniform(-1.0, 1.0)('{}_cosinc'.format(name))
    # periapsis position angle [radians]
    gamma_0 = parameter.Uniform(0.0, np.pi)('{}_gamma0'.format(name))

    # Earth-term eccentricity
    if ecc is None:
        e_0 = parameter.Uniform(0.0, 0.99)('{}_e0'.format(name))
    elif isinstance(ecc, list):
        e_0 = parameter.Uniform(ecc[0], ecc[1])('{}_e0'.format(name))
    else:
        e_0 = parameter.Constant(ecc)('{}_e0'.format(name))

    # initial mean anomaly [radians]
    l_0 = parameter.Uniform(0.0, 2.0*np.pi)('{}_l0'.format(name))
    # mass ratio = M_2/M_1
    q = parameter.Constant(1.0)('{}_q'.format(name))

    # polarization
    pol_name = '{}_pol'.format(name)
    pol = parameter.Uniform(0, np.pi)(pol_name)

    # sky location
    costh_name = '{}_costheta'.format(name)
    phi_name = '{}_phi'.format(name)
    if skyloc is None:
        if patch is None:
            costh = parameter.Uniform(-1, 1)(costh_name)
            phi = parameter.Uniform(0, 2*np.pi)(phi_name)
        else:
            costh = parameter.Uniform(patch[0, 0], patch[0, 1])(costh_name)
            phi = parameter.Uniform(patch[1, 0], patch[1, 1])(phi_name)
    else:
        costh = parameter.Constant(skyloc[0])(costh_name)
        phi = parameter.Constant(skyloc[1])(phi_name)

    # continuous wave signal
    wf = compute_eccentric_residuals(cos_gwtheta=costh, gwphi=phi,
                                     log10_mc=log10_Mc, log10_dist=log10_dL,
                                     log10_h=log10_h, log10_F=log10_Forb,
                                     cos_inc=cosinc, psi=pol, gamma0=gamma_0,
                                     e0=e_0, l0=l_0, q=q, nmax=nmax,
                                     pdist=None, pphase=None, pgam=None,
                                     tref=tref, check=False)
    cw = deterministic.CWSignal(wf, ecc=True, psrTerm=psrTerm)

    return cw







################################################################
################################################################
################################################################


#Because 1713 is a special snowflake
dm_expdip = True
num_dmdips = 2
dm_expdip_tmin = [54500, 57300]
dm_expdip_tmax = [54900, 57800]
dm_expdip_idx = [2, 2]
dmexp_sign = 'negative'


def PTA(psrs, 
        noisemodels, noise_dict,
        wn_vary=False, wideband=False):
    
    
    # compute tspan
    Tspan = model_utils.get_tspan(psrs)
    
    tmin_glob = [p.toas.min() for p in psrs]
    tmax_glob = [p.toas.max() for p in psrs]
    
    
    models = [] # store the model for each pulsar in here
    s_common = []
    string_common = 'PTA common signal:'
    
    # define the signal common to all pulsars
    if "common" in noisemodels.keys():
        
        if "CRN" in noisemodels["common"].keys():
            curn_dict = noisemodels["common"]["CRN"]
            
            CRN = blocks.common_red_noise_block(psd=curn_dict["psd"],
                                                 components=int(curn_dict["ncomp"]),
                                                 gamma_val=curn_dict["gamma"],
                                                 orf=curn_dict["orf"]
                                                 )
            string_common += ' CRN -- ORF: {}; PSD: {}; comp: {}; gamma {}) '.format(curn_dict["orf"], curn_dict["psd"], curn_dict["ncomp"], curn_dict["gamma"])
                
            s_common.append(CRN)
                
            
        
        
            
    
        elif 'CGW' in noisemodels["common"].keys():
            '''
            amp_prior = 'uniform' if upper_limit else 'log-uniform'
            prefixes = [prefix]
            for ngw in range(n_cgw-1):
                prefixes.append(prefix+'_'+str(ngw+1))
            CGW = []
            for ngw in range(n_cgw):
                if not ecc:
                    CGW = cw_block_circ(amp_prior=amp_prior, patch=patch,
                            skyloc=skyloc, log10_fgw=log10_F, drop=drop, drop_psr=drop_psr,
                            psrTerm=psrTerm, tref=tmin_glob, name=prefixes[ngw])
                    s_common.append(CGW)
                    del CGW
 
                else:
                    if type(ecc) is not float:
                        e0 = None
                    if isinstance(ecc, list):
                        e0 = ecc
                    ECGW = cw_block_ecc(amp_prior=amp_prior, patch=patch,
                                    skyloc=skyloc, log10_F=log10_F, ecc=e0,
                                    psrTerm=psrTerm, tref=tmin_glob, name=prefixes[ngw])
                    s_common.append(ECGW)
                    del ECGW
            '''
            cos_gwtheta = parameter.Uniform(-1, 1)('cos_gwtheta')   # position of source
            gwphi = parameter.Uniform(0, 2*np.pi)('gwphi')          # position of source
            log10_mc = parameter.Uniform(6.5, 10.)('log10_Mc')      # chirp mass
            log10_h = parameter.Uniform(-16, -11)('log10_h')        # strain amplitude
            log10_fgw = parameter.Uniform(-9.0, -7.0)('log10_fgw')    # gw frequency
            phase0 = parameter.Uniform(0, 2*np.pi)('phase0')        # gw phase
            psi = parameter.Uniform(0, np.pi)('psi')                # gw polarization 
            cos_inc = parameter.Uniform(-1, 1)('cos_inc')           # inclination of binary with respect to Earth
            
            cw_wf = deterministic.cw_delay(cos_gwtheta=cos_gwtheta, gwphi=gwphi, log10_mc=log10_mc, 
                                           log10_h=log10_h, log10_fgw=log10_fgw, phase0=phase0, 
                                           psi=psi, cos_inc=cos_inc)
            CGW = deterministic.CWSignal(cw_wf, ecc=False, psrTerm=False)
            s_common.append(CGW)
            string_common += ' CGW '
        
        
        
        elif 'EPH' in noisemodels["common"].keys():
            EPH = deterministic_signals.PhysicalEphemerisSignal()
            s_common.append(EPH)
            string_common += ' EPH '
            
            
            
        else:
            string_common += "No common signal given"

    
    
    # individual model        
    for psr in psrs:
        
        # get the individual models for each pulsar
        s = [x for x in s_common]
        SPNA_dict = noisemodels[psr.name]
        psr_modelcomps = SPNA_dict.keys()   
        
        string = psr.name + ': '
        
        
        # helper variables
        tmin = psr.toas.min()
        tmax = psr.toas.max()
        Tspan=(tmax-tmin)
        freqs = np.linspace(1/Tspan,30/Tspan,30)
        
        
        for comp in psr_modelcomps:
            if comp == 'TM':
                #TM = gp_signals.MarginalizingTimingModel(use_svd=True)
                TM = gp_signals.TimingModel(use_svd=True)
                s.append(TM)
                string += 'TM\t'
            
            elif comp == 'WN':
                if SPNA_dict['WN'] == 'fix':
                    WN = blocks.white_noise_block(vary=False, inc_ecorr=False,
                                           gp_ecorr=False, tnequad=True)
                elif SPNA_dict['WN'] == 'n':
                    WN = white_signals.MeasurementNoise(efac = parameter.Constant(1.))
                else:
                    WN = blocks.white_noise_block(vary=True, inc_ecorr=False,
                                       gp_ecorr=False, tnequad=True)
                                    
                s.append(WN)
                string += 'WN\t'

            
            elif comp == 'RN':
                RN = blocks.red_noise_block(psd = SPNA_dict['RN']['psd'], components = int(SPNA_dict['RN']['ncomp']))
                s.append(RN)
                string += 'RN\t'

        
            elif comp == 'DM':
                log10_A_dm = parameter.Uniform(-20, -11)
                
                if "wideband" in SPNA_dict['DM'].keys():
                    gamma_dm = parameter.Uniform(0, 14)
                    string += 'DM_WIDE\t'
                else:
                    gamma_dm = parameter.Uniform(0, 7)
                    string += 'DM\t'
                
                DM = blocks.dm_noise_block(log10_A=log10_A_dm,gamma=gamma_dm,Tspan=Tspan,
                              components=int(SPNA_dict['DM']['ncomp']), option=SPNA_dict['DM']['psd'])
                s.append(DM)
            
            
            elif comp == 'DMEXP':
                tmin = (SPNA_dict['DMEXP']['tmin'] if isinstance(SPNA_dict['DMEXP']['tmin'], list)
                    else [SPNA_dict['DMEXP']['tmin']])
                tmax = (SPNA_dict['DMEXP']['tmax'] if isinstance(SPNA_dict['DMEXP']['tmax'], list)
                    else [SPNA_dict['DMEXP']['tmax']])
            
                dmdipname_base = ['dmexp_{0}'.format(ii+1)
                                    for ii in range(SPNA_dict['DMEXP']['n_dips'])]
            
                dm_expdip_idx = (SPNA_dict['DMEXP']['idx'] if isinstance(SPNA_dict['DMEXP']['idx'], list)
                                 else [SPNA_dict['DMEXP']['idx']])
                
                for dd in range(SPNA_dict['DMEXP']['n_dips']):
                    DMEXP = chromatic.dm_exponential_dip(tmin=tmin[dd],
                                                      tmax=tmax[dd],
                                                      idx=dm_expdip_idx[dd],
                                                      sign=SPNA_dict['DMEXP']['sign'],
                                                      name=dmdipname_base[dd])
                    s.append(DMEXP)
                    del DMEXP
                string += '{} exp dip\t'.format(num_dmdips)
            
        
            elif comp == 'CHROM':
                log10_A_chrom_prior = parameter.Uniform(-20, -11)
                
                if "wideband" in SPNA_dict['CHROM'].keys() and "idx" in SPNA_dict['CHROM'].keys():
                    gamma_chrom_prior = parameter.Uniform(0, 14)
                    chrom_gp_idx = parameter.Constant(float(SPNA_dict['CHROM']['idx']))
                    string += 'CHROM_WIDE,idx4\t'
                    
                elif "wideband" in SPNA_dict['CHROM'].keys() and "idx" not in SPNA_dict['CHROM'].keys():
                    gamma_chrom_prior = parameter.Uniform(0, 14)
                    chrom_gp_idx = parameter.Uniform(0,14)
                    string += 'CHROM_WIDE,free\t'
                    
                elif "wideband" not in SPNA_dict['CHROM'].keys() and "idx" in SPNA_dict['CHROM'].keys():
                    gamma_chrom_prior = parameter.Uniform(0, 7)
                    chrom_gp_idx = parameter.Constant(float(SPNA_dict['CHROM']['idx']))
                    string += 'CHROM,idx4\t'
                    
                else:
                    gamma_chrom_prior = parameter.Uniform(0, 7)
                    chrom_gp_idx = parameter.Uniform(0,14)
                    string += 'CHROM,free\t'
                                    
                chrom_model = utils.powerlaw(log10_A=log10_A_chrom_prior, gamma=gamma_chrom_prior)
                chrom_basis = gp_bases.createfourierdesignmatrix_chromatic(nmodes=int(SPNA_dict['CHROM']['nb']),
                                                                           idx=chrom_gp_idx)
                CHROM = gp_signals.BasisGP(chrom_model, chrom_basis, name='chrom_gp')
                s.append(CHROM)

  

            elif comp == 'CHROMANNUAL':
                log10_Amp_chrom1yr = parameter.Uniform(-20, -5)
                phase_chrom1yr = parameter.Uniform(0, 2*np.pi)
                idx_chrom1yr = parameter.Uniform(0, 14)
                wf = chrom_yearly_sinusoid(log10_Amp=log10_Amp_chrom1yr, phase=phase_chrom1yr, idx=idx_chrom1yr)
                CHROMANNUAL = deterministic_signals.Deterministic(wf, name="chrom1yr")
                s.append(CHROMANNUAL)
                string += 'ChromAnnual\t'


            elif comp=='CHROMBUMP':
                CHROMBUMP = dm_gaussian_bump(tmin_glob, tmax_glob, idx="vary", sigma_min=604800,
                                         sigma_max=tmax_glob[0]-tmin_glob[0], log10_A_low=-10, log10_A_high=-1)
                s.append(CHROMBUMP)
                string += 'ChromBump\t'


        
            elif comp=='SW' :
                if SPNA_dict['SW']['type']=='det':
                    n_earth = parameter.Uniform(0, 20)
                    deter_sw = solar_wind(n_earth=n_earth)
                    SWDET = deterministic_signals.Deterministic(deter_sw, name='n_earth')
                    s.append(SWDET)
                    string += 'SW,det\t'
                    
                elif SPNA_dict['SW']['type']=='vary':
                    n_earth = parameter.Uniform(0, 20)
                    log10_A_sw = parameter.Uniform(-10, 1)
                    gamma_sw = parameter.Uniform(-4, 4)
                    Tspan = psr.toas.max() - psr.toas.min()
                    deter_sw = solar_wind(n_earth=n_earth)
                    mean_sw = deterministic_signals.Deterministic(deter_sw, name='n_earth')
                    sw_prior = utils.powerlaw(log10_A=log10_A_sw, gamma=gamma_sw)
                    sw_basis = createfourierdesignmatrix_solar_dm(nmodes=int(SPNA_dict['SW']['nb']), Tspan=Tspan)
                    SWVARY = mean_sw + gp_signals.BasisGP(sw_prior, sw_basis, name='gp_sw')
                    s.append(SWVARY)
                    string += 'SW,det+stoch\t'
                
                elif SPNA_dict['SW']['type']=='const':
                    n_earth = parameter.Constant(SPNA_dict['SW']['n_earth'])
                    deter_sw = solar_wind(n_earth=n_earth)
                    SWDET = deterministic_signals.Deterministic(deter_sw, name='sw')
                    s.append(SWDET)
                    string += 'SW,const\t'
        
                else:
                    string += 'unknown SW model'
            
            else:
                string += 'unknown SPNA model'
                    
        
        full_model = s[0]
        for pt in s[1:]:
            full_model += pt
        
        models.append(full_model(psr))
        
        del full_model
        del s
        
        print(string, flush=True)
    # put all pulsar models together into a PTA
    print(string_common)
    print(50*'-', '\n')

    pta = signal_base.PTA(models)
    
    return pta




            

