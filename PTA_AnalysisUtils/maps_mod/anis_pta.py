import numpy as np
import scipy.optimize as sopt

import pickle, healpy as hp

import numpy.random as nr, scipy.stats as scst
from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc

from enterprise.signals import anis_coefficients as ac


import scipy.linalg as sl

from . import clebschGordan as CG, utils

from scipy.interpolate import interp1d
from astroML.linear_model import LinearRegression

import lmfit
from lmfit import minimize, Parameters

class anis_pta():

    def __init__(self, psrs_theta, psrs_phi,
                 xi = [], rho = [], sig = [], os = 1, os_cov = [],
                 l_max = 6, nside = 2,
                 mode = 'power_basis', use_physical_prior = False, include_pta_monopole = False):

        self.psrs_theta = psrs_theta
        self.psrs_phi = psrs_phi
        if len(xi) != 0:
            self.xi = xi
        else:
            self.xi = self._get_xi()

        #Read in OS and normalize cross-correlations by OS.
        self.os = os
        if len(rho) > 0 and len(sig) > 0:
            self.rho = np.array(rho) / self.os
            self.sig = np.array(sig) / self.os
            self.os_cov = os_cov / (self.os)**2
            
        else:
            self.rho = rho / self.os
            self.sig = sig / self.os
            self.os_cov = os_cov / (self.os)**2

        self.l_max = l_max
        self.nside = nside
        self.use_physical_prior = use_physical_prior
        self.include_pta_monopole = include_pta_monopole

        self.npsrs = len(self.psrs_theta)
        self.npairs = int(np.math.factorial(self.npsrs) / (np.math.factorial(2) * np.math.factorial(self.npsrs - 2)))

        self.npix = hp.nside2npix(self.nside)
        self.gw_theta, self.gw_phi = hp.pix2ang(nside = self.nside, ipix = np.arange(self.npix))

        self.mode = mode

        #Some configuration for spherical harmonic basis runs
        #clm refers to normal spherical harmonic basis
        #blm refers to sqrt power spherical harmonic basis
        self.blmax = int(self.l_max / 2.)
        self.clm_size = (self.l_max + 1) ** 2
        self.blm_size = hp.Alm.getsize(self.blmax)

        self.sqrt_basis_helper = CG.clebschGordan(l_max = self.l_max)
        #self.reorder, self.neg_idx, self.zero_idx, self.pos_idx = self.reorder_hp_ylm()

        if self.mode == 'power_basis':
            self.ndim = 1 + (self.l_max + 1) ** 2
        elif self.mode == 'sqrt_power_basis':
            #self.ndim = 1 + (2 * (hp.Alm.getsize(int(self.blmax)) - self.blmax))
            self.ndim = 1 + (self.blmax + 1) ** 2

        self.F_mat = self.antenna_response()

        if self.mode == 'power_basis' or self.mode == 'sqrt_power_basis':

            #The spherical harmonic basis for \Gamma_lm
            #shape (nclm, npsrs, npsrs)
            Gamma_lm = ac.anis_basis(np.dstack((self.psrs_phi, self.psrs_theta))[0], lmax = self.l_max,
                                     nside = self.nside)
            
            uti = np.triu_indices(n = self.npsrs, k = 1)

            self.Gamma_lm = np.full((Gamma_lm.shape[0], self.npairs), 0.0)

            for ii in range(Gamma_lm.shape[0]):

                self.Gamma_lm[ii] = Gamma_lm[ii][uti]
        
        elif self.mode == 'sph_basis':
            # The spherical harmonic basis for \Gamma_lm
            # shape (nclm, npsrs, npsrs)
            # changed --> use consistent polarisation basis with libstempo
            Gamma_lm = ac.anis_basis(np.dstack((self.psrs_phi, self.psrs_theta))[0], lmax = self.l_max,
                                     nside = self.nside, sph_complex = True, to_source=True)

            uti = np.triu_indices(n = self.npsrs, k = 1)

            self.Gamma_lm = np.full((Gamma_lm.shape[0], self.npairs), 0.0, dtype='complex_')

            for ii in range(Gamma_lm.shape[0]):

                self.Gamma_lm[ii] = Gamma_lm[ii][uti]
            print(np.shape(self.Gamma_lm))
        return None




    def _get_radec(self):

        psr_ra = []
        psr_dec = []

        for ptheta, pphi in zip(self.psrs_theta, self.psrs_phi):

            psr_ra.append(pphi)
            psr_dec.append(np.pi / 2 - ptheta)

        return psr_ra, psr_dec




    def _get_xi(self):

        psrs_ra, psrs_dec = self._get_radec()

        pos_vectors = np.array(
                [np.cos(psrs_ra) * np.cos(psrs_dec), np.sin(psrs_ra) * np.cos(psrs_dec), np.sin(psrs_dec)])

        xi = []

        for ii in range(len(self.psrs_theta)):
                for jj in range(ii+1, len(self.psrs_theta)):

                    xi.append(np.arccos(np.dot(pos_vectors[:, ii], pos_vectors[:, jj])))

        return np.array(xi)




    def _fplus_fcross(self, psrtheta, psrphi, gwtheta, gwphi):
        """
        Compute gravitational-wave quadrupolar antenna pattern.
        (From NX01)
        :param psr: pulsar object
        :param gwtheta: Polar angle of GW source in celestial coords [radians]
        :param gwphi: Azimuthal angle of GW source in celestial coords [radians]
        :returns: fplus, fcross
        """

        # define variable for later use
        cosgwtheta, cosgwphi = np.cos(gwtheta), np.cos(gwphi)
        singwtheta, singwphi = np.sin(gwtheta), np.sin(gwphi)

        # unit vectors to GW source
        m = np.array([singwphi, -cosgwphi, 0.0])
        n = np.array([-cosgwtheta*cosgwphi, -cosgwtheta*singwphi, singwtheta])
        omhat = np.array([-singwtheta*cosgwphi, -singwtheta*singwphi, -cosgwtheta])

        # pulsar location
        ptheta = psrtheta
        pphi = psrphi

        # use definition from Sesana et al 2010 and Ellis et al 2012
        phat = np.array([np.sin(ptheta)*np.cos(pphi), np.sin(ptheta)*np.sin(pphi),\
                np.cos(ptheta)])

        fplus = 0.5 * (np.dot(m, phat)**2 - np.dot(n, phat)**2) / (1 - np.dot(omhat, phat))
        fcross = (np.dot(m, phat)*np.dot(n, phat)) / (1 - np.dot(omhat, phat))

        return fplus, fcross




    def antenna_response(self):

        F_mat = np.zeros((self.npairs, self.npix))
        pair_no = 0

        for ii in range(self.npsrs):

            for jj in range(ii + 1, self.npsrs):

                for kk in range(self.npix):

                    pp_1, pc_1 = self._fplus_fcross(self.psrs_theta[ii], self.psrs_phi[ii],
                                              self.gw_theta[kk], self.gw_phi[kk])
                    pp_2, pc_2 = self._fplus_fcross(self.psrs_theta[jj], self.psrs_phi[jj],
                                              self.gw_theta[kk], self.gw_phi[kk])

                    F_mat[pair_no][kk] =  (pp_1 * pp_2 + pc_1 * pc_2)  * 1.5 / (self.npix)

                pair_no += 1

        return F_mat




    def get_pure_HD(self):
        #Return the theoretical HD curve given xi

        xx = (1 - np.cos(self.xi)) / 2.
        hd_curve = 1.5 * xx * np.log(xx) - xx / 4 + 0.5

        return hd_curve





    def orf_from_clm(self, params):
        #Using supplied clm values, calculate the corresponding power map
        #and calculate the ORF from that power map (convoluted, I know)

        amp2 = 10 ** params[0]
        clm = params[1:]

        sh_map = ac.mapFromClm(clm, nside = self.nside, sph_complex=True)

        orf = amp2 * np.dot(self.F_mat, sh_map)

        return np.longdouble(orf)


    

    def fisher_matrix_sph(self, cross_corr = False):

        if cross_corr == False:
            N_mat = np.zeros((len(self.rho), len(self.rho)))
            N_mat_inv = np.zeros((len(self.rho), len(self.rho)))
           
            N_mat[np.diag_indices(N_mat.shape[0])] = self.sig ** 2
            N_mat_inv[np.diag_indices(N_mat_inv.shape[0])] = 1 / self.sig ** 2
           
        else:
            N_mat = self.os_cov
            N_mat_inv = sl.inv(N_mat)
        F_mat_clm = self.Gamma_lm.transpose() # R-matrix

        fisher_mat = np.matmul(F_mat_clm.transpose().conj(), np.matmul(N_mat_inv, F_mat_clm))

        return fisher_mat

    

    def max_lkl_clm(self, cutoff = None, use_svd_reg = False, reg_type = 'l2',
                    alpha = 0, cross_corr = False, U_0=None, Vh_0=None):

        if cross_corr == False:
            N_mat = np.zeros((len(self.rho), len(self.rho)))
            N_mat_inv = np.zeros((len(self.rho), len(self.rho)))
            
            N_mat[np.diag_indices(N_mat.shape[0])] = self.sig ** 2
            N_mat_inv[np.diag_indices(N_mat_inv.shape[0])] = 1 / self.sig ** 2
            
        else:
            N_mat = self.os_cov
            N_mat_inv = sl.inv(N_mat)

        F_mat_clm = self.Gamma_lm.transpose() # R-matrix
        M = np.matmul(F_mat_clm.transpose().conj(), np.matmul(N_mat_inv, F_mat_clm))
        
        #SVD of the Fisher matrix
        U, sv, Vh = sl.svd(M, compute_uv = True)
        cn = np.max(sv) / np.min(sv)

        if use_svd_reg:
            if U_0 is not None and Vh_0 is not None:
                #calculate inverse fisher using V and U from iteration 0
                V_0 = Vh_0.transpose().conj()
                Uh_0 = U_0.transpose().conj()

                fisher0 = np.matmul(Vh_0, np.matmul(M, U_0))
                fisher0_segment = fisher0[:cutoff, :cutoff]
                # invert the segment
                invfisher0_segment = sl.inv(fisher0_segment)
                # fill upper left corner in empty matrix with inverted segment
                invfisher0 = np.zeros(np.shape(fisher0))
                invfisher0[:cutoff, :cutoff] = invfisher0_segment
                
                # transform back into original basis
                Minv = np.matmul(V_0, np.matmul(invfisher0, Uh_0))
                #print(Minv)

            else:
                if cutoff is not None:
                    abs_cutoff = sv[cutoff]
                    print('\t max_lkl_clm regularised at threshold {}, corresponding to sv {}'.format(sv[cutoff], cutoff) )
                else:
                    abs_cutoff = 0.
                    print('\t max_lkl_clm is not regularised')

                # inverse of regularised Fisher matrix
                Minv = sl.pinvh(M, cond = abs_cutoff)
                # X matrix
            
            X = np.matmul(F_mat_clm.transpose().conj(), np.matmul(N_mat_inv, self.rho))
            # clean map coefficients
            clms = np.matmul(Minv, X)

            clm_err = np.sqrt(np.diag(Minv))


        else:

            diag_identity = np.diag(np.full(F_mat_clm.shape[1], 1.0))

            fac1r = sl.pinvh(np.matmul(F_mat_clm.transpose(), np.matmul(N_mat_inv, F_mat_clm)) + alpha * diag_identity)

            clf = LinearRegression(regularization = reg_type, fit_intercept = False, kwds = dict(alpha = alpha))

            clf.fit(F_mat_clm, self.rho, self.sig)

            clms = clf.coef_

            clm_err = np.sqrt(np.diag(fac1r))

        return clms, clm_err, cn, sv, Minv, X, U, Vh

    
