# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 15:04:10 2026

@author: kgrun
"""


import os
import numpy as np
from scipy.special import lpmv, factorial

cartography_path = os.path.dirname(__file__) +"/"


def lm_to_index(l, m):
    """0-based index of (l,m) in PTA ordering: l^2 + l + m"""
    return l * l + l + m




def _clm2clmreal(clm):
    """
    Convert complex SpH coefficient vector (PTA ordering) to real lmcosi array.
    
    Input
    -----
    clm : 1D complex array, length (lmax+1)**2, PTA ordering.
    
    Output
    ------
    lmcosi : (Ncoef, 4) float array [l, m, cos_coeff, sin_coeff]
             Ncoef = (lmax+1)(lmax+2)/2, ordered l=0..lmax, m=0..l.
    """
    clm   = np.asarray(clm, dtype=complex).ravel()
    N     = len(clm)
    lmax  = int(round(np.sqrt(N) - 1))
    ncoef = (lmax + 1) * (lmax + 2) // 2
    
    l_out = np.zeros(ncoef, dtype=int)
    m_out = np.zeros(ncoef, dtype=int)
    co    = np.zeros(ncoef, dtype=float)
    si    = np.zeros(ncoef, dtype=float)
    
    idx_out = 0
    for l in range(0, lmax + 1):
        for m in range(0, l + 1):
            l_out[idx_out] = l
            m_out[idx_out] = m
            
            c_pos = clm[lm_to_index(l,  m)]
            c_neg = clm[lm_to_index(l, -m)] if m > 0 else 0.0
            
            if m == 0:
                co[idx_out] = c_pos.real
                si[idx_out] = 0.0
            else:
                phase         = (-1) ** m
                co[idx_out]   = (phase * c_pos + c_neg).real         / np.sqrt(2)
                si[idx_out]   = (1j * (phase * c_pos - c_neg)).real  / np.sqrt(2)
            
            idx_out += 1
    
    return np.column_stack([l_out.astype(float), m_out.astype(float), co, si])






def _clm2xyz(lmcosi, degres):
    """
    Evaluate a real SpH expansion on a regular (lat, lon) grid.
    
    Parameters
    ----------
    lmcosi : (Ncoef, 4) array [l, m, cos_coeff, sin_coeff]
    degres : pixel size in degrees
    
    Returns
    -------
    r   : (nlat, nlon) field values
    lon : (nlon,) longitude in degrees [0, 360]
    lat : (nlat,) latitude  in degrees [-90, 90]
    """
    lmcosi = np.asarray(lmcosi, dtype=float)
    lmin   = int(lmcosi[0, 0])
    L      = int(lmcosi[-1, 0])
    
    nlon = int(np.ceil(360.0 / degres + 1))
    nlat = int(np.ceil(180.0 / degres + 1))
    
    phi       = np.linspace(0, 2 * np.pi, nlon)
    theta     = np.linspace(0, np.pi,     nlat)
    r         = np.zeros((nlat, nlon))
    cos_theta = np.cos(theta)
    
    for l in range(lmin, L + 1):
        ms    = np.arange(0, l + 1)
        clm_l = _schmidt_legendre(l, cos_theta)              # (l+1, nlat)
        mphi  = np.outer(ms, phi)                            # (l+1, nlon)
        
        mask  = lmcosi[:, 0].astype(int) == l
        clm   = lmcosi[mask, 2].reshape(-1, 1)
        slm   = lmcosi[mask, 3].reshape(-1, 1)
        
        r    += clm_l.T @ (clm * np.cos(mphi) + slm * np.sin(mphi))
    
    lon = phi   * 180.0 / np.pi
    lat = 90.0  - theta * 180.0 / np.pi
    return r, lon, lat



def _schmidt_legendre(l, x):
    """
    Schmidt semi-normalized associated Legendre functions scaled by sqrt(2l+1),
    for m = 0..l evaluated at x = cos(theta).
    Returns shape (l+1, len(x)).
    """
    x    = np.asarray(x)
    ms   = np.arange(0, l + 1)
    P    = np.stack([lpmv(m, l, x) for m in ms], axis=0)

    # remove the Condon-Shortley phase, which is automatically includedin scipy.lpmv, because it was added above already
    cs_phase = ((-1) ** ms)[:, None]

    norm = np.sqrt((2.0 - (ms == 0)) * factorial(l - ms) / factorial(l + ms))
    Pbar = norm[:, None] * P * cs_phase * np.sqrt(2 * l + 1)
    return Pbar





# =============================================================================
# |  Old implementation, direct translation from matlab                       |
# =============================================================================

def build_pixel_conversion(lmax, res):
    """
    Build (or return cached) the SpH->pixel transform matrix U,
    shape (Npix, N), where N = (lmax+1)**2, columns in PTA ordering.
    """
    
    cache_folder = cartography_path + "PIXELCONVERSION_CACHE/"
    # check if there the cache folder exists
    if not os.path.exists(cache_folder):
        os.makedirs(cache_folder)
        
    N = (lmax + 1) ** 2
    e0 = np.zeros(N, dtype=complex)
    e0[0] = 1.0
    col0, RA, DECL, dOmg = makemap(e0, res, cmplx=True)
    
    
    # check if the folder already contains the conversion matrix
    if os.path.exists(cache_folder + "LSSM_{}_{}.txt".format(lmax, res)):
        print("\t load existing conversion matrix")
        U = np.loadtxt(cache_folder + "LSSM_{}_{}.txt".format(lmax, res), dtype=np.complex_)       
        pc = {'lmax': lmax, 'res': res, 'U': U, 'RA': RA, 'DECL': DECL, 'dOmg': dOmg}
        
    
    else:
        Npix = col0.shape[0]
        U = np.zeros((Npix, N), dtype=complex)  # real, spans cos AND sin

        U[:, 0]  = col0
        for ii in range(1, N):
            e            = np.zeros(N, dtype=complex)
            e[ii]        = 1.0
            col, _, _, _ = makemap(e, res, cmplx=True)        # <-- use complex variant
            U[:, ii]     = col

        pc = {'lmax': lmax, 'res': res, 'U': U, 'RA': RA, 'DECL': DECL, 'dOmg': dOmg}
        np.savetxt(cache_folder + "LSSM_{}_{}.txt".format(lmax, res), U)

    return pc









def makemap(clm, deg=1.0, cmplx = False):
    """
    Convert complex SpH coefficient vector to a pixel map.
    
    Parameters
    ----------
    clm : 1D complex array, length (lmax+1)**2, PTA ordering.
    deg : pixel resolution in degrees.
    
    Returns
    -------
    map_ : 1D array of pixel values (real)
    RA   : 1D array of right ascensions (hours)
    DECL : 1D array of declinations (degrees)
    dOmg : 1D array of solid-angle elements
    """
    clm_vec       = np.asarray(clm).ravel()
    
    # real part map
    lmcosi        = _clm2clmreal(clm_vec)
    Ymn, lon, lat = _clm2xyz(lmcosi, deg)
    
    # imaginary part map
    lmcosi_i        = _clm2clmreal(-1j * clm_vec)
    Ymn_i, _, _     = _clm2xyz(lmcosi_i, deg)
    Ymn             = Ymn.astype(complex) + 1j * Ymn_i
    
    Nlon = len(lon) - 1
    Nlat = len(lat)
    
    Ymn = Ymn[:, :-1] / np.sqrt(4 * np.pi)
    lon = lon[:-1]
    lat = lat[::-1]
    Ymn = np.flipud(Ymn)
    
    ra   = np.outer(np.ones(Nlat), lon / 15.0)
    decl = np.outer(lat, np.ones(Nlon))
    
    dOmega = np.cos(decl / 180.0 * np.pi) * (2 * np.pi**2 / (Nlon * (Nlat - 1)))
    
    RA   = ra.reshape(Nlon * Nlat, order='F')
    DECL = decl.reshape(Nlon * Nlat, order='F')
    
    if cmplx is True:
        map_ = Ymn.reshape(Nlon * Nlat, order='F')
    elif cmplx is False:
        map_ = np.real( Ymn.reshape(Nlon * Nlat, order='F') )
    else:
        print("Error")
    dOmg = dOmega.reshape(Nlon * Nlat, order='F')
    
    return map_, RA, DECL, dOmg






def diag_pixel(M, res):
    """
    Compute diagonal of U @ M @ U^H in pixel space.
    Returns (d, RA, DECL, dOmg, U).
    """
    lmax    = int(round(np.sqrt(M.shape[0]) - 1))
    pc      = build_pixel_conversion(lmax, res)
    U       = pc['U']
    d       = np.einsum('ij,jk,ik->i', U, M, np.conj(U) )
    return d, pc['RA'], pc['DECL'], pc['dOmg'], U




def get_sigma_map(covarM, res):
    """
    Pixel-domain standard deviation map.
    Returns (sigmaPix, RA, DECL, dOmg, U).
    """
    d, RA, DECL, dOmg, U = diag_pixel(covarM, res)
    sigmaPix = np.sqrt(np.real(d))
    return sigmaPix, RA, DECL, dOmg, U
# =============================================================================








# =============================================================================
# =============================================================================
#   new implementation, faster conversion                                    
# =============================================================================
# =============================================================================



def build_pixel_conversion_fast(lmax, res):
    """
    Build (or return cached) the SpH->pixel transform matrix U,
    shape (Npix, N), where N = (lmax+1)**2, columns in PTA ordering.

    The transform is constructed directly from the spherical-harmonic
    basis rather than calling makemap() once for every basis vector.
    """
    cache_folder = cartography_path + "PIXELCONVERSION_CACHE/"
    os.makedirs(cache_folder, exist_ok=True)
    
    cache_file = cache_folder + "LSSM_{}_{}.npz".format(lmax, res)
    
    if os.path.exists(cache_file):
        print("\t load existing conversion matrix")
        pc = np.load(cache_file)
        
    else:
        # -- create grid on the sphere ----------------------------------------
        
        # drop duplicate longitude at 360 degrees.
        nlon_full = int(np.ceil(360.0 / res + 1))
        phi = np.linspace(0.0, 2.0 * np.pi, nlon_full)[:-1]
        nlon = len(phi)
        
        nlat = int(np.ceil(180.0 / res + 1))
        theta = np.linspace(0.0, np.pi, nlat)

        # transform to radian
        lon = phi * 180.0 / np.pi
        lat = 90.0 - theta * 180.0 / np.pi
        lat = lat[::-1]
        
        # transform to RA [h], DEC [deg]
        RA = np.outer(np.ones(nlat), lon / 15.0).reshape(nlon * nlat, order='F')
        DEC = np.outer(lat, np.ones(nlon)).reshape(nlon * nlat, order='F')
        
        # number of pixels and anular resolution of pixel
        npix = nlat * nlon
        dOmega = np.cos(lon) * (2*np.pi**2 / (nlon*(nlat - 1)))
        # ---------------------------------------------------------------------
   

        # -- calculate conversion matrix U -------------------------------------
        U = np.zeros((npix, (lmax + 1) ** 2), dtype=complex)
        cos_theta = np.cos(theta)
        
        

        # evaluate Legendre polynomials once per l
        # construct all m-columns from that 
        
        for l in range(lmax + 1):
            # flip the latitude axis to match RA orientation
            Pbar = _schmidt_legendre(l, cos_theta)[:, ::-1]      # Shape: (l+1, nlat)

            
            # .. m = 0 ........................................................
            n_m0 = 1/np.sqrt(4*np.pi)
            U[:, lm_to_index(l, 0)] = np.broadcast_to(n_m0 * Pbar[0, :, None], (nlat, nlon)).reshape(npix, order='F')
            
            if l == 0:
                continue
            
            
            # .. m != 0 .......................................................
            ms = np.arange(1, l + 1)
            n_m = 1.0 / np.sqrt(8.0 * np.pi)
            # get indices for the pos and neg m
            
            
            
            # angular part of the matrix entry
            ang = Pbar[1:, :, None] * np.exp(1j * ms[:, None, None] * phi[None, None, :]) # Shape: (m, nlat, nlon)
            ang_pix = ang.transpose(1,2,0).reshape(npix, l, order='F')
            
            
            # . positive m values .......... 
            pos_idx = np.array([lm_to_index(l, m) for m in ms])
            U[:, pos_idx] = ( n_m * ((-1.0)**ms) * ang_pix )[None, :]
            
            # . negative m values .......... 
            neg_idx = np.array([lm_to_index(l, -m) for m in ms])
            U[:, neg_idx] = ( n_m * np.conj(ang_pix) )[None, :]
            
            
            
        pc = {'lmax': lmax,
              'res': res,
              'U': U,
              'RA': RA,
              'DEC': DEC,
              'dOmg': dOmega.reshape(-1, order='F')
              }
        
        np.savez(cache_file, U=U, RA=RA, DEC=DEC, dOmg=dOmega, lmax=lmax, res=res)

    return pc

