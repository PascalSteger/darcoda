#!/usr/bin/env ipython-python3.2
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
# all functions working on spherical coordinates

import pdb
import numpy as np
import gl_params as gp
import gl_funs as gfun
import gl_plot as gpl
import gl_helper as gh
from gl_int import *
from gl_project import rho_INT_Rho
from gl_analytic import *

def rhodm_hernquist(r,rho0,r_DM,alpha_DM,beta_DM,gamma_DM):
    'return hernquist DM density from generalized Hernquist profile'
    return rho0*(r/r_DM)**(-gamma_DM)*\
      (1+(r/r_DM)**alpha_DM)**((gamma_DM-beta_DM)/alpha_DM)

def densdefault(denspars):
    'return density from denspars, assuming default gp.xipol as radius'
    return dens(gp.xipol,denspars)        # [TODO]

def dens(xipol, denspars):
    # take denspars for polynomial coefficients, calculate polynomial, 
    # give back density,
    # [pc], [1] or [msun/pc^3]'
    if gp.checkint:
        return gp.ipol.densdat          # [Msun/pc^3]
    if not gp.poly:                     # for sure after init
        if gp.denslog:
            return 10.**denspars        # [Msun/pc^3]
        else:
            return denspars             # [Msun/pc^3]

    scale = gp.scaledens*max(xipol)     # [pc]

    tmp = np.zeros(len(xipol))
    for i in range(0,len(denspars)):
        tmp += denspars[i]*((scale-xipol)/scale)**i # [log10(msun/pc^3)]

    # gpl.plot(gp.ipol.densr,gp.ipol.densdat); 
    # gpl.plot(gp.xipol,10.**tmp); gpl.yscale('log')
    
    dout = np.power(10.,tmp)   # [Msun/pc^3] # always in denslog case!
    gh.checknan(dout)
    return dout                               # [msun/pc^3]

def delta(dpars):
    '''calculate cumulative sum for representation of delta'''
    return np.cumsum(dpars)

def invdelta(delta):
    '''calculate delta parameters corresponding to a given delta profile'''
    return np.hstack([delta[0], np.diff(delta)])


def calculate_dens(r, M):
    'take enclosed mass, radii, and compute 3D density in shells'
    r0 = np.hstack([0,r])                         # [lunit]
    M0 = np.hstack([0,M])

    deltaM   = M0[1:]-M0[:-1]                     # [munit]
    if min(deltaM)<0.:
        print('unphysical negative mass increment encountered')
        exit(1)

    deltavol = 4./3.*np.pi*(r0[1:]**3-r0[:-1]**3) # [lunit^3]
    dens     = deltaM / deltavol                  # [munit/lunit^3]
    gh.checknan(dens)
    return dens                                   # [munit/lunit^3]

def calculate_surfdens(r, M):
    'take mass(<r) in bins, calc mass in annuli, get surface density'
    r0 = np.hstack([0,r])                         # [lunit]
    M0 = np.hstack([0,M])                         # [munit]

    deltaM = M0[1:]-M0[:-1]                       # [munit]
    if min(deltaM)<0.:
        print('unphysical negative mass increment encountered')
        exit(1)
    
    deltavol = np.pi*(r0[1:]**2 - r0[:-1]**2)        # [lunit^2]
    dens = deltaM/deltavol                           # [munit/lunit^2]
    gh.checknan(dens)
    return dens                                      # [munit/lunit^2]

def nu(nupars):
    'General function to describe the density profile., [(log10) Msun/pc^3]'
    if gp.nulog:
        nuout = np.power(10.0, nupars) # [Msun/pc^3]
    else:
        nuout = nupars[:]       # [Musn/pc^3]

    gh.checknan(nuout)
    if gp.geom == 'disc':
        print('nu in spherical phys taken! bug!')
        pdb.set_trace()
        nuout = nuout/max(nuout)
    return nuout                # [Msun/pc^3]


def sig_kap_los(dens_r, nu_r, beta_r):
    # General function to calculate \sigma_{los} 
    # with analytic integral over fitting polynomial'
    # [munit, 3D], [munis/pc^3], [1]
    r0 = gp.xipol[:]

    # TODO: fix spline warnings
    surfden = rho_INT_Rho(r0, nu_r) # or INTIPOL
    # takes [pc], [munit/pc^3], gives back [munit/pc^2]
    # gpl.start()
    # pdb.set_trace()
    # gpl.plot(r0, surfden,'red')
    # gpl.plot(r0, Sigma_anf(r0))

    # Calculate density and M force:
    intbeta   = ant_intbeta(r0, beta_r) # [1]
    # TODO:pdb.set_trace() get splines correctly
    siglos2surf,kaplos4surf = ant_sigkaplos2surf(r0,beta_r,intbeta,dens_r,nu_r)
    # takes [pc], [1*pc], [munit], [munit/pc^3], gives back [(km/s)^2], [1]

    siglos2 = siglos2surf/surfden   # [(km/s)^2]
    siglos  = np.sqrt(siglos2)      # [km/s]

    # gpl.start(); pdb.set_trace()
    # gpl.plot(r0,sig_los_anf(r0),lw=2)
    
    kaplos4     = kaplos4surf/surfden
    # takes [munit/pc^2 (km/s)^2], gives back [(km/s)^2]
    
    kaplos = kaplos4/(siglos2**2)
    # - 3.0 # subtract 3.0 for Gaussian distribution in Fisher version.

    return siglos,kaplos                                 # [km/s], [1]
