#!/usr/bin/env python3

##
# @file
# all functions working on spherical coordinates

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import pdb
import numpy as np
import gl_params as gp
import gl_chi as gc
import gl_helper as gh
from gl_int import *
from gl_project import rho_INT_Rho, rho_param_INT_Rho
from gl_analytic import *

def rhodm_hernquist(r, rho0, r_DM, alpha_DM, beta_DM, gamma_DM):
    return rho0*(r/r_DM)**(-gamma_DM)*\
      (1+(r/r_DM)**alpha_DM)**((gamma_DM-beta_DM)/alpha_DM)
## \fn rhodm_hernquist(r,rho0,r_DM,alpha_DM,beta_DM,gamma_DM)
# return hernquist DM density from generalized Hernquist profile
# @param r radius in [pc]
# @param rho0 central 3D density in [Msun/pc^3]
# @param r_DM scale radius of dark matter, [pc]
# @param alpha_DM [1]
# @param beta_DM [1]
# @param gamma_DM [1] central density slope

      
def nr(r0, dlogrhodlogr):
    # extend asymptotes to 0, and high radius
    r0 = np.hstack([r0[0]/1e4, r0[0]/2., r0, gp.rinfty*max(r0)])
    logr0 = np.log(r0/gp.rstarhalf)
    dlogrhodlogr = np.hstack([dlogrhodlogr[0], dlogrhodlogr]) # dlogrhodlogr[-1]])
    dlogrhodlogr *= -1.
    
    # use linear spline interpolation in r
    spline_n = splrep(logr0, dlogrhodlogr, k=1)
    
    # evaluate spline at any points in between
    return spline_n #, splev(r0, spline_n)
## \fn nr(r0, dlogrhodlogr, defrext)
# calculate n(r) at any given radius, as linear interpolation with two asymptotes
# @param dlogrhodlogr : asymptote at 0, n(r) for all bins, asymptote at infinity
# @parma defr0  corresponding radii
# @param r0     new radii to interpolate to


def rho(r0, arr):
    rhoathalf = arr[0]
    arr = arr[1:]
    spline_n = nr(gp.xepol, arr) # fix on gp.xepol, as this is where definition is on

    rs =  np.log(r0/gp.rstarhalf) # have to integrate in d log(r)

    # check assumption r_min < rstarhalf < r_max
    if min(rs)>0. or max(rs)<0.:
        print('assumption r_min < rstarhalf < r_max violated')
        pdb.set_trace()
    logrright = rs[(rs>=0.)]
    logrleft  = rs[(rs<0.)]
    logrleft  = logrleft[::-1]
    logrhoright = []
    for i in np.arange(0, len(logrright)):
        logrhoright.append(np.log(rhoathalf) + \
                           splint(0., logrright[i], spline_n))
                           # TODO: integration along dlog(r) instead of dr?!

    logrholeft = []
    for i in np.arange(0, len(logrleft)):
        logrholeft.append(np.log(rhoathalf) + \
                          splint(0., logrleft[i], spline_n))

    tmp = np.exp(np.hstack([logrholeft[::-1], logrhoright])) # still defined on log(r)
    gh.checkpositive(tmp, 'rho()')
    return tmp
## \fn rho(r0, arr)
# calculate density, from interpolated n(r) = -log(rho(r))
# using interpolation to left and right of r=r_{*, 1/2}
# @param arr: rho(rstarhalf), asymptote_0, nr(gp.xipol), asymptote_infty
# @param r0: radii to calculate density for, in physical units (pc)


# def rho_plug_formula(r0, arr, rarr, rstarhalf):
#     rhohalf = arr[0]
#     narr    = arr[1:]
#     return rhohalf * (r0/rstarhalf)**(-nr(r0, narr, rarr))
## \fn rho_plug_formula(r0, arr, rarr, rstarhalf)
# scale [0,1] to [0, maxrange], monotonically decreasing for all radii
# thus setting n = - arr = - d ln(rho)/ d ln(r)
# rho(r) = rho(rhalf) * (r/rhalf)^{-n(r/rhalf)}
# n(r/rhalf) is defined at Nbin points, and two asymptotic values for r\to0,\infty
# and linearly interpolated in between
# @param arr array of floats, [0..1]
# @param maxrange float of maximum
    
# old version:
# logrho = np.zeros(len(arr))
# logrho[0] = arr[0] * np.log(maxrange)
# for i in np.arange(1,len(arr)):
#     logrho[i] = logrho[i-1]-arr[i]*(np.log(r0[i])-np.log(r0[i-1]))
# return np.exp(logrho)


def beta_star_to_beta(betastar):
    # \beta^* = \frac{\sigma_r^2-\sigma_t^2}{\sigma_r^2+\sigma_t^2}
    # betastar is in [-1,1]
    # and symmetric
    return 2.*betastar/(1.+betastar)
## \fn beta_star_to_beta(betastar)
# map beta* to beta


def mapping_beta_star_poly(r0, arr):
    betatmp = 0.
    for i in range(gp.nbeta):
        betatmp += arr[i] * (r0/max(gp.xipol))**i
    # clipping beta* to the range [-1,1]
    # thus not allowing any unphysical beta,
    # but still allowing parameters to go the the max. value
    for i in range(len(r0)):
        if betatmp[i] > 1.:
            betatmp[i] = 1.
        if betatmp[i] <= -1.:
            betatmp[i] = -0.99999999999 # excluding -inf values in beta
    return betatmp
## \fn mapping_beta_star_poly(r0, arr)
# map [0,1] to [-1,1] with a polynomial
# @param arr normalized a_i, s.t. abs(sum(a_i)) = 1


def beta(r0, arr):
    betastar = mapping_beta_star_poly(r0, arr)
    betatmp  = beta_star_to_beta(betastar)
    return betatmp
## \fn beta(r0, arr)
# return beta calculated from parameters a_i of polynomial representation
# @param arr parameters a_i


def invdelta(delta):
    return np.hstack([delta[0], np.diff(delta)])
## \fn invdelta(delta)
# calculate delta parameters corresponding to a given delta profile
# @param delta [1] velocity anisotropy profile

    
def calculate_rho(r, M):
    r0 = np.hstack([0,r])                         # [lunit]
    M0 = np.hstack([0,M])

    deltaM   = M0[1:]-M0[:-1]                     # [munit]
    gh.checkpositive(deltaM,  'unphysical negative mass increment encountered')

    deltavol = 4./3.*np.pi*(r0[1:]**3-r0[:-1]**3) # [lunit^3]
    rho     = deltaM / deltavol                  # [munit/lunit^3]
    gh.checkpositive(rho, 'rho in calculate_rho')
    return rho                                   # [munit/lunit^3]
## \fn calculate_rho(r, M)
# take enclosed mass, radii, and compute 3D density in shells
# @param r radius in [pc]
# @param M enclosed mass profile, [Msun]


def calculate_surfdens(r, M):
    r0 = np.hstack([0,r])                         # [lunit]
    M0 = np.hstack([0,M])                         # [munit]

    deltaM = M0[1:]-M0[:-1]                       # [munit]
    gh.checkpositive(deltaM, 'unphysical negative mass increment encountered')
    
    deltavol = np.pi*(r0[1:]**2 - r0[:-1]**2)        # [lunit^2]
    Rho = deltaM/deltavol                           # [munit/lunit^2]
    gh.checkpositive(Rho, 'Rho in calculate_surfdens')
    return Rho                                      # [munit/lunit^2]
## \fn calculate_surfdens(r, M)
# take mass(<r) in bins, calc mass in annuli, get surface density
# @param r radius in [pc]
# @param M 3D mass profile


def sig_kap_los(r0, pop, rho_param, nu_param, beta_param):
    surfden = rho_param_INT_Rho(r0, rho_param) # total surface density
    siglos2surf,kaplos4surf = ant_sigkaplos2surf(r0, beta_param,\
                                                 rho_param, nu_param)
    # takes [pc], [1*pc], [munit], [munit/pc^3], returns [(km/s)^2], [1]

    siglos2 = siglos2surf[:-3]/surfden   # [(km/s)^2]
    siglos  = np.sqrt(siglos2)           # [km/s]

    if gp.usekappa:
        kaplos4     = kaplos4surf/surfden
        # takes [munit/pc^2 (km/s)^2], gives back [(km/s)^2]
        
        kaplos = kaplos4/(siglos2**2)
        # - 3.0 # subtract 3.0 for Gaussian distribution in Fisher version.
    else:
        kaplos = 3.*np.ones(len(siglos))

    return siglos, kaplos                                 # [km/s], [1]
## \fn sig_kap_los(r0, pop, rho_param, nu_param, beta_param)
# General function to calculate \sigma_{los} 
# with analytic integral over fitting polynomial'
# @param r0 radius [pc]
# @param pop [0 (all), 1, 2]
# @param rho_param [1]
# @param nu_param [1]
# @param beta_param [1] velocity anisotropy
