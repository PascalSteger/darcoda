#!/usr/bin/env ipython3

##
# @file
# all functions working on spherical coordinates

# (c) GPL v3 2015 ETHZ Pascal Steger, pascal@steger.aero

import pdb
import numpy as np

from scipy.interpolate import splrep, splint
from pylab import *
import gi_helper as gh
import gi_int as gi

def rhodm_hernquist(r, rho0, r_DM, alpha_DM, beta_DM, gamma_DM):
    return rho0*(r/r_DM)**(-gamma_DM)*\
      (1+(r/r_DM)**alpha_DM)**((gamma_DM-beta_DM)/alpha_DM)
## \fn rhodm_hernquist(r,rho0,r_DM,alpha_DM,beta_DM,gamma_DM)
# return hernquist DM density from generalized Hernquist profile
# @param r radius in [pc]
# @param rho0 central 3D density in [Munit/pc^3]
# @param r_DM scale radius of dark matter, [pc]
# @param alpha_DM [1]
# @param beta_DM [1]
# @param gamma_DM [1] central density slope

def nr(r0, dlr, pop, gp):
    # extend asymptotes to 0, and high radius
    rnu = np.hstack([r0[0]/gp.rinfty, r0, gp.rinfty*r0[-1]])
    # up: common radii r0, but different scale radius for each pop
    logrnu = np.log(rnu/gp.dat.rhalf[pop])
    dlrnu = -1.*dlr
    # use linear spline interpolation for dn/dlogr (logr)
    #pdb.set_trace()
    spline_n = splrep(logrnu, dlrnu, k=1)
    # output to evaluate spline at any points in between
    return spline_n
## \fn nr(r0, dlr, pop, gp)
# calculate n(r) at any given radius, as linear interpolation with two asymptotes
# @param r0 radii [pc]
# @param dlr d log rho/ d log r: asymptote at 0, n(r) for all bins, asymptote at infinity
# @param pop int for population (both, 1, 2, ...)
# @param gp global parameters

def rho(r0, rhodmpar, pop, gp):
    gh.sanitize_vector(rhodmpar, gp.nrho, 0, 1e30, gp.debug)
    vec = 1.*rhodmpar # make a new copy so we do not overwrite rhodmpar
    rho_at_rhalf = vec[0]
    vec = vec[1:]
    # get spline representation on gp.xepol, where rhodmpar are defined on
    spline_n = nr(gp.xepol, vec, pop, gp)
    # and apply it to these radii, which may be anything in between
    # self.Xscale is determined in 2D, not 3D
    # use gp.dat.rhalf[0]
    rs =  np.log(r0/gp.dat.rhalf[pop]) # have to integrate in d log(r)
    lnrright = []
    lnrleft = []
    if np.rank(rs) == 0:
        if rs>0:
            lnrright.append(rs)
        else:
            lnrleft.append(rs)
    else:
        lnrright = rs[(rs>=0.)]
        lnrleft  = rs[(rs<0.)]
        lnrleft  = lnrleft[::-1] # inverse order

    lnrhoright = []
    for i in np.arange(0, len(lnrright)):
        lnrhoright.append(np.log(rho_at_rhalf) + \
                           splint(0., lnrright[i], spline_n))
                           # integration along dlog(r) instead of dr

    lnrholeft = []
    for i in np.arange(0, len(lnrleft)):
        lnrholeft.append(np.log(rho_at_rhalf) + \
                          splint(0., lnrleft[i], spline_n))

    tmp = np.exp(np.hstack([lnrholeft[::-1], lnrhoright])) # still defined on ln(r)
    gh.checkpositive(tmp, 'rho()')
    return tmp
## \fn rho(r0, rhodmpar, pop, gp)
# calculate density, from interpolated n(r) = -d ln(rho(r))/ d ln(r)
# using interpolation to left and right of r=r_{*, 1/2}
# @param rhodmpar: rho(rstarhalf), asymptote nr 0, nr(xepol), asymptote nr infty
# @param r0 radii to calculate density for, in physical units (pc)
# @param pop int for population, 0 all or DM, 1, 2, ...
# @param gp global parameters


def betastar2beta(betastar):
    # beta^* = \frac{sig_r^2-sig_t^2}{sig_r^2+sig_t^2}
    # betastar is in [-1,1]
    # and symmetric
    return 2.*betastar/(1.+betastar)
## \fn betastar2beta(betastar)
# map beta* to beta
# @param betastar float array

def beta2betastar(beta):
    return beta/(2.-beta)
## \fn beta2betastar(beta)
# get back betastar from a beta
# @param beta [1]

def betastar_4(r0, params, gp):
    gh.sanitize_vector(params, gp.nbeta, -1, max(gp.xepol), gp.debug)
    a0 = params[0]
    a1 = params[1]
    alpha = params[2]
    s0 = np.log(r0/np.exp(params[3])) # r_s=params[3] is given in log
    betatmp = (a0-a1)/(1+np.exp(alpha*s0))+a1
    return betatmp
## \fn betastar_4(r0, params, gp)
# calculate betastar from sigmoid with 4 parameters, using exp directly, with explicit meaning
# @param r0 radii [pc]
# @param params 4 parameters: asymptote of beta at r->0, asymptote of beta at r->infty, speed of change, scale radius at which change takes place
# @param gp global parameters

def betastar(r0, params, gp):
    bs = betastar_4(r0, params, gp)
    for k in range(len(r0)):
        if bs[k] < min(gp.minbetastar_0, gp.minbetastar_inf):
            bs[k] = min(gp.minbetastar_0, gp.minbetastar_inf)
        if bs[k] > max(gp.maxbetastar_0, gp.maxbetastar_inf):
            bs[k] = max(gp.maxbetastar_0, gp.maxbetastar_inf)
    return bs
## \fn betastar(r0, params, gp)
# calculate betastar from 4 parameters
# @param r0 radii [pc]
# @param params 4 parameters
# @param gp global parameters

def beta(r0, params, gp):
    bstar = betastar(r0, params, gp)
    betatmp = betastar2beta(bstar)
    return betatmp, bstar
## \fn beta(r0, params, gp)
# beta and beta* from beta parameter array
# @param r0 radii [pc]
# @param params array, see gi_class_cube
# @param gp global parameters

def calculate_rho(r, M):
    r0 = np.hstack([0,r])                         # [lunit]
    M0 = np.hstack([0,M])
    deltaM   = M0[1:]-M0[:-1]                     # [Munit]
    gh.checkpositive(deltaM,  'unphysical negative mass increment encountered')
    deltavol = 4./3.*np.pi*(r0[1:]**3-r0[:-1]**3) # [lunit^3]
    rho     = deltaM / deltavol                  # [Munit/lunit^3]
    gh.checkpositive(rho, 'rho in calculate_rho')
    return rho                                   # [Munit/lunit^3]
## \fn calculate_rho(r, M)
# take enclosed mass, radii, and compute 3D density in shells
# @param r radius in [pc]
# @param M enclosed mass profile, [Munit]

def calculate_surfdens(r, M):
    r0 = np.hstack([0,r])                         # [lunit]
    M0 = np.hstack([0,M])                         # [Munit]
    deltaM = M0[1:]-M0[:-1]                       # [Munit]
    gh.checkpositive(deltaM, 'unphysical negative mass increment encountered')
    deltavol = np.pi*(r0[1:]**2 - r0[:-1]**2)        # [lunit^2]
    Sig = deltaM/deltavol                           # [Munit/lunit^2]
    gh.checkpositive(Sig, 'Sig in calculate_surfdens')
    return Sig                                      # [Munit/lunit^2]
## \fn calculate_surfdens(r, M)
# take mass(<r) in bins, calc mass in annuli, get surface density
# @param r radius in [pc]
# @param M 3D mass profile

def sig_kap_zet(r0, rhodmpar, lbaryonpar, MtoL, nupar, betapar, pop, gp):
    siglos2, kaplos4, zetaa, zetab = gi.ant_sigkaplos(r0, rhodmpar, lbaryonpar, MtoL, nupar, betapar, pop, gp)
    siglos  = np.sqrt(siglos2)           # [km/s]

    if gp.usekappa:
        kaplos4 = kaplos4surf/Sig
        # takes [Munit/pc^2 (km/s)^2], gives back [(km/s)^2]

        kaplos  = kaplos4/(siglos2**2)
        # - 3.0 # subtract 3.0 for Gaussian distribution in Fisher version.
    else:
        kaplos = 3.*np.ones(len(siglos))

    return siglos, kaplos, zetaa, zetab  # [km/s], [1]
## \fn sig_kap_zet(r0, rhodmpar, lbaryonpar, MtoL, nupar, betapar, pop, gp)
# General function to calculate sig_los
# with analytic integral over fitting polynomial'
# @param r0 radius [pc]
# @param rhodmpar [1]
# @param lbaryonpar [1]
# @param MtoL mass-to-light ratio [Msun/Lsun]
# @param nupar 1D array 3D tracer density parameters [1]
# @param betapar 1D array 3D velocity anisotropy parameters [1]
# @param pop int population to take halflight radius from
# @param gp global parameters


if __name__=="__main__":
    # check whether backward-calculation of parameters is done correctly
    import gi_params
    gp = gi_params.Params()
    annu = rhodm_hernquist(gp.xipol, 1.e0, max(gp.xipol)/2., 3, 2, 1)
