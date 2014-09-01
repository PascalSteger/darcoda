#!/usr/bin/env ipython3

## @file
# define log likelihood function to be called by pyMultinest
# spherical version

import numpy as np
import pdb
from pylab import *

import gl_physics as phys

from gl_class_profiles import Profiles
from gl_priors import check_bprior, check_beta
from gl_chi import calc_chi2
import gl_helper as gh
import gl_project as glp


def geom_loglike(cube, ndim, nparams, gp):
    tmp_profs = Profiles(gp.pops, gp.nipol)
    off = 0

    rhopar = np.array(cube[off:off+gp.nrho])
    tmp_profs.set_prof('nr', 1.*rhopar[1+1+gp.nexp:-gp.nexp-1], 0, gp)

    tmp_rho = phys.rho(gp.xepol, rhopar, 0, gp)
    # rhopar hold [rho(rhalf), nr to be used for integration
    # from halflight radius, defined on gp.xepol]
    tmp_profs.set_prof('rho', tmp_rho[gp.nexp:-gp.nexp], 0, gp)

    # (only calculate) M, check
    tmp_M = glp.rho_SUM_Mr(gp.xepol, tmp_rho)
    tmp_profs.set_prof('M', tmp_M[gp.nexp:-gp.nexp], 0, gp)
    off += gp.nrho

    # get profile for rho*
    rhostarpar = np.array(cube[off:off+gp.nrho])
    rhostar = phys.rho(gp.xepol, rhostarpar, 0, gp)
    tmp_profs.set_prof('nu', rhostar[gp.nexp:-gp.nexp], 0, gp)
    off += gp.nrho

    Signu = glp.rho_INTIPOL_Rho(gp.xepol, rhostar, gp) # [Munit/pc^2]
    Sig = gh.linipollog(gp.xepol, Signu, gp.xipol)
    tmp_profs.set_prof('Sig', Sig, 0, gp)

    MtoL = cube[off]
    off += 1

    for pop in np.arange(1, gp.pops+1):  # [1, 2, ..., gp.pops]
        nupar = np.array(cube[off:off+gp.nrho])
        tmp_nrnu = 1.*nupar[1+1+gp.nexp:-gp.nexp-1]
        tmp_profs.set_prof('nrnu', tmp_nrnu, pop, gp)
        tmp_nu = phys.rho(gp.xepol, nupar, pop, gp)
        tmp_profs.set_prof('nu', tmp_nu[gp.nexp:-gp.nexp], pop, gp)
        off += gp.nrho

        betapar = np.array(cube[off:off+gp.nbeta])
        off += gp.nbeta

        tmp_Signu = glp.rho_INTIPOL_Rho(gp.xepol, tmp_nu, gp) # [Munit/pc^2]
        tmp_Sig = gh.linipollog(gp.xepol, tmp_Signu, gp.xipol)
        tmp_profs.set_prof('Sig', tmp_Sig, pop, gp)

        tmp_beta, tmp_betastar = phys.beta(gp.xipol, betapar, gp)

        if check_beta(tmp_beta, gp):
            print('beta error')
            tmp_profs.chi2 = gh.err(1., gp)
            return tmp_profs
        tmp_profs.set_prof('beta', tmp_beta, pop, gp)
        tmp_profs.set_prof('betastar', tmp_betastar, pop, gp)
        
        try:
            sig,kap,zetaa,zetab=phys.sig_kap_zet(gp.xepol, rhopar, rhostarpar, MtoL, nupar, betapar, pop, gp)
        
        except Exception as detail:
            tmp_profs.chi2 = gh.err(2., gp)
            return tmp_profs
        tmp_profs.set_prof('sig', sig, pop, gp)
        tmp_profs.set_prof('kap', kap, pop, gp)
        tmp_profs.set_prof('zetaa', zetaa, pop, gp)
        tmp_profs.set_prof('zetab', zetab, pop, gp)
        

    # determine log likelihood
    chi2 = calc_chi2(tmp_profs, gp)
    print('found log likelihood = ', -chi2/2.)
    tmp_profs.chi2 = chi2

    if gp.checksig:
        import gl_analytic as ga
        rho0check, rho1check = ga.rho_gaia(gp.xepol, gp)
        dummy, beta1check = ga.beta_gaia(gp.xepol, gp)
        betastar1check = phys.beta2betastar(beta1check)
        pdb.set_trace()

    return tmp_profs
## \fn geom_loglike(cube, ndim, nparams, gp)
# define log likelihood function to be called by pyMultinest and plot_profiles
# spherical version
# @param cube parameter cube as defined by gl_class_cube, in physical space already (not [0,1] cube anymore)
# @param ndim number of dimensions, needed as argument by multinest
# @param nparams number of parameters, needed as argument by multinest
# @param gp global parameters
