#!/usr/bin/env ipython3

## @file
# define log likelihood function to be called by pyMultinest
# disc version

import numpy as np
import ipdb
import gl_helper as gh
from gl_class_profiles import Profiles
from gl_priors import check_bprior, check_tilt
from gl_chi import calc_chi2
import gl_physics as phys
from pylab import *
ion()

def geom_loglike(cube, ndim, nparams, gp):
    tmp_profs = Profiles(gp.pops, gp.nipol)
    off = 0
    offstep = 1
    norm = cube[off]
    off += offstep

    offstep = gp.nrho
    rhopar = np.array(cube[off:off+offstep])
    tmp_rho = phys.rho(gp.xepol, rhopar, 0, gp)
    tmp_profs.set_prof('rho', tmp_rho[gp.nexp:-gp.nexp], 0, gp)
    off += offstep

    # TODO: M
    # tmp_M = glp.rho_SUM_Mr(gp.xepol, tmp_rho)
    # tmp_profs.set_prof('M', tmp_M[gp.nexp:-gp.nexp], 0, gp)

    offstep = gp.nrho
    rhostarpar = np.array(cube[off:off+offstep])
    tmp_rhostar = phys.rho(gp.xepol, rhostarpar, 0, gp)[gp.nexp:-gp.nexp]
    tmp_profs.set_prof('nu', tmp_rhostar, 0, gp) # [Munit/pc^3]
    Sigstar = phys.nu_SUM_Sig(gp.dat.binmin, gp.dat.binmax, tmp_rhostar) # [Munit/pc^2]
    tmp_profs.set_prof('Sig', Sigstar, 0, gp)
    off += offstep

    MtoL = cube[off]
    off += 1

    for pop in np.arange(1, gp.pops+1):
        offstep = gp.nrho
        nupar = np.array(cube[off:off+offstep])
        tmp_nu = phys.rho(gp.xepol, nupar, pop, gp)[gp.nexp:-gp.nexp]
        tmp_profs.set_prof('nu', tmp_nu, pop, gp) # [Munit/pc^3]
        tmp_Sig = phys.nu_SUM_Sig(gp.dat.binmin, gp.dat.binmax, tmp_nu) # [Munit/pc^2]
        tmp_profs.set_prof('Sig', tmp_Sig, pop, gp)
        off += offstep

        if gp.checksig:
            ipdb.set_trace()

        offstep = gp.nbeta
        if gp.chi2_nu_converged:
            tiltpar = np.array(cube[off:off+offstep])
            tmp_tilt = phys.tilt(gp.xipol, tiltpar, gp)
            if check_tilt(tmp_tilt, gp):
                gh.LOG(1, 'tilt error')
                tmp_profs.chi2 = gh.err(2., gp)
                return tmp_profs
            tmp_profs.set_prof('tilt', tmp_tilt, pop, gp)
            sig = phys.sigz(gp.xepol, rhopar, rhostarpar, MtoL, nupar, norm, tiltpar, pop, gp)
            tmp_profs.set_prof('sig', sig, pop, gp)
            # tmp_profs.set_prof('kap', kap, pop, gp)
        off += offstep # add also in case Sig has not yet converged
        # to get the right variables

    if off != gp.ndim:
        gh.LOG(1,'wrong subscripts in gl_class_cube')
        raise Exception('wrong subscripts in gl_class_cube')

    # determine log likelihood
    chi2 = calc_chi2(tmp_profs, gp)
    gh.LOG(1, '   log L = ', -chi2/2.)
    tmp_profs.chi2 = chi2

    return tmp_profs   # from   likelihood L = exp(-\chi^2/2), want log of that
## \fn geom_loglike(cube, ndim, nparams, gp)
# define log likelihood function to be called by pyMultinest and plot_profiles
# disc version
# @param cube parameter cube as defined by gl_class_cube, in physical space
# @param ndim number of dimensions, needed as argument by multinest
# @param nparams number of parameters, needed as argument by multinest
# @param gp global parameters
