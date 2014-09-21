#!/usr/bin/env ipython3

## @file
# define log likelihood function to be called by pyMultinest
# disc version

import numpy as np
import pdb
import gl_helper as gh
from gl_class_profiles import Profiles
from gl_priors import check_bprior, check_tilt
from gl_chi import calc_chi2
import gl_physics as phys

    
def geom_loglike(cube, ndim, nparams, gp):
    tmp_profs = Profiles(gp.pops, gp.nipol)
    off = 0
    norm = cube[off]
    off += 1

    rhopar = np.array(cube[off:off+gp.nrho])
    tmp_rho = phys.rho(gp.xepol, rhopar, 0, gp)
    tmp_profs.set_prof('rho', tmp_rho[3:-3], 0, gp)
    off += gp.nrho
    # TODO: M
    # tmp_M = glp.rho_SUM_Mr(gp.xepol, tmp_rho)
    # tmp_profs.set_prof('M', tmp_M[:gp.nipol], 0, gp)


    for pop in np.arange(gp.pops)+1:
        nupar = np.array(cube[off:off+gp.nepol])
        tmp_profs.set_prof('nu', nupar, pop, gp) # [Munit/pc^3]
        off += gp.nepol

        tiltpar = np.array(cube[off:off+gp.nbeta])
        tmp_tilt = phys.tilt(gp.xipol, tiltpar, gp)
        if check_tilt(tmp_tilt, gp):
            print('tilt error')
            tmp_profs.chi2 = gh.err(2., gp)
            return tmp_profs
        tmp_profs.set_prof('tilt', tmp_tilt, pop, gp)
        off += gp.nbeta

        #try:
        sig = phys.sigz(gp.xepol, rhopar, nupar, norm, tiltpar, pop, gp)
        #except Exception as detail:
        #    print('sig kap exception')
        #    tmp_profs.chi2 = gh.err(3., gp)
        #    return tmp_profs
        tmp_profs.set_prof('sig', sig, pop, gp)
        # tmp_profs.set_prof('kap', kap, pop, gp)
    
    # determine log likelihood (*not* reduced chi2)
    chi2 = calc_chi2(tmp_profs, gp)
    print('found log likelihood = ', -chi2/2.)
    tmp_profs.chi2 = chi2
    return tmp_profs   # from   likelihood L = exp(-\chi^2/2), want log of that
## \fn geom_loglike(cube, ndim, nparams, gp)
# define log likelihood function to be called by pyMultinest and plot_profiles
# disc version
# @param cube parameter cube as defined by gl_class_cube, in physical space
# @param ndim number of dimensions, needed as argument by multinest
# @param nparams number of parameters, needed as argument by multinest
# @param gp global parameters
