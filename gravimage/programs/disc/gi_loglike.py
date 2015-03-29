#!/usr/bin/env ipython3

## @file
# define log likelihood function to be called by Multinest
# disc version

# (c) GPL v3 2015 ETHZ, Pascal Steger pascal@steger.aero

import numpy as np
import pdb
import gi_helper as gh
from gi_class_profiles import Profiles
from gi_priors import check_tilt
from gi_chi import calc_chi2
import gi_physics as phys
#from pylab import *
#ion()

def geom_loglike(cube, ndim, nparams, gp):
    tmp_profs = Profiles(gp.pops, gp.nipol)
    off = 0
    offstep = 1
    norm = cube[off]
    off += offstep

    offstep = gp.nrho
    rhodmpar = np.array(cube[off:off+offstep])  #SS cube[1:1+nrho]
    tmp_rho = phys.rho(gp.xepol, rhodmpar, 0, gp)
    tmp_profs.set_prof('rho', tmp_rho[gp.nexp:-gp.nexp], 0, gp)
    off += offstep

    offstep = gp.nrho
    lbaryonpar = np.array(cube[off:off+offstep]) #SS cube[1+nrho:1+2*nrho]
    tmp_rhostar = phys.rho(gp.xepol, lbaryonpar, 0, gp)[gp.nexp:-gp.nexp]
    tmp_profs.set_prof('nu', tmp_rhostar, 0, gp) # [Munit/pc^3]
    Sigstar = phys.nu_SUM_Sig(gp.dat.binmin, gp.dat.binmax, tmp_rhostar) # [Munit/pc^2]
    tmp_profs.set_prof('Sig', Sigstar, 0, gp)
    off += offstep

    MtoL = cube[off]  #SS cube[1+2*nrho]
    off += 1

    for pop in np.arange(1, gp.pops+1):
        offstep = gp.nrho
        nupar = np.array(cube[off:off+offstep])  #SS 1 cube[2+2*nrho:2+3*nrho]
        tmp_nu = phys.rho(gp.xepol, nupar, pop, gp)[gp.nexp:-gp.nexp]
        tmp_profs.set_prof('nu', tmp_nu, pop, gp) # [Munit/pc^3]
        tmp_Sig = phys.nu_SUM_Sig(gp.dat.binmin, gp.dat.binmax, tmp_nu) # [Munit/pc^2]
        tmp_profs.set_prof('Sig', tmp_Sig, pop, gp)
        off += offstep

        if gp.checksig:
            pdb.set_trace()

        offstep = gp.nbeta
        if gp.chi2_nu_converged:
            tiltpar = np.array(cube[off:off+offstep])#SS cube[2+3*nrho:..+nbeta]
            tmp_tilt = phys.tilt(gp.xipol, tiltpar, gp)
            if check_tilt(tmp_tilt, gp):
                gh.LOG(1, 'tilt error')
                tmp_profs.chi2 = gh.err(2., gp)
                return tmp_profs
            tmp_profs.set_prof('tilt', tmp_tilt, pop, gp)
            sig = phys.sigz(gp.xepol, rhodmpar, lbaryonpar, MtoL, nupar, norm, tiltpar, pop, gp)
            tmp_profs.set_prof('sig', sig, pop, gp)
            # tmp_profs.set_prof('kap', kap, pop, gp)
        off += offstep # add also in case Sig has not yet converged
        # to get the right variables

    if off != gp.ndim:
        gh.LOG(1,'wrong subscripts in gi_class_cube')
        raise Exception('wrong subscripts in gi_class_cube')

    # determine log likelihood
    chi2 = calc_chi2(tmp_profs, gp)
    gh.LOG(1, '   log L = ', -chi2/2.)
    tmp_profs.chi2 = chi2

    return tmp_profs   # from   likelihood L = exp(-\chi^2/2), want log of that
## \fn geom_loglike(cube, ndim, nparams, gp)
# define log likelihood function to be called by Multinest and plot_profiles
# disc version
# @param cube parameter cube as defined by gi_class_cube, in physical space
# @param ndim number of dimensions, needed as argument by multinest
# @param nparams number of parameters, needed as argument by multinest
# @param gp global parameters
