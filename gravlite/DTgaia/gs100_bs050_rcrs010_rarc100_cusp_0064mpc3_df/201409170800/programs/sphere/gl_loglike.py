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
from gl_helper import LOG
import gl_project as glp


def geom_loglike(cube, ndim, nparams, gp):
    tmp_profs = Profiles(gp.pops, gp.nipol)
    off = 0

    offstep = gp.nrho
    if gp.chi2_Sig_converged:
        rhopar = np.array(cube[off:off+offstep])
        tmp_profs.set_prof('nr', 1.*rhopar[1+1+gp.nexp:-gp.nexp-1], 0, gp)
        tmp_rho = phys.rho(gp.xepol, rhopar, 0, gp)
        # rhopar hold [rho(rhalf), nr to be used for integration
        # from halflight radius, defined on gp.xepol]
        tmp_profs.set_prof('rho', tmp_rho[gp.nexp:-gp.nexp], 0, gp)
        # (only calculate) M, check
        tmp_M = glp.rho_SUM_Mr(gp.xepol, tmp_rho)
        tmp_profs.set_prof('M', tmp_M[gp.nexp:-gp.nexp], 0, gp)
    else:
        rhopar = np.zeros(gp.nrho)
    off += offstep # anyhow, even if Sig not yet converged

    # get profile for rho*
    if gp.investigate == 'obs':
        offstep = gp.nrho
        rhostarpar = np.array(cube[off:off+offstep])
        rhostar = phys.rho(gp.xepol, rhostarpar, 0, gp)
        tmp_profs.set_prof('nu', rhostar[gp.nexp:-gp.nexp], 0, gp)
        off += offstep

        Signu = glp.rho_INTIPOL_Sig(gp.xepol, rhostar, gp) # [Munit/pc^2]
        Sig = gh.linipollog(gp.xepol, Signu, gp.xipol)
        tmp_profs.set_prof('Sig', Sig, 0, gp)

        MtoL = cube[off]
        off += 1
    else:
        rhostarpar = np.zeros(gp.nrho)
        MtoL = 0.

    for pop in np.arange(1, gp.pops+1):  # [1, 2, ..., gp.pops]
        offstep = gp.nrho
        nupar = np.array(cube[off:off+offstep])
        tmp_nrnu = 1.*nupar[1+1+gp.nexp:-gp.nexp-1]
        tmp_profs.set_prof('nrnu', tmp_nrnu, pop, gp)
        tmp_nu = phys.rho(gp.xepol, nupar, pop, gp)
        tmp_profs.set_prof('nu', tmp_nu[gp.nexp:-gp.nexp], pop, gp)
        tmp_Signu = glp.rho_INTIPOL_Sig(gp.xepol, tmp_nu, gp) # [Munit/pc^2]
        tmp_Sig = gh.linipollog(gp.xepol, tmp_Signu, gp.xipol)
        tmp_profs.set_prof('Sig', tmp_Sig, pop, gp)
        off += offstep

        offstep = gp.nbeta
        if gp.chi2_Sig_converged:
            betapar = np.array(cube[off:off+offstep])
            tmp_beta, tmp_betastar = phys.beta(gp.xipol, gp.x0turn, betapar, gp)
            
            if check_beta(tmp_beta, gp):
                LOG(2, 'beta error')
                tmp_profs.chi2 = gh.err(1., gp)
                return tmp_profs
            tmp_profs.set_prof('beta', tmp_beta, pop, gp)
            tmp_profs.set_prof('betastar', tmp_betastar, pop, gp)

            #try:
            if gp.checksig:
                import gl_analytic as ga
                anrho = ga.rho_gaia(gp.xfine, gp)[0]
                rhopar_half = anrho[np.argmin(np.abs(gp.xfine-gp.Xscale[0]))]
                nr = -gh.derivipol(np.log(anrho), np.log(gp.xfine))
                dlr = np.hstack([nr[0], nr, nr[-1]])
                rhopar = np.hstack([rhopar_half, dlr])
                rhostarpar = 0.0*rhopar
                MtoL = 0.0
                betapar = np.array([ -4.29471374e-02,   1.92468972e+00,  -1.39137160e+00,\
                                     4.93307888e-01,  -9.10293049e-02,   8.37852613e-03,\
                                     -3.03557008e-04])
                # we get the right betastar profile again
                #ba = ga.beta_gaia(gp.xfine, gp)[pop]
                #ba = phys.beta2betastar(ba)
                #plot(gp.xfine, ba, 'b.-')
                #plot(gp.xfine, phys.betastar(gp.xfine, gp.x0turn, betapar, gp), 'r.-')
                #pdb.set_trace()

                # rho
                #loglog(gp.xfine, anrho, 'b.-')
                #loglog(gp.xfine, phys.rho(gp.xfine, rhopar, 0, gp), 'r.-')
                #pdb.set_trace()

                #Sig = glp.rho_param_INT_Sig(gp.xfine, nupar, pop, gp)
                #loglog(gp.xipol, gp.dat.Sig[1], 'b.-')
                #loglog(gp.xfine[gp.nexp:-gp.nexp], Sig[gp.nexp:], 'r.-')

                # nupar
                annu = ga.rho_gaia(gp.xfine, gp)[1]
                nupar_half = annu[np.argmin(np.abs(gp.xfine-gp.Xscale[1]))]
                nrnu = -gh.derivipol(np.log(annu), np.log(gp.xfine))
                dlrnu = np.hstack([nrnu[0], nrnu, nrnu[-1]])
                nupar = np.hstack([nupar_half, dlrnu])
                #loglog(gp.xfine, phys.rho(gp.xfine, nupar, 1, gp), 'r.-')
                #loglog(gp.xipol, gp.dat.nu[0], 'b.-')
                

            # TODO change gp.xfine back to gp.xepol

            sig,kap,zetaa,zetab=phys.sig_kap_zet(gp.xfine, rhopar, rhostarpar, MtoL, nupar, betapar, pop, gp)
            plot(gp.xipol, gp.dat.sig[pop], 'b.-')
            plot(gp.xfine[gp.nexp:-gp.nexp], sig, 'r.-')
            pdb.set_trace()
            
            #except Exception as detail:
            #    LOG('sigma error')
            #    tmp_profs.chi2 = gh.err(2., gp)
            #    return tmp_profs
            tmp_profs.set_prof('sig', sig, pop, gp)
            tmp_profs.set_prof('kap', kap, pop, gp)
            tmp_profs.set_zeta(zetaa, zetab, pop)
        off += offstep # still do this even if gp.chi2_Sig_converged is False

    if off != gp.ndim:
        LOG(1, 'wrong subscripts in gl_loglike')
        pdb.set_trace()

    # determine log likelihood
    chi2 = calc_chi2(tmp_profs, gp)
    LOG(1, 'found log likelihood = ', -chi2/2.)
    tmp_profs.chi2 = chi2

    return tmp_profs
## \fn geom_loglike(cube, ndim, nparams, gp)
# define log likelihood function to be called by pyMultinest and plot_profiles
# spherical version
# @param cube parameter cube as defined by gl_class_cube, in physical space already (not [0,1] cube anymore)
# @param ndim number of dimensions, needed as argument by multinest
# @param nparams number of parameters, needed as argument by multinest
# @param gp global parameters
