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
from pylab import *
ion()

def geom_loglike(cube, ndim, nparams, gp):
    tmp_profs = Profiles(gp.ntracer_pops, gp.nbins, gp.nrhonu, gp.nbaryon_pops, gp.nbaryon_params)

    #Normalisation constant C for sigz calculation
    off = 0
    offstep = 1
    norm = cube[off]
    off += offstep

    #Dark Matter rho parameters (rho_C, kz_C, kz_vector, kz_LS)
    offstep = gp.nrhonu + 1
    rho_DM_params = np.array(cube[off:off+offstep])
    rho_DM_C = rho_DM_params[0] #rho_C
    kz_rho_DM_allz = rho_DM_params[1:] #kz for rho across all z points [0, bin_centres, LS]
    tmp_rho_DM_allz = phys.rho(gp.z_all_pts, kz_rho_DM_allz, rho_DM_C) #outputs rho across all points

    tmp_profs.rho_DM_C = tmp_rho_DM_allz[0]
    tmp_profs.set_prof('rho_DM_vec', tmp_rho_DM_allz[1:-1], 0, gp)
    tmp_profs.rho_DM_LS = tmp_rho_DM_allz[-1]
    off += offstep

    #Baryons
    for bary_pop in range(0, gp.nbaryon_pops):
        offstep = gp.nbaryon_params
        bary_params = np.array(cube[off:off+offstep])
        off += offstep

    #Tracer params, nu_C, kz_nu_C, kz_nu_vector, kz_nu_LS
    for tracer_pop in range(0, gp.ntracer_pops):
        offstep = gp.nrhonu + 1
        tracer_params = np.array(cube[off:off+offstep])
        nu_C = tracer_params[0]
        kz_nu_allz = tracer_params[1:] #kz for rho across all z points [0, bin_centres, LS]
        tmp_nu_allz = phys.rho(gp.z_all_pts, kz_nu_allz, nu_C) #outputs nu across all z points
        tmp_profs.nu_C = tmp_nu_allz[0]
        tmp_profs.set_prof('nu_vec', tmp_nu_allz[1:-1], tracer_pop, gp)
        tmp_profs.nu_LS = tmp_nu_allz[-1]
        off += offstep

    if off != gp.ndim:
        gh.LOG(1,'wrong subscripts in gl_class_cube')
        raise Exception('wrong subscripts in gl_class_cube')

    #Calculate Sigma (surface density)
    Sig_DM_allz = phys.Sig(gp.z_all_pts, tmp_rho_DM_allz)
    tmp_profs.Sig_DM_C = Sig_DM_allz[0]
    tmp_profs.set_prof('Sig_DM_vec', Sig_DM_allz[1:-1], 0, gp)
    tmp_profs.Sig_DM_LS = Sig_DM_allz[-1]

    #Calculate sigma (velocity dispersion)
    sigz_vecLS = phys.sigz(gp.z_all_pts, Sig_DM_allz, tmp_nu_allz, norm)
    tmp_profs.set_prof('sig_vec', sigz_vecLS[0:-1], 0, gp)
    tmp_profs.sig_LS = sigz_vecLS[-1]

    # determine log likelihood
    chi2 = calc_chi2(tmp_profs, gp) #HS currently rewriting calc_chi2
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
