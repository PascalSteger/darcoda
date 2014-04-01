#!/usr/bin/env python3

##
# @file
# all functions called directly from gravlite
# disc version, done

# (c) ETHZ 2013 Pascal Steger, psteger@phys.ethz.ch

from types import *
import pdb
import numpy.random as npr

# from gl_analytic import *
# from gl_project import rho_INT_Rho, rho_param_INT_Rho
from gl_class_profiles import Profiles

def chi2red(model, data, sig, dof):
    # if Degrees Of Freedom = 1, return non-reduced chi2
    model = np.array(model)
    data  = np.array(data)
    sig   = np.array(sig)
    return np.sum(((model-data)**2./sig**2.)/dof)
## \fn chi2red(model, data, sig, dof)
# determine 'reduced chi2'
# @param model
# @param data
# @param sig spread
# @param dof Degrees Of Freedom


def calc_chi2(profs, nuparstore, gp):
    chi2 = 0.
    off = 0
    for pop in np.arange(gp.pops)+1: # look at 0 (==1) for pop=1, and (1,2) for pop==2
        nuparams = nuparstore[pop-1]
        Numodel  = rho_param_INT_Rho(gp.xepol, nuparams, gp) # [Msun/pc^2], on nipol bins
        # Numodel = rho_INT_Rho(gp.xipol, profs.get_nu(pop)) # [Msun/pc^2]
        Nudata  = gp.dat.Nudat[pop]                      # [Msun/pc^2]
        Nuerr   = gp.dat.Nuerr[pop]                      # [Msun/pc^2]
        
        chi2_nu  = chi2red(Numodel, Nudata, Nuerr, gp.dof) # [1]
        chi2 += chi2_nu                 # [1]
        
        sigdat  = gp.dat.sigdat[pop]    # [km/s]
        sigerr  = gp.dat.sigerr[pop]    # [km/s]
        chi2_sig = chi2red(profs.get_sig(pop), sigdat, sigerr, gp.dof) # [1]
        chi2 += chi2_sig                # [1]
        # print('profs: ', profs.get_nu(pop)[0],', data: ', Nudata[0])
        print('chi2_nu, chi2_sig = ', chi2_nu, ' ', chi2_sig)
        if gp.usekappa:
            kapdat  = gp.dat.kapdat[pop]       # [1]
            kaperr  = gp.dat.kaperr[pop]       # [1]
            chi2_kap = chi2red(profs.get_kap(pop), kapdat, kaperr, gp.dof) # [1]
            chi2 += chi2_kap                                     # [1]
    return chi2
## \fn calc_chi2(profs)
# Calculate chi^2
# @param profs profiles for rho, M, nu_i, beta_i, sig_i, kap_i

