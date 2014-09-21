#!/usr/bin/env ipython3

##
# @file
# all functions called directly from gravlite
# disc version, done

# (c) ETHZ 2013 Pascal Steger, psteger@phys.ethz.ch

from types import *
import pdb
import numpy as np

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


def calc_chi2(profs, gp):
    chi2 = 0.
    off = 0
    # TODO: include rho*?
    for pop in np.arange(gp.pops)+1: # look at 0 (==1) for pop=1, and (1,2) for pop==2
        nudat   = gp.dat.nu[pop]     # [Munit/pc^3]
        nuerr   = gp.dat.nuerr[pop]  # [Munit/pc^3]
        numodel = profs.get_prof('nu', pop)
        chi2_nu  = chi2red(numodel, nudat, nuerr, gp.nipol)
        chi2 += chi2_nu              # [1]
        
        sigdat  = gp.dat.sig[pop]    # [km/s]
        sigerr  = gp.dat.sigerr[pop] # [km/s]
        sigmodel= profs.get_prof('sig', pop)
        chi2_sig = chi2red(sigmodel, sigdat, sigerr, gp.nipol) # [1]
        chi2 += chi2_sig             # [1]
        print('chi2_nu, chi2_sig = ', chi2_nu, ' ', chi2_sig)
        if gp.usekappa:
            kapdat  = gp.dat.kap[pop]       # [1]
            kaperr  = gp.dat.kaperr[pop]    # [1]
            chi2_kap = chi2red(profs.get_kap(pop), kap, kaperr, gp.nipol) # [1]
            chi2 += chi2_kap                                     # [1]
    return chi2
## \fn calc_chi2(profs)
# Calculate chi^2
# @param profs profiles for rho, M, nu_i, beta_i, sig_i, kap_i

