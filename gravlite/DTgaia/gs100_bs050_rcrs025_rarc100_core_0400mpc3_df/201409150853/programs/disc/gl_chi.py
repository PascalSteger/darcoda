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

    # include rho*
    Sigdat   = gp.dat.Sig[0]     # [Munit/pc^3]
    Sigerr   = gp.dat.Sigerr[0]  # [Munit/pc^3]
    Sigmodel = profs.get_prof('Sig', 0)
    chi2_Sig  = chi2red(Sigmodel, Sigdat, Sigerr, gp.nipol)
    print('chi2_Sig = ',chi2_Sig)
    chi2 += chi2_Sig              # [1]


    for pop in np.arange(gp.pops)+1: # look at 0 (==1) for pop=1, and (1,2) for pop==2
        Sigdat   = gp.dat.Sig[pop]     # [Munit/pc^3]
        Sigerr   = gp.dat.Sigerr[pop]  # [Munit/pc^3]
        Sigmodel = profs.get_prof('Sig', pop)
        chi2_Sig  = chi2red(Sigmodel, Sigdat, Sigerr, gp.nipol)
        chi2 += chi2_Sig              # [1]

        if not gp.chi2_Sig_converged:
            continue # with pop loop

        sigdat  = gp.dat.sig[pop]    # [km/s]
        sigerr  = gp.dat.sigerr[pop] # [km/s]
        sigmodel= profs.get_prof('sig', pop)
        chi2_sig = chi2red(sigmodel, sigdat, sigerr, gp.nipol) # [1]
        chi2 += chi2_sig             # [1]
        print('chi2_Sig, chi2_sig = ', chi2_Sig, ' ', chi2_sig)
        if gp.usekappa:
            kapdat  = gp.dat.kap[pop]       # [1]
            kaperr  = gp.dat.kaperr[pop]    # [1]
            chi2_kap = chi2red(profs.get_kap(pop), kap, kaperr, gp.nipol) # [1]
            chi2 += chi2_kap                                     # [1]

    # switch to chi2_sig calculation too, if converged on Sig
    if not gp.chi2_Sig_converged:
        chi2 *= 10
        if chi2 < gp.chi2_switch:
            print('Sig burn-in finished, switching on sigma')
            gp.chi2_Sig_converged = True

    return chi2
## \fn calc_chi2(profs)
# Calculate chi^2
# @param profs profiles for rho, M, nu_i, beta_i, sig_i, kap_i

