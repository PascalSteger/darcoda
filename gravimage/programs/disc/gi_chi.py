#!/usr/bin/env ipython3

##
# @file
# all functions called directly from gravimage
# disc version, done

# (c) GPL v3 ETHZ 2014 Pascal Steger, pascal@steger.aero

#from types import *
import pdb
import numpy as np
import gi_helper as gh

def chi2red(model, data, sig, dof):
    # if Degrees Of Freedom = 1, return non-reduced chi2
    model = np.array(model)
    data  = np.array(data)
    sig   = np.array(sig)
    chired = np.sum(((model-data)**2./sig**2.)/dof)
    return chired
## \fn chi2red(model, data, sig, dof)
# determine 'reduced chi2'
# @param model
# @param data
# @param sig spread
# @param dof Degrees Of Freedom


def calc_chi2(profs, gp):
    chi2 = 0.

    # include rho* in chi^2 calculation
    nudat    = gp.dat.nu[0]
    nuerr    = gp.dat.nuerr[0]
    numodel  = profs.get_prof('nu', 0)
    chi2_nu = chi2red(numodel, nudat, nuerr, gp.nipol)
    gh.LOG(2, ' chi2_nu0 = ', chi2_nu)
    chi2 +=chi2_nu

    for pop in np.arange(1,gp.pops+1): # look at pops 1, 2, ...
        nudat    = gp.dat.nu[pop]
        nuerr    = gp.dat.nuerr[pop]
        numodel  = profs.get_prof('nu', pop)
        chi2_nu = chi2red(numodel, nudat, nuerr, gp.nipol)
        gh.LOG(2, ' chi2_nu['+str(pop)+'] = ', chi2_nu)
        chi2 +=chi2_nu
        pdb.set_trace()
        if not gp.chi2_nu_converged:
            continue # with pop loop

        sigdat  = gp.dat.sig[pop]    # [km/s]
        sigerr  = gp.dat.sigerr[pop] # [km/s]
        sigmodel= profs.get_prof('sig', pop)
        chi2_sig = chi2red(sigmodel, sigdat, sigerr, gp.nipol) # [1]
        if chi2_sig == np.inf:
            print('chi2_sig has become infinite')
            pdb.set_trace()
        chi2 += chi2_sig             # [1]
        gh.LOG(1, '  chi2_sig  = ', chi2_sig)

        if gp.usekappa:
            kapdat  = gp.dat.kap[pop]       # [1]
            kaperr  = gp.dat.kaperr[pop]    # [1]
            chi2_kap = chi2red(profs.get_kap(pop), kapdat, kaperr, gp.nipol) # [1]
            chi2 += chi2_kap                                     # [1]

    # switch to chi2_sig calculation too, if converged on Sig
    if not gp.chi2_nu_converged:
        chi2 *= 10
        if chi2 < gp.chi2_switch:
            gh.LOG(1, 'nu burn-in finished, switching on sigma')
            gp.chi2_nu_converged = True

    return chi2
## \fn calc_chi2(profs)
# Calculate chi^2
# @param profs profiles for rho, M, nu_i, beta_i, sig_i, kap_i
