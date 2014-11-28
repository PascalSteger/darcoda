#!/usr/bin/env ipython3

##
# @file
# all functions called directly from gravlite
# disc version, done

# (c) ETHZ 2013 Pascal Steger, psteger@phys.ethz.ch

#from types import *
import pdb
import numpy as np
import gl_helper as gh

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

    #Tracer population comparison

    for pop in range(0, gp.ntracer_pops):
        nudat    = gp.dat.nu[pop]
        nuerr    = gp.dat.nuerr[pop]
        numodel  = profs.get_prof('nu_vec', pop)
        chi2_nu = chi2red(numodel, nudat, nuerr, gp.nbins)
        gh.LOG(1, ' chi2_nu0 = ', chi2_nu)
        chi2 +=chi2_nu

        if not gp.chi2_nu_converged:
            continue # with pop loop

        sigdat  = gp.dat.sig[pop]    # [km/s]
        sigerr  = gp.dat.sigerr[pop] # [km/s]
        sigmodel= profs.get_prof('sig_vec', pop)
        chi2_sig = chi2red(sigmodel, sigdat, sigerr, gp.nbins) # [1]
        if chi2_sig == np.inf:
            print('chi2_sig has become infinite')
            pdb.set_trace()
        chi2 += chi2_sig             # [1]
        gh.LOG(1, '  chi2_sig  = ', chi2_sig)

    # switch to chi2_sig calculation too, if converged on Sig
    if not gp.chi2_nu_converged:
        chi2 *= 10
        if chi2 < gp.chi2_switch:
            gp.chi2_switch_counter +=1
            print('chi2 less than switch found #######################################')
        if gp.chi2_switch_counter>= gp.chi2_switch_mincount:
            gh.LOG(1, 'nu burn-in finished, switching on sigma')
            gp.chi2_nu_converged = True

        #if chi2 < gp.chi2_switch:
        #    pdb.set_trace()
        #    gh.LOG(1, 'nu burn-in finished, switching on sigma')
        #    gp.chi2_nu_converged = True

    return chi2
## \fn calc_chi2(profs)
# Calculate chi^2
# @param profs profiles for rho, M, nu_i, beta_i, sig_i, kap_i
