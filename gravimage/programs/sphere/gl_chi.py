#!/usr/bin/env python3

##
# @file
# all functions called directly from gravimage

# (c) GPL v3 ETHZ 2014 Pascal Steger, psteger@phys.ethz.ch

import types
# TODO find missing modules

import ipdb
import numpy as np
import gl_helper as gh

def chi2red(model, data, sig, dof):
    # if Degrees Of Freedom = 1, return non-reduced chi2
    model = np.array(model)
    data  = np.array(data)
    sig   = np.array(sig)
    return np.sum(((model-data)**2./sig**2.)/dof)
## \fn chi2red(model, data, sig, dof)
# determine 'reduced chi2'
# @param model profile
# @param data profile
# @param sig spread
# @param dof Degrees Of Freedom

def calc_chi2(profs, gp):
    chi2 = 0.
    # calc chi^2 from overall baryonic tracers
    # is not needed, as photometric data has very small errors,
    # and the corresponding M/L uncertainty is automatically
    # penalized from wrong sigma_LOS

    # now run through the stellar tracers
    for pop in np.arange(1, gp.pops+1): # [1, 2, ... , pops]
        Sigdat   = gp.dat.Sig[pop]      # [Munit/pc^2]
        Sigerr   = gp.dat.Sigerr[pop]   # [Munit/pc^2]
        Sigmodel = profs.get_prof('Sig', pop)
        chi2_Sig  = chi2red(Sigmodel, Sigdat, Sigerr, gp.nipol) # [1]
        chi2 += chi2_Sig                 # [1]
        gh.LOG(1, ' chi2_Sig   = ', chi2_Sig)

        # use the following only if chi2_nu_converged used rather than Sig_converged
        #nudat   = gp.dat.nu[pop]      # [Munit/pc^2]
        #nuerr   = gp.dat.nuerr[pop]   # [Munit/pc^2]
        #numodel = profs.get_prof('nu', pop)
        #chi2_nu  = chi2red(numodel, nudat, nuerr, gp.nipol) # [1]
        #chi2 += chi2_nu                 # [1]
        #gh.LOG(1, ' chi2_nu   = ', chi2_nu)

        if gp.chi2_Sig_converged > 0:
            continue

        sigdat  = gp.dat.sig[pop]    # [km/s]
        sigerr  = gp.dat.sigerr[pop]    # [km/s]
        chi2_sig = chi2red(profs.get_prof('sig', pop), sigdat, sigerr, gp.nipol) # [1]
        chi2 += chi2_sig                # [1]
        gh.LOG(1, '  chi2_sig  = ', chi2_sig)
        if gp.usekappa:
            kapdat  = gp.dat.kap[pop] # [1]
            kaperr  = gp.dat.kaperr[pop] # [1]
            chi2_kap = chi2red(profs.get_kap(pop), kapdat, kaperr, gp.nipol) # [1]
            chi2 += chi2_kap            # [1]

        if gp.usezeta:
            zetaadat = gp.dat.zetaadat[pop]
            zetabdat = gp.dat.zetabdat[pop]
            zetaaerr = gp.dat.zetaaerr[pop]
            zetaberr = gp.dat.zetaberr[pop]
            zetaa_model, zetab_model = profs.get_zeta(pop)
            chi2_zetaa = chi2red(zetaa_model, zetaadat, zetaaerr, 1)
            chi2_zetab = chi2red(zetab_model, zetabdat, zetaberr, 1)
            chi2 += (chi2_zetaa + chi2_zetab)

    if gp.chi2_Sig_converged > 0:
        chi2 *= 10 # overamplify chi2 to get better models after switch
        if chi2 < gp.chi2_switch:
            gp.chi2_Sig_converged -= 1
            gh.LOG(1, 'Sig finished burn-in, waiting to get stable, ', gp.chi2_Sig_converged)
    return chi2
## \fn calc_chi2(profs, gp)
# Calculate chi^2
# @param profs profiles for rho, M, rho*, nu_i, beta_i, sig_i, kap_i
# @param gp global parameters
