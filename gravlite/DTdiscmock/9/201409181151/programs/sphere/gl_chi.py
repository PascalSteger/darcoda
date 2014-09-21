#!/usr/bin/env ipython3

##
# @file
# all functions called directly from gravlite

# (c) ETHZ 2013 Pascal Steger, psteger@phys.ethz.ch

from types import *
import pdb
import numpy as np
import numpy.random as npr

import gl_analytic as ga
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

    if gp.investigate =='obs':
        Sigdat = gp.dat.Sig[0]              # [Munit/pc^2] from rho*
        Sigerr = gp.dat.Sigerr[0]           # [Munit/pc^2]
        Sigmodel = profs.get_prof('Sig', 0)
        chi2_Sigstar = chi2red(Sigmodel, Sigdat, Sigerr, gp.nipol) # [1]
        chi2 += chi2_Sigstar
        gh.LOG(1, 'chi2_Sigstar  = ', chi2_Sigstar)

    # now run through the stellar tracers
    for pop in np.arange(1, gp.pops+1): # [1, 2, ... , pops]
        Sigdat   = gp.dat.Sig[pop]      # [Munit/pc^2]
        Sigerr   = gp.dat.Sigerr[pop]   # [Munit/pc^2]
        Sigmodel = profs.get_prof('Sig', pop)
        chi2_Sig  = chi2red(Sigmodel, Sigdat, Sigerr, gp.nipol) # [1]
        chi2 += chi2_Sig                 # [1]
        gh.LOG(1, ' chi2_Sig   = ', chi2_Sig)

        if not gp.chi2_Sig_converged:
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

    if not gp.chi2_Sig_converged:
        chi2 *= 10
        if chi2 < gp.chi2_switch:
            gh.LOG(1, 'Sig finished burn-in, switching on sigma!')
            gp.chi2_Sig_converged = True
    return chi2
## \fn calc_chi2(profs, gp)
# Calculate chi^2
# @param profs profiles for rho, M, rho*, nu_i, beta_i, sig_i, kap_i
# @param gp global parameters

