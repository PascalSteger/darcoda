#!/usr/bin/env python3

##
# @file
# check all parameters for prior constraints
# disc version

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import pdb
import numpy as np
from scipy.interpolate import splrep, splev
from gl_int import g
import gl_plot as gpl

def check_nr(nr, gp):
    r0 = gp.xepol/gp.rstarhalf
    if max(np.abs((nr[:-1]-nr[1:])/(r0[:-1]-r0[1:]))) > 20:
        return True
    return False
## \fn check_nr(nr, gp)
# check that n(r) is not jumping wildly
# @param nr -d ln(rho)/d ln(r)
# @param gp


def check_rho(rho, rprior, rhotol):
    if rprior:
        rightrho = rho[1:]
        leftrho  = rho[:-1]
        if sum(rightrho/leftrho > rhotol) > 0:
            return True
    return False
## \fn check_rho(rho)
# check that density is not jumping wildly
# @param rho density profile in [Msun/pc^3]


def check_tilt(tmp_tilt, gp):
    return False
## \fn check_tilt(tmp_tilt, gp)
# check that tilt is physically possible
# TODO: sensible values, e.g. percentage of K_z
# @param tmp_tilt tilt profile
# @param gp


def check_bprior(rhocheck, nucheck):
    for jj in range(len(rhocheck)):
        if rhocheck[jj] < nucheck[jj]:
            # gpl.clf(); gpl.yscale('log'); gpl.xscale('log')
            # gpl.plot(gp.xepol, rhocheck, 'k')
            # gpl.plot(gp.xepol, nucheck, 'b')
            # gpl.axvline(gp.rstarhalf)
            return True
    return False
## \fn check_bprior(rhocheck, nucheck)
# check that observed tracer mass is less than total mass
# @param rhocheck density profile in [Msun/pc^3]
# @param nucheck density profile in [Msun/pc^3]
