#!/usr/bin/env python3

##
# @file
# check all parameters for prior constraints

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
# @param nr -d ln(\rho)/d ln(r)


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


def check_nu(nu):
    return check_rho(nu)
## \fn check_nu(nu)
# check that tracer density is not jumping wildly


def check_beta(beta, gp):
    # now checking beta <= 1
    if max(beta)>1.:
        print('max beta!')              # TODO: remove, should be done by my_prior already
        return True
    # TODO: check smoothness of beta

    # now checking physical kappa: g(rvar, rfix, beta, dbetadr) >= 0
    if gp.usekappa == False:
        return False
    r0 = gp.xipol
    dR = r0[1:]-r0[:-1]
    r0extl = np.array([r0[0]/6., r0[0]/5., r0[0]/4., r0[0]/3., r0[0]/2., r0[0]/1.5])
    
    # extrapolation to the right (attention, could overshoot)
    dr0 = (r0[-1]-r0[-2])/8.
    r0extr = np.hstack([r0[-1]+dr0, r0[-1]+2*dr0, r0[-1]+3*dr0, r0[-1]+4*dr0])

    r0nu = np.hstack([r0extl, r0, dR/2.+r0[:-1], r0extr])
    r0nu.sort()
    tck0 = splrep(r0,beta*(r0**2+np.median(r0)**2), k=1, s=0.) # previous: k=2, s=0.1
    betanu = splev(r0nu,tck0)/(r0nu**2+np.median(r0)**2)

    drspl = splev(r0nu,tck0, der=1)
    dbetanudr = (drspl-betanu*2*r0nu)/(r0nu**2+np.median(r0)**2)
    for i in range(len(r0nu)-4):
        for j in range(i+1,len(r0nu)):
            if g(r0nu[j], r0nu[i], betanu[j], dbetanudr[j]) < 0:
                return True

    return False
## \fn check_beta(beta, gp)
# check that beta is bound and not jumping
# @param beta physical beta, to be in the range ]-infty,1]

def check_bprior(rhocheck, nucheck):
    for jj in range(len(rhocheck)):
        if rhocheck[jj] < nucheck[jj]:
            # gpl.clf(); gpl.yscale('log'); gpl.xscale('log')
            # gpl.plot(gp.xepol, rhocheck, 'k')
            # gpl.plot(gp.xepol, nucheck, 'b')
            # gpl.axvline(gp.rstarhalf)
            # pdb.set_trace()
            return True
    #print('NO bprior!')
    return False
## \fn check_bprior(rhocheck, nucheck)
# check that observed tracer mass is less than total mass
# @param rhocheck density profile in [Msun/pc^3]
# @param nucheck density profile in [Msun/pc^3]
