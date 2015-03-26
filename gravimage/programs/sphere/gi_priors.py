#!/usr/bin/env ipython3

##
# @file
# check all parameters for prior constraints

# (c) GPL v3 2015 Pascal Steger, pascal@steger.aero

import pdb
import numpy as np
from scipy.interpolate import splrep, splev
from gi_int import g
from gi_helper import LOG


def check_beta(beta, gp):
    # now checking beta <= 1
    if max(beta)>1.:
        LOG(2, 'max beta!')
        return True

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
    splpar_bint = splrep(r0, beta*(r0**2+np.median(r0)**2), k=1, s=0.) # previous: k=2, s=0.1
    betanu = splev(r0nu, splpar_bint)/(r0nu**2+np.median(r0)**2)

    drspl = splev(r0nu, splpar_bint, der=1)
    dbetanudr = (drspl-betanu*2*r0nu)/(r0nu**2+np.median(r0)**2)
    for i in range(len(r0nu)-4):
        for j in range(i+1,len(r0nu)):
            if g(r0nu[j], r0nu[i], betanu[j], dbetanudr[j]) < 0:
                return True

    return False
## \fn check_beta(beta, gp)
# check that beta is bound and not jumping
# @param beta physical beta, to be in the range ]-infty,1]
# @param gp global parameters
