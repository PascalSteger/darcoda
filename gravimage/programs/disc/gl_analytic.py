#!/usr/bin/env python3

##
# @file
# @ingroup gravimage
# all analytic profiles from gl_physics

# (c) GPL v3 2014 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import pdb
from scipy.integrate import simps
#import numpy as np
import numpy.random as npr
import gl_helper as gh

asech = lambda x: np.arccosh(1./x)
asec  = lambda x: np.arccos(1./x)

def get_prof(prof, pop, gp):
    zmin = 100.                               # [pc], first bin center
    zmax = 1300.                              # [pc], last bin center
    # get Stuetzpunkte for theoretical profiles (not yet stars, finer spacing in real space)
    nth = gp.nipol                            # [1] number of bins
    zth = 1.* np.arange(nth) * (zmax-zmin)/(nth-1.) + zmin # [pc] bin centers
    z0  = 240.                                # [pc], scaleheight of first population
    z02 = 200.                                # [pc], scaleheight of second population
    D   = 250.                                # [pc], scaleheight of all stellar tracers
    K   = 1.65                                # [TODO]
    F   = 1.65e-4                             # [TODO]
    C   = 17.**2.                             # [km/s] integration constant in sig

    # Draw mock data from exponential disk:
    nu_zth = np.exp(-zth/z0)                                 # [nu0] = [Msun/A/pc] 3D tracer density
    if prof == 'nu' and pop==1:
        return zth, nu_zth
    Kz_zth = -(K*zth/np.sqrt(zth**2.+D**2.) + 2.0 * F * zth) # [TODO]

    if gp.adddarkdisc:
        DD = 600                                         # [pc] scaleheight of dark disc
        KD = 0.15 * 1.650                                # [TODO]
        Kz_zth = Kz_zth - KD*zth/np.sqrt(zth**2. + DD**2.) # [TODO]

    # calculate sig_z^2
    inti = np.zeros(nth)
    for i in range(1, nth):
        inti[i] = simps(Kz_zth[:i]*nu_zth[:i], zth[:i])

    sigzth = np.sqrt((inti + C) / nu_zth)
    if prof == 'sig' and pop == 1:
        return zth, sigzth
    # project back to positions of stars
    ran = npr.uniform(size=int(gp.ntracer[1-1]))                 # [1]
    zstar = -z0 * np.log(1.0 - ran)           # [pc] stellar positions, exponential falloff

    sigzstar = gh.ipol(zth, sigzth, zstar)
    # > 0 ((IDL, Justin)) stellar velocity dispersion

    # assign [0,1] * maxsig
    ran2 = npr.normal(size=int(gp.ntracer[2-1]))  # [1]
    vzstar = ran2 * sigzstar                      # [km/s]

    # Add second population [thick-disc like]:
    if gp.pops == 2:
        nu_zth2 = gp.ntracer[2-1]/gp.ntracer[1-1]*np.exp(-zth/z02)
        if prof == 'nu' and pop == 2:
            return zth, nu_zth2
        # [nu0,2] = [Msun/A/pc], 3D tracer density, exponentially falling
        # no normalization to 1 done here
        inti    = np.zeros(nth)
        for i in range(1, nth):
            inti[i] = simps(Kz_zth[:i]*nu_zth2[:i], zth[:i])
        sigzth2 = np.sqrt((inti + C) / nu_zth2) # same integration constant
        if prof == 'sig' and pop == 2:
            return zth, sigzth2
        ran = npr.uniform(-1., 1., gp.ntracer[2-1])            # [1]
        zstar2 = -z02 * np.log(1.0 - ran)                      # [pc]
        #zstarobs = np.hstack([zstar, zstar2]) # concat pop1, pop2 for all stars
        sigzstar2 = gh.ipol(zth, sigzth2, zstar2)
        ran2 = npr.normal(-1., 1, gp.ntracer[2-1])        # [1]
        vzstar2 = ran2 * sigzstar2                        # [(km/2)^2]
## \fn get_prof(prof, pop, gp)
# get analytic profile prof for global parameters gp
# @param prof profile string, rho, nu, sig
# @param pop population int 1, 2
# @param gp global parameters
# @return z, profile
