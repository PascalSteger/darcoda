#!/usr/bin/env python3

##
# @ file
# generate simple disc data ourselves

# (c) GPL v3 2014 Pascal Steger, psteger@phys.ethz.ch


import numpy as np
import numpy.random as npr
import pdb
from scipy.integrate import simps

import gl_units as gu
import gl_helper as gh


def disc_mock(gp):
    global K,C,D,F, zth, zp_kz, zmin, zmax, z0, z02
    # Set up simple population here using analytic formulae:
    zmin = 100.                               # [pc], first bin center
    zmax = 1300.                              # [pc], last bin center
    # get Stuetzpunkte for theoretical profiles (not yet stars, finer spacing in real space)
    nth = 3*gp.nipol
    zth = 1.* np.arange(nth) * (zmax-zmin)/(nth-1.) + zmin
    z0  = 240.                                # [pc], scaleheight of first population
    z02 = 200.                                # [pc], scaleheight of second population
    D   = 250.                                # [pc], scaleheight of all stellar tracers
    K   = 1.65                                # [TODO]
    F   = 1.65e-4                             # [TODO]
    C   = 17.**2.                             # [km/s] integration constant in sig

    # Draw mock data:
    nu_zth = np.exp(-zth/z0)                                 # [1]
    Kz_zth = -(K*zth/np.sqrt(zth**2.+D**2.) + 2.0 * F * zth) # [TODO]

    if gp.adddarkdisc:
        DD = 600                                         # [pc]
        KD = 0.15 * 1.650                                # [TODO]
        Kz_zth = Kz_zth - KD*zth/np.sqrt(zth**2. + DD**2.) # [TODO]

    # calculate sig_z^2
    inti = np.zeros(nth)
    for i in range(1, nth):
        inti[i] = simps(Kz_zth[:i]*nu_zth[:i], zth[:i])

    sigzth = np.sqrt((inti + C) / nu_zth)

    # project back to positions of stars
    ran = npr.uniform(size=int(gp.ntracer[1-1]))                 # [1]
    zstar = -z0 * np.log(1.0 - ran)           # [pc] stellar positions
    sigzstar = gh.ipol(zth, sigzth, zstar)
    # > 0 ((IDL, Justin)) stellar velocity dispersion

    # assign [0,1] * maxsig
    ran2 = npr.normal(size=int(gp.ntracer[1-1]))  # [1]
    vzstar = ran2 * sigzstar                      # [km/s]

    # Add second population [thick-disc like]:
    if gp.pops == 2:
        nu_zth2 = gp.ntracer[2-1]/gp.ntracer[1-1]*np.exp(-zth/z02)
        # no normalization to 1
        inti    = np.zeros(nth)
        for i in range(1, nth):
            inti[i] = simps(Kz_zth[:i]*nu_zth2[:i], zth[:i])
        sigzth2 = np.sqrt((inti + C) / nu_zth2) # same integration constant
        ran = npr.uniform(-1., 1., gp.ntracer[2-1])            # [1]
        zstar2 = -z02 * np.log(1.0 - ran)                      # [pc]
        zstarobs = np.hstack([zstar, zstar2]) # concat pop1, pop2 for all stars
        sigzstar2 = gh.ipol(zth, sigzth2, zstar2)
        ran2 = npr.normal(-1., 1, gp.ntracer[2-1])        # [1]
        vzstar2 = ran2 * sigzstar2                        # [(km/2)^2]

    # enforce observational cut on zmax:
    sel = (zstar < zmax)
    print('fraction of z<zmax selected elements: ', 1.*sum(sel)/(1.*len(sel)))
    z_dat  = zstar[sel]
    vz_dat = vzstar[sel]

    # throw away velocities of value zero (unstable?):
    sel = (abs(vz_dat) > 0)
    print('fraction of vz_dat>0 selected elements: ', 1.*sum(sel)/(1.*len(sel)))
    z_dat  = z_dat[sel]
    vz_dat = vz_dat[sel]

    # Calulate binned data (for plots/binned anal.). old way, linear spacings, no const #particles/bin
    binmin, binmax, z_dat_bin, sig_dat_bin, count_bin = gh.binsmooth(z_dat, vz_dat, \
                                                                     zmin, zmax, gp.nipol, 0.)
    sig_dat_err_bin = sig_dat_bin / np.sqrt(count_bin)

    nu_dat_bin, count_bin = gh.bincount(z_dat, binmax)
    nu_dat_err_bin = nu_dat_bin / np.sqrt(count_bin)
    renorm = max(nu_dat_bin)
    nu_dat_bin = nu_dat_bin / renorm
    nu_dat_err_bin = nu_dat_err_bin / renorm

    # if only 1 pop, use 0 for all components
    binmin0 = binmin; binmax0 = binmax; z_dat_bin0 = z_dat_bin
    sig_dat_bin0 = sig_dat_bin; sig_dat_err_bin0 = sig_dat_err_bin
    nu_dat_bin0  = nu_dat_bin;  nu_dat_err_bin0  = nu_dat_err_bin

    if gp.pops == 2:
        # enforce observational constraint on z<z_max
        sel = (zstar2 < zmax)
        z_dat2  = zstar2[sel]
        vz_dat2 = vzstar2[sel]

        # cut zero velocities:
        sel = (abs(vz_dat2) > 0)
        z_dat2  = z_dat2[sel]
        vz_dat2 = vz_dat2[sel]

        # Calulate binned data (for plots/binned analysis):
        binmin2, binmax2, z_dat_bin2, sig_dat_bin2, count_bin2 = gh.binsmooth(z_dat2, vz_dat2, \
                                                                              zmin, zmax, gp.nipol, 0.)
        sig_dat_err_bin2 = sig_dat_bin2 / np.sqrt(count_bin2)

        nu_dat_bin2, count_bin2 = gh.bincount(z_dat2, binmax2)
        nu_dat_err_bin2 = nu_dat_bin2 / np.sqrt(count_bin2)
        renorm2 = max(nu_dat_bin2) # normalize by max density of first bin, rather
        nu_dat_bin2 = nu_dat_bin2 / renorm2
        nu_dat_err_bin2 = nu_dat_err_bin2 / renorm2

        # CALCULATE PROPERTIES FOR ALL POP TOGETHER
        z_dat0 = np.hstack([z_dat, z_dat2])
        vz_dat0 = np.hstack([vz_dat, vz_dat2])

        # Calulate binned data (for plots/binned anal.). old way, linear spacings, no const #particles/bin
        binmin0, binmax0, z_dat_bin0, sig_dat_bin0, count_bin0 = gh.binsmooth(z_dat0, vz_dat0, \
                                                                              zmin, zmax, gp.nipol, 0.)
        sig_dat_err_bin0 = sig_dat_bin0 / np.sqrt(count_bin0)
        # binmin, binmax, z_dat_bin = gh.bin_r_const_tracers(z_dat, gp.nipol) # TODO: enable, get sig2

        nu_dat_bin0, count_bin0 = gh.bincount(z_dat0, binmax0)
        nu_dat_err_bin0 = nu_dat_bin0 / np.sqrt(count_bin0)
        renorm0 = max(nu_dat_bin0)
        nu_dat_bin0 = nu_dat_bin0 / renorm0
        nu_dat_err_bin0 = nu_dat_err_bin0 / renorm0


    xip = np.copy(z_dat_bin0)                        # [pc]
    gp.dat.binmin = binmin0
    gp.dat.rbin   = xip
    gp.xipol      = gp.dat.rbin
    gp.Rscale.append(D)  # [pc]
    gp.Rscale.append(z0) # [pc]
    maxr          = max(gp.dat.rbin)
    gp.xepol      = np.hstack([gp.dat.rbin, 2*maxr, 4*maxr, 8*maxr])
    gp.dat.binmax = binmax0
    gp.dat.Mrdat   = K*xip/np.sqrt(xip**2.+D**2.) / (2.0*np.pi*gu.G1__pcMsun_1km2s_2)
    gp.dat.Mrerr   = gp.dat.Mrdat*nu_dat_err_bin/nu_dat_bin

    gp.dat.nu.append(nu_dat_bin0)        # [Msun/pc^3], normalized to 1 at center
    gp.dat.nuerr.append(nu_dat_err_bin0) # [Msun/pc^3], normalized

    gp.dat.sig.append(sig_dat_bin0)       # [km/s]
    gp.dat.sigerr.append(sig_dat_err_bin0)# [km/s]

    # TODO: define as rho* = (M/L) * L* as in spherical case?
    gp.dat.nu.append(nu_dat_bin)         # [Msun/pc^3], normalized to 1
    gp.dat.nuerr.append(nu_dat_err_bin)  # [Msun/pc^3], normalized

    gp.dat.sig.append(sig_dat_bin)        # [km/s]
    gp.dat.sigerr.append(sig_dat_err_bin) # [km/s]

    if gp.pops == 2:
        gp.Rscale.append(z02)                 # [pc]
        gp.dat.nu.append(nu_dat_bin2)        # [Msun/pc^3], normalized to 1
        gp.dat.nuerr.append(nu_dat_err_bin2) # [Msun/pc^3], normalized

        gp.dat.sig.append(sig_dat_bin2)       # [km/s]
        gp.dat.sigerr.append(sig_dat_err_bin2)# [km/s]
    return gp.dat

## \fn disc_mock(gp)
# generate disc data from analytic form, return 3D densities, delta (=tilt)
# @param gp global parameters
