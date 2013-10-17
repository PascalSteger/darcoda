#!/usr/bin/python3.2

##
# @ file
# generate simple disc data ourselves
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch


import numpy as np
import numpy.random as npr
import random
import scipy
import scipy.integrate
from scipy.integrate import simps,trapz
import pdb

import gl_params as gp
if gp.geom == 'disc':
    import physics_disc as phys
else:
    import physics_sphere as phys
import gl_helper as gh
from binsmooth import binsmoo
from bincount import bincou

global K,C,D,F, zth, zp_kz, zmin, zmax, z0, z02

# Set up simple population here using analytic formulae: 
zmin = 100. ; zmax = 1300.; zpnts = gp.nipol
zth = 1.* np.arange(zpnts) * (zmax-zmin)/(zpnts-1.) + zmin 
z0  = 240.
z02 = 200.
D   = 250.

K   = 1.65
F   = 1.65e-4
C   = 17.**2.



## generate disc data from analytic form, return 3D densities, delta (=tilt), ...
def disc_simple():
    # trick to speed up things considerably, after first data generated, if zpnts > 50 (or other low number)
    # gp.dat.load(gp.files.dir+'pp')
    # return gp.dat

    # Draw mock data: 
    nu_zth = np.exp(-zth/z0)              # [1]
    Kz_zth = -(K*zth/np.sqrt(zth**2.+D**2.) + 2.0 * F * zth)   # [TODO]

    if gp.adddarkdisc: # False
        DD = 600       # [pc]
        KD = 0.15 * 1.650                                # [TODO]
        Kz_zth = Kz_zth - KD*zth/np.sqrt(zth**2.+DD**2.) # [TODO]

    inti = np.zeros(zpnts) 
    for i in range(1,zpnts):
        inti[i] = simps(Kz_zth[:i]*nu_zth[:i],zth[:i])
    sigzth = np.sqrt((inti + C) / nu_zth)
    nstars = 10000.                                # [1]
    ran = npr.uniform(size=int(nstars))           # [1]
    zstar = -z0 * np.log(1.0 - ran)               # [pc]
    sigzstar = gh.ipol(zth,sigzth,zstar) # > 0 #selection? no, IDL GT means ">", so perhaps '>' is a shift
    ran2 = npr.normal(size=int(nstars))  # [1]
    vzstar = ran2 * sigzstar             # [TODO]

    # Add second population [thick-disc like]: 
    fac2 = 1.
    if gp.pops == 2:
        nu_zth2 = fac2*np.exp(-zth/z02)
        inti    = np.zeros(zpnts)
        for i in range(1,zpnts):
            inti[i] = simps(Kz_zth[:i]*nu_zth2[:i],zth[:i])
        sigzth2 = np.sqrt((inti + C) / nu_zth2)
        nstars2 = nstars*fac2             # [1]
        ran = npr.uniform(-1.,1.,nstars2) # [1]
        zstar2 = -z02 * np.log(1.0 - ran) # [pc]
        # zstar = [zstar,zstar2]
        sigzstar2 = gh.ipol(zth,sigzth2,zstar2)   # TODO: > 0?
        ran2 = npr.normal(-1.,1,nstars2)          # [1]
        vzstar2 = ran2 * sigzstar2                # [(km/2)^2]

    # Cut on zmax:
    # z_dat = [zstar((zstar < zmax)),zstar2((zstar2 < zmax))] 
    # vz_dat = [vzstar((zstar < zmax)),vzstar2((zstar2 < zmax))]
    sel = (zstar < zmax)
    z_dat = zstar[sel];   vz_dat = vzstar[sel]
    
    # Cut zero velocities:
    sel = (abs(vz_dat) > 0)
    z_dat = z_dat[sel];   vz_dat = vz_dat[sel]
    
    # Calulate binned data (for plots/binned anal.): 
    index = np.argsort(z_dat)

    z_dat_bin, sig_dat_bin, count_bin = binsmoo(z_dat[index], vz_dat[index], zmin, zmax, gp.nipol, 0.)
    sig_dat_bin = np.sqrt(sig_dat_bin)
    sig_dat_err_bin = sig_dat_bin / np.sqrt(count_bin)

    z_dat_bin, nu_dat_bin, count_bin = bincou(z_dat[index], zmin, zmax, gp.nipol)
    print(nu_dat_bin)
    nu_dat_err_bin = nu_dat_bin / np.sqrt(count_bin)
    renorm = max(nu_dat_bin)
    nu_dat_bin = nu_dat_bin / renorm
    nu_dat_err_bin = nu_dat_err_bin / renorm

    if gp.pops == 2:
        sel = (zstar2 < zmax)
        z_dat2  = zstar2[sel];       vz_dat2 = vzstar2[sel]
        
        # Cut zero velocities:
        sel = (abs(vz_dat2) > 0)
        z_dat2 = z_dat2[sel];        vz_dat2 = vz_dat2[sel]
        
        # Calulate binned data (for plots/binned anal.): 
        index2 = np.argsort(z_dat2)
        z_dat_bin2, sig_dat_bin2, count_bin2 = binsmoo(z_dat2[index2],vz_dat2[index2],zmin,zmax,gp.nipol,0.)
        sig_dat_bin2 = np.sqrt(sig_dat_bin2)
        sig_dat_err_bin2 = sig_dat_bin2 / np.sqrt(count_bin2)
        
        z_dat_bin2, nu_dat_bin2, count_bin2 = bincou(z_dat2[index2], zmin, zmax, gp.nipol)
        nu_dat_err_bin2 = nu_dat_bin2 / np.sqrt(count_bin2)
        renorm2 = max(nu_dat_bin2) # normalize by max density of first bin, rather
        nu_dat_bin2 = nu_dat_bin2 / renorm2
        nu_dat_err_bin2 = nu_dat_err_bin2 / renorm2


    if not gp.uselike:
        sel = (z_dat_bin>0)
        xip = z_dat_bin[sel]              # [pc]

        gp.dat.Mx = xip                   # [pc]
        gp.dat.Mdat = K*xip/np.sqrt(xip**2.+D**2.) / (2.0*np.pi*gp.G1)
        gp.dat.Merr      = gp.dat.Mdat*nu_dat_err_bin/nu_dat_bin

        gp.Mmodel = (K*xip/np.sqrt(xip**2.+D**2.)+2.*F*xip) / (2.0*np.pi*gp.G1)

        
        Kz_zstar = -(K*xip/np.sqrt(xip**2.+D**2.) + 2.0 * F * xip)
        if gp.adddarkdisc:
            DD = 0.6
            KD = 0.15 * 1.650
            Kz_zstar = Kz_zstar - KD*zstar/np.sqrt(xip**2.+DD**2.)   # [TODO]
        
        gp.dat.densx     = xip          # TODO: needed somewhere?
        gp.dat.densdat   = -Kz_zstar/(2.*np.pi*gp.G1) # TODO: stellar surface density, store elsewhere
        gp.dat.denserr   = gp.dat.densdat * nu_dat_err_bin/nu_dat_bin

        gp.dat.nux1 = xip
        gp.dat.nudat1 = nu_dat_bin[sel]
        gp.dat.nuerr1 = nu_dat_err_bin[sel]
        
        gp.dat.sigx1 = xip
        gp.dat.sigdat1 = sig_dat_bin[sel]
        gp.dat.sigerr1 = sig_dat_err_bin[sel]

        if gp.pops == 2:
            gp.dat.nux2 = xip
            gp.dat.nudat2 = nu_dat_bin2[sel]
            gp.dat.nuerr2 = nu_dat_err_bin2[sel]

            gp.dat.sigx2 = xip
            gp.dat.sigdat2 = sig_dat_bin2[sel]
            gp.dat.sigerr2 = sig_dat_err_bin2[sel]

        # gp.dat.output()
    gp.dat.save(gp.files.dir+'pp') # pickle
    return gp.dat


## calculate K_z parameters (see Justin Read's paper for definition)
def get_kzpars():
    Kz_zthd = -2.0 * F * zth              # [with pc]
    Sigz_zth = abs(Kz_zthd) / (2.0*np.pi*gp.G1)
    denth = gh.derivcoarse(Sigz_zth,zth)
    kzpars = gh.ipol(zth,abs(denth),gp.xipol)*(2.0*np.pi*gp.G1)
    return kzpars
