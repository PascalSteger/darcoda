#!/usr/bin/python
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
'''generate simple disc data ourselves'''
import gl_params as gp
import gl_helper as gh
import numpy as np
import scipy
import pdb
import scipy.integrate
from scipy.integrate import simps,trapz
import numpy.random as npr
import random
from binsmooth import binsmoo
from bincount import bincou
if gp.geom == 'disc':
    import physics_disc as phys
else:
    import physics_sphere as phys

global K,C,D,F, zth, zp_kz, zmin, zmax, z0, z02

# Set up simple population here using analytic formulae: 
zmin = 100. ; zmax = 1300.; zpnts = 60 #00
zth = 1.* np.arange(zpnts) * (zmax-zmin)/(zpnts-1.) + zmin 
z0 = 240.
K = 1.0 * 1650.
F = 0.1 * 1.650
C = 17.**2.*1000.
D = 250.
z02 = 200.


def disc_simple():
    '''generate disc data from analytic form, return 3D densities, delta, ...'''

    # trick to speed up things considerably, after first data generated, if zpnts > 50 (or other low number)
    # gp.dat.load(gp.files.dir+'pp')
    # return gp.dat

    
    # Set up simple population here using analytic formulae: 
    # zmin = 50. ; zmax = 1000.; zpnts = 50 # 2*[pc], 1 (zpnts = 5000 previously, too many bins)
    # zth = 1.* np.arange(zpnts) * (zmax-zmin)/(zpnts-1.) + zmin 
    # z0 = 240.                             # [pc], scaleheigth
    # K = 1.0 * 1650.                       # [TODO]
    # F = 0.1 * 1.650                       # [TODO]
    # C = 17.**2.*1000.                     # [TODO]
    # D = 250.                              # [pc] TODO
    # z02 = 200.                            # [pc], scaleheight for second population
    
    # Draw mock data: 
    nu_zth = np.exp(-zth/z0)              # [1]
    Kz_zth = -(K*zth/np.sqrt(zth**2.+D**2.) + 2.0 * F * zth)   # [TODO]

    if gp.adddarkdisc: # False
        DD = 600       # [pc]
        KD = 0.15 * 1650.                                # [TODO]
        Kz_zth = Kz_zth - KD*zth/np.sqrt(zth**2.+DD**2.) # [TODO]

    inti = np.zeros(zpnts) 
    for i in range(1,zpnts):
        inti[i] = simps(Kz_zth[:i]*nu_zth[:i],zth[:i])
    sigzth = np.sqrt((inti + C) / nu_zth / 1000.)
    nstars = 7e3                          # [1]
    ran = npr.uniform(size=int(nstars))   # [1]
    zstar = -z0 * np.log(1.0 - ran)       # [pc]
    sigzstar = gh.ipol(zth,sigzth,zstar) # > 0 #selection? no, IDL GT means ">", so perhaps '>' is a shift
    ran2 = npr.normal(size=int(nstars))  # [1]
    vzstar = ran2 * sigzstar             # [TODO]

    # Add second population [thick-disc like]: 
    fac2 = 1.
    if gp.addpoptwo:
        # nu_zth = np.exp(-zth/z0) + fac2 * np.exp(-zth/z02)
        nu_zth2 = np.exp(-zth/z02)
        inti    = np.zeros(zpnts)
        for i in range(1,zpnts):
            inti[i] = simps(Kz_zth[:i]*nu_zth2[:i],zth[:i])
        sigzth2 = np.sqrt(inti / nu_zth2 + C / nu_zth2)
        nstars2 = nstars*fac2             # [1]
        ran = npr.uniform(-1.,1.,nstars2) # [1]
        zstar2 = -z02 * np.log(1.0 - ran) # [pc]
        # zstar = [zstar,zstar2]
        sigzstar2 = gh.ipol(zth,sigzth2,zstar2) > 0   # TODO: > 0?
        ran2 = npr.normal(-1.,1,nstars2)              # [1]
        vzstar2 = ran2 * sigzstar2                    # [TODO]

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
    print nu_dat_bin
    nu_dat_err_bin = nu_dat_bin / np.sqrt(count_bin)
    renorm = max(nu_dat_bin)
    nu_dat_bin = nu_dat_bin / renorm
    nu_dat_err_bin = nu_dat_err_bin / renorm

    if gp.addpoptwo:
        z_dat2  = zstar2[(zstar2 < zmax)]
        vz_dat2 = vzstar2[(zstar2 < zmax)]
        
        # Cut zero velocities: 
        z_dat2 = z_dat2[(abs(vz_dat2) > 0)]
        vz_dat2 = vz_dat2[(abs(vz_dat2) > 0)]
        
        # Calulate binned data (for plots/binned anal.): 
        index2 = np.argsort(z_dat2)
        z_dat_bin2, sig_dat_bin2, count_bin2 = binsmoo(z_dat2[index2],vz_dat2[index2],zmin,zmax,gp.nipol,0.)
        sig_dat_bin2 = np.sqrt(sig_dat_bin2)
        sig_dat_err_bin2 = sig_dat_bin2 / np.sqrt(count_bin2)
        
        z_dat_bin2, nu_dat_bin2, count_bin2 = bincou(z_dat2[index2], zmin, zmax, gp.nipol)
        nu_dat_err_bin2 = nu_dat_bin2 / np.sqrt(count_bin2)
        renorm2 = max(nu_dat_bin2)
        nu_dat_bin2 = nu_dat_bin2 / renorm2
        nu_dat_err_bin2 = nu_dat_err_bin2 / renorm2


    if not gp.uselike:
        sel = (z_dat_bin>0)
        xip = z_dat_bin[sel]              # [pc]
        gp.dat.Mx = xip                   # [pc]
        gp.dat.Mdat = (K*xip/np.sqrt(xip**2.+D**2.)+2.*F*xip) / (2.0*np.pi*gp.G1) / 1000.
        # TODO: xip, also for zstar
        gp.blow = K*xip/np.sqrt(xip**2.+D**2.) / (2.0*np.pi*gp.G1) / 1000.
        # gp.blow = K*xip/np.sqrt(xip**2.+D**2.)
        # TODO: check length

        gp.dat.Merr      = gp.dat.Mdat*0.01
        gp.dat.densx     = xip          # TODO: needed somewhere?

        Kz_zstar = -(K*xip/np.sqrt(xip**2.+D**2.) + 2.0 * F * xip)
        if gp.adddarkdisc:
            DD = 0.6
            KD = 0.15 * 1650.
            Kz_zstar = Kz_zstar - KD*zstar/np.sqrt(xip**2.+DD**2.)   # [TODO]
        
        gp.dat.densdat   = -Kz_zstar/(2.*np.pi*gp.G1) # TODO: surface density, must be stored somewhere else
        gp.dat.denserr   = gp.dat.Merr/gp.dat.Mdat*gp.dat.densdat

        gp.dat.nux1 = xip
        gp.dat.nudat1 = nu_dat_bin[sel]
        gp.dat.nuerr1 = nu_dat_err_bin[sel]
        
        gp.dat.sigx1 = xip
        gp.dat.sigdat1 = sig_dat_bin[sel]
        gp.dat.sigerr1 = sig_dat_err_bin[sel]

        if gp.addpoptwo:
            gp.dat.nux2 = xip
            gp.dat.nudat2 = nu_dat_bin2[sel]
            gp.dat.nuerr2 = nu_dat_err_bin2[sel]

            gp.dat.sigx2 = xip
            gp.dat.sigdat2 = sig_dat_bin2[sel]
            gp.dat.sigerr2 = sig_dat_err_bin2[sel]

        gp.dat.output()
    gp.dat.save(gp.files.dir+'pp') # pickle
    return gp.dat


def get_kzpars():
    Kz_zthd = -2.0 * F * zth              # [with pc]
    Sigz_zth = abs(Kz_zthd) / (2.0*np.pi*gp.G1)
    denth = gh.derivcoarse(Sigz_zth,zth)
    kzpars = gh.ipol(zth,abs(denth),gp.xipol)*(2.0*np.pi*gp.G1)
    return kzpars
