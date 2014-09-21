#!/usr/bin/env ipython3

##
# @ file
# generate simple disc data ourselves

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch


import numpy as np
import numpy.random as npr
import pdb
from scipy.integrate import simps,trapz

import gl_helper as gh


def write_disc_output_files(Rbin, Binmin, Binmax, nudat, nuerr, Sigdat, Sigerr, Mdat, Merr, sigdat, sigerr, scales, gp):
    for pop in range(gp.pops+1):
        # write scales
        crscale = open(gp.files.get_scale_file(pop), 'w')
        print('# Zscale in [pc], surfdens_central (=dens0) in [Munit/Zscale^2], and in [Munit/pc^2], and totmass [Munit], and max(v_LOS) in [km/s]', file=crscale)
        for b in range(4):
            # Zscale, Dens0pc, totmass, maxsiglos
            # Dens0pc is *not* density at halflight radius as in spherical
            # code, but rather 3D tracer density at center
            print(scales[pop][b], file=crscale) 
        crscale.close()

        # write tracer densities
        de = open(gp.files.Sigfiles[pop], 'w')
        print('Rbin [Zscale];','Binmin [Zscale];','Binmax [Zscale];',\
              'Sig(R) [Munit/pc^2];','error [1]', file=de)
        for b in range(gp.nipol):
            print(Rbin[b], Binmin[b], Binmax[b], \
                  nudat[pop][b], nuerr[pop][b], file=de)
        de.close()
        
        # write velocity dispersion
        sigfil = open(gp.files.sigfiles[pop],'w')
        print('R [Zscale];','Binmin [Zscale];','Binmax [Zscale];',\
              'sigma_p(R) [maxsiglos];','error [km/s]', file=sigfil)
        for b in range(gp.nipol):
            print(Rbin[b], Binmin[b], Binmax[b], \
                  sigdat[pop][b], sigerr[pop][b], file=sigfil)
        sigfil.close()
        
        # write enclosed mass
        em = open(gp.files.massfiles[pop],'w')
        print('R [Zscale];','Binmin [Zscale];','Binmax [Zscale];',\
              'M(<Binmax) [Munit];','error [Munit]', file=em)
        for b in range(gp.nipol):
            print(Rbin[b], Binmin[b], Binmax[b], \
                  Mdat[pop][b], Merr[pop][b], file=em)
        em.close()
        

    

## \fn write_disc_output_files(Rbin, Binmin, Binmax, nudat, nuerr, Sigdat, Sigerr, Mdat, Merr, sigdat, sigerr, scales, gp)
# for permanent and consistent data handling
# @param Rbin
# @param Binmin
# @param Binmax
# @param nudat
# @param nuerr
# @param Sigdat
# @param Sigerr
# @param Mdat
# @param Merr
# @param sigdat
# @param sigerr
# @param scales
# @param gp global parameters


def run(gp):
    global K,C,D,F, zth, zp_kz, zmin, zmax, z0, z02
    # Set up simple population here using analytic formulae: 
    zmin = 100.                               # [pc], first bin center
    zmax = 1300.                              # [pc], last bin center
    # get Stuetzpunkte for theoretical profiles (not yet stars, finer spacing in real space)
    nth = gp.nipol
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
    ran2 = npr.normal(size=int(gp.ntracer[2-1]))  # [1]
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
    z_dat1  = zstar[sel]
    vz_dat1 = vzstar[sel]
    
    # throw away velocities of value zero (unstable?):
    sel = (abs(vz_dat1) > 0)
    print('fraction of vz_dat>0 selected elements: ', 1.*sum(sel)/(1.*len(sel)))
    z_dat1  = z_dat1[sel]
    vz_dat1 = vz_dat1[sel]
    
    # Calulate binned data (for plots/binned anal.). old way, linear spacings, no const #particles/bin
    binmin1, binmax1, z_dat_bin1, sig_dat_bin1, count_bin1 = gh.binsmooth(z_dat1, vz_dat1, \
                                                                     zmin, zmax, gp.nipol, 0.)
    sig_dat_err_bin1 = np.sqrt(sig_dat_bin1) # Poisson errors
    
    nu_dat_bin1, nu_dat_err_bin1 = gh.bincount(z_dat1, binmax1)
    Mrdat1 = np.cumsum(nu_dat_bin1)
    Mrerr1 = Mrdat1*nu_dat_err_bin1/nu_dat_bin1

    scales=[[],[],[]]
    scales[1].append(z0) # [pc]
    scales[1].append(max(nu_dat_bin1))
    scales[1].append(Mrdat1[-1])
    scales[1].append(max(sig_dat_bin1))

    # start analysis of "all stars" with only component 1,
    # append to it later if more populations required
    z_dat0 = z_dat1   # [pc]
    vz_dat0 = vz_dat1 # [km/s]

    if gp.pops == 2:
        # enforce observational constraints on z<z_max
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
        sig_dat_err_bin2 = np.sqrt(sig_dat_bin2) # Poissonian errors
        
        nu_dat_bin2, nu_dat_err_bin2 = gh.bincount(z_dat2, binmax2)
        Mrdat2 = np.cumsum(nu_dat_bin2)
        Mrerr2 = Mrdat2*nu_dat_err_bin2/nu_dat_bin2

        scales[2].append(z02) # [pc]
        scales[2].append(max(nu_dat_bin2)) # normalize by max density of first bin, rather
        scales[2].append(Mrdat2[-1])
        scales[2].append(max(sig_dat_bin2))
 
        # calculate properties for all pop together with stacked values
        z_dat0 = np.hstack([z_dat1, z_dat2])
        vz_dat0 = np.hstack([vz_dat1, vz_dat2])
        
    # Calulate binned data (for plots/binned anal.). old way, linear spacings, no const #particles/bin
    binmin0, binmax0, z_dat_bin0, sig_dat_bin0, count_bin0 = gh.binsmooth(z_dat0, vz_dat0, \
                                                                          zmin, zmax, gp.nipol, 0.)
    sig_dat_err_bin0 = sig_dat_bin0 / np.sqrt(count_bin0)
    # binmin, binmax, z_dat_bin = gh.bin_r_const_tracers(z_dat, gp.nipol) # TODO: enable, get sig2
    
    nu_dat_bin0, nu_dat_err_bin0 = gh.bincount(z_dat0, binmax0)
    renorm0 = max(nu_dat_bin0)
 
    xip = np.copy(z_dat_bin0)                        # [pc]
    Mrdat0   = K*xip/np.sqrt(xip**2.+D**2.) / (2.0*np.pi*gp.G1)
    Mrerr0   = Mrdat0*nu_dat_err_bin0/nu_dat_bin0

    scales[0].append(D)  # [pc]
    scales[0].append(max(nu_dat_bin0))
    scales[0].append(Mrdat0[-1])
    scales[0].append(max(sig_dat_bin0))

    rmin = binmin0/scales[0][0] # [pc]
    rbin = xip/scales[0][0]     # [pc]
    rmax = binmax0/scales[0][0] # [pc]
    
    # store parameters for output
    # normalized by scale values
    nudat = []
    nudat.append(nu_dat_bin0/scales[0][1])       # [Msun/pc^3]
    nudat.append(nu_dat_bin1/scales[1][1])
    nudat.append(nu_dat_bin2/scales[2][1])

    nuerr = []
    nuerr.append(nu_dat_err_bin0/scales[0][1])   # [Msun/pc^3]
    nuerr.append(nu_dat_err_bin1/scales[1][1])
    nuerr.append(nu_dat_err_bin2/scales[2][1])

    Sigdat = []
    # TODO
    #Sigdat0 = gl_project.nu_param_INT_Sig_disc(xip, nu_dat_bin0, 0, gp)
    #Sigdat1 = gl_project.nu_param_INT_Sig_disc(xip, nu_dat_bin1, 1, gp)
    #Sigdat2 = gl_project.nu_param_INT_Sig_disc(xip, nu_dat_bin2, 2, gp)
    #Sigdat.append(Sigdat0)
    #Sigdat.append(Sigdat1)
    #Sigdat.append(Sigdat2)

    Sigerr = []
    #Sigerr.append(Sigdat0/np.sqrt(len(Sigdat0)))
    #Sigerr.append(Sigdat1/np.sqrt(len(Sigdat1)))
    #Sigerr.append(Sigdat1/np.sqrt(len(Sigdat1)))

    Mrdat = []
    Mrdat.append(Mrdat0/scales[0][2])            # [Msun]
    Mrdat.append(Mrdat1/scales[1][2])
    Mrdat.append(Mrdat2/scales[2][2])
    
    Mrerr = []
    Mrerr.append(Mrerr0/scales[0][2])            # [Msun]
    Mrerr.append(Mrerr1/scales[1][2])
    Mrerr.append(Mrerr2/scales[2][2])

    sigdat = []
    sigdat.append(sig_dat_bin0/scales[0][3])     # [km/s]
    sigdat.append(sig_dat_bin1/scales[1][3])     
    sigdat.append(sig_dat_bin2/scales[2][3])

    sigerr = []
    sigerr.append(sig_dat_err_bin0/scales[0][3]) # [km/s]
    sigerr.append(sig_dat_err_bin1/scales[1][3])
    sigerr.append(sig_dat_err_bin2/scales[2][3])

    write_disc_output_files(rbin, rmin, rmax, nudat, nuerr, \
                            Sigdat, Sigerr, Mrdat, Mrerr,\
                            sigdat, sigerr, scales, gp)
    return gp.dat

## \fn run(gp)
# generate disc data from analytic form, return 3D densities, delta (=tilt)
# write it to files
# @param gp global parameters

