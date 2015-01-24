#!/usr/bin/env ipython3

##
# @file
# read data from simulation

# (c) GPL v3 2014 Pascal Steger, pascal@steger.aero

import numpy as np
import numpy.random as npr
import pdb

import gi_units as gu
import gi_helper as gh
import physics_disc as phys

import binsmooth as bs
import bincount as bc


def disc_sim(gp):
    gp.zpmin = -1
    gp.zpmax = -1

    #import all data from files
    if gp.importdata:
        z_nu1_raw,nu1_dat_raw,nu1_dat_err_raw = gh.readcoln(gp.files.Sigfiles[0])
        z_sig1_raw,sig1_dat_raw,sig1_dat_err_raw = gh.readcoln(gp.files.sigfiles[0])
        #z_surf_raw,surftot_dat_raw,surftot_dat_err_raw = gh.readcoln(gp.files.surfdenfiles[0])
        z_surf_raw,surfbar_dat_raw,surfbar_dat_err_raw = gh.readcoln(gp.files.surfdenfiles[0])
        z_surf_raw,surfdm_dat_raw,surfdm_dat_err_raw = gh.readcoln(gp.files.surfdenfiles[1])

        selnu1  = (z_nu1_raw > 0)
        selsig1 = (z_sig1_raw > 0)
        selsurf = (z_surf_raw > 0)

        if gp.pops == 2:
            z_nu2_raw,nu2_dat_raw,nu2_dat_err_raw = gh.readcoln(gp.files.Sigfiles[1])
            z_sig2_raw,sig2_dat_raw,sig2_dat_err_raw = gh.readcoln(gp.files.sigfiles[2])

            selnu2 = (z_nu2_raw > 0)
            selsig2 = (z_sig2_raw > 0)

        # baryonic surface density
        gp.dat.Mx   = z_surf_raw[selsurf]*1000.      # [pc]
        gp.dat.Mrdat = surfbar_dat_raw[selsurf]      # [Munit/pc^2]
        gp.dat.Mrerr = surfbar_dat_err_raw[selsurf]  # [Munit/pc^2]

        # total surface density
        gp.Mmodel = surftot_dat_raw[selsusrf]        # [Munit/pc^2]
        Kz_zstar = -gp.Mmodel * (2.*np.pi*gu.G1__pcMsun_1km2s_2)

        # should be kappa data (not sure whether this is necessary)
        gp.dat.densx     = z_surf_raw[selsusrf]*1000.         # [pc]
        gp.dat.densdat   = phys.kappa(gp.dat.densx, Kz_zstar)   #should be the total kappa, not sure about x array though
        gp.dat.denserr   = gp.dat.densdat  #not necessary for a certainty

        gp.dat.nux1 = z_nu1_raw[selnu1]*1000.   # [pc]
        gp.dat.nu1 = nu1_dat_raw[selnu1]
        gp.dat.nuerr1 = nu1_dat_err_raw[selnu1]

        gp.dat.sigx1 = z_sig1_raw[selsig1]*1000.
        gp.dat.sig1 = sig1_dat_raw[selsig1]
        gp.dat.sigerr1 = sig1_dat_err_raw[selsig1]

        if gp.pops == 2:
            gp.dat.nux2 = z_nu2_raw[selnu2]*1000.
            gp.dat.nu2 = nu2_dat_raw[selnu2]
            gp.dat.nuerr2 = nu2_dat_err_raw[selnu2]

            gp.dat.sigx2 = z_sig2_raw[selsig2]*1000.
            gp.dat.sig2 = z_sig2_raw[selsig2]
            gp.dat.sigerr2 = z_sig2_raw[selsig2]

        gp.dat.output()
        return gp.dat



    else:
        # import simulation datapoints and calculate nu, sig
        zmin = 100. ; zmax = 1300.          # [pc]
        zbinbndry = np.linspace(zmin, zmax, gp.nipol+1)   # [pc] assuming linear spacing of bins
        zbinmin = zbinbndry[:-1]                          # [pc]
        zbinmax = zbinbndry[1:]                           # [pc]
        gp.xipol= zbinmin + (zbinmax-zbinmin)/2.          # [pc]

        # Read in the data:
        mass, x_dat,y_dat,z_dat, vx_dat,vy_dat,vz_dat, pot_dat = gh.readcoln(gp.files.posvelfiles[0])
        # assume units: Munit, 3*kpc, 3*km/s, [pot] <= last one not needed
        # [Dave] v is in units [100 km/s] <= not possible?!
        if max(mass) != min(mass):
            print('**** Multimass data not yet supported ****')
            exit(1)

        # change to [pc]
        x_dat *= 1000.;    y_dat *= 1000.;        z_dat *= 1000. # [pc]
        z_mean  = np.sum(mass*z_dat)/np.sum(mass)                  # [pc]
        vz_mean = np.sum(mass*vz_dat)/np.sum(mass)                 # [km/s]

        # center on coordinate, if also negative z values read in
        if min(z_dat)<0:
            z_dat = z_dat - z_mean        # [pc]
            vz_dat = vz_dat - vz_mean     # [km/s]

        # Add errors:
        if gp.adderrors:
            # Assume normal errors for now:
            xerrfac  = 10.0;  yerrfac = 10.0;  zerrfac = 10.0
            vxerrfac = 10.0; vyerrfac = 10.0; vzerrfac = 10.0
            x_dat_err = abs(x_dat/xerrfac)          # [pc]
            y_dat_err = abs(y_dat/yerrfac)          # [pc]
            z_dat_err = abs(z_dat/zerrfac)          # [pc]
            vx_dat_err = abs(vx_dat/vxerrfac)       # [km/s]
            vy_dat_err = abs(vy_dat/vyerrfac)       # [km/s]
            vz_dat_err = abs(vz_dat/vzerrfac)       # [km/s]

            x_dat = x_dat + npr.normal(-1.,1.,len(z_dat)) * x_dat_err # [pc]
            y_dat = y_dat + npr.normal(-1.,1.,len(z_dat)) * y_dat_err # [pc]
            z_dat = z_dat + npr.normal(-1.,1.,len(z_dat)) * z_dat_err # [pc]
            vx_dat = vx_dat + npr.normal(-1.,1.,len(z_dat)) * vx_dat_err # [km/s]
            vy_dat = vy_dat + npr.normal(-1.,1.,len(z_dat)) * vy_dat_err # [km/s]
            vz_dat = vz_dat + npr.normal(-1.,1.,len(z_dat)) * vz_dat_err # [km/s]

        # Cut on zmax, cut zero velocities
        sel = (z_dat < zmax) * (abs(vz_dat) >= 0.)   # [bool]
        z_dat  = z_dat[sel]               # [pc]
        vz_dat = vz_dat[sel]              # [km/s]

        # determine sigma_v
        sig_dat_bin = np.zeros(gp.nipol)
        sig_dat_err_bin = np.zeros(gp.nipol)
        for i in range(gp.nipol):
            sel = (z_dat > zbinmin[i])*(z_dat < zbinmax[i])   # select bin
            vtemp = np.array(vz_dat[sel])                     # [km/s]
            sig_dat_bin[i] = np.sqrt(np.mean(vtemp**2) - np.mean(vtemp)**2) # [km/s]
            sig_dat_err_bin[i] = sig_dat_bin[i]/(1.*np.sum(sel))   # [km/s]


        nu_dat_bin = np.zeros(gp.nipol)
        nu_dat_err_bin = np.zeros(gp.nipol)
        for i in range(gp.nipol):
            sel = (z_dat > zbinmin[i])*(z_dat < zbinmax[i])   # select bin
            nu_dat_bin[i] = 1.*np.sum(sel)/(1.*(zbinmax[i]-zbinmin[i])) # [1/tot. area/pc]
            nu_dat_err_bin[i] = nu_dat_bin[i] / np.sqrt(np.sum(sel)) # [1/tot. area/pc], Poisson distributed

        renorm = 1.*max(nu_dat_bin)       # [1/tot.area/pc]
        nu_dat_bin = nu_dat_bin / renorm  # [1]
        nu_dat_err_bin = nu_dat_err_bin / renorm   # [1]

        if gp.pops == 2:
            mass2, x_dat2,y_dat2,z_dat2, vx_dat2,vy_dat2,vz_dat2, pot_dat2 = gh.readcoln(gp.files.posvelfiles[1])
            if max(mass2) > min(mass2):
                print('**** Multimass data not yet supported ****')
                exit(1)

            # change to [pc]
            x_dat2 *= 1000.;    y_dat2 *= 1000.;        z_dat2 *= 1000. # [pc]
            z_mean2  = np.sum(mass2*z_dat2)/np.sum(mass2)                  # [pc]
            vz_mean2 = np.sum(mass2*vz_dat2)/np.sum(mass2)                 # [km/s]

            # center on coordinate, if also negative z values read in
            if min(z_dat2)<0:
                z_dat2 = z_dat2 - z_mean2        # [pc]
                vz_dat2 = vz_dat2 - vz_mean2     # [km/s]

            # Add errors:
            if gp.adderrors:
                # Assume normal errors for now:
                x_dat_err2 = abs(x_dat2/xerrfac)          # [pc]
                y_dat_err2 = abs(y_dat2/yerrfac)          # [pc]
                z_dat_err2 = abs(z_dat2/zerrfac)          # [pc]
                vx_dat_err2 = abs(vx_dat2/vxerrfac)       # [km/s]
                vy_dat_err2 = abs(vy_dat2/vyerrfac)       # [km/s]
                vz_dat_err2 = abs(vz_dat2/vzerrfac)       # [km/s]

                x_dat2 = x_dat2 + npr.normal(-1.,1.,len(z_dat2)) * x_dat_err2 # [pc]
                y_dat2 = y_dat2 + npr.normal(-1.,1.,len(z_dat2)) * y_dat_err2 # [pc]
                z_dat2 = z_dat2 + npr.normal(-1.,1.,len(z_dat2)) * z_dat_err2 # [pc]
                vx_dat2 = vx_dat2 + npr.normal(-1.,1.,len(z_dat2)) * vx_dat_err2 # [km/s]
                vy_dat2 = vy_dat2 + npr.normal(-1.,1.,len(z_dat2)) * vy_dat_err2 # [km/s]
                vz_dat2 = vz_dat2 + npr.normal(-1.,1.,len(z_dat2)) * vz_dat_err2 # [km/s]

            # Cut on zmax, cut zero velocities
            sel = (z_dat2 < zmax) * (abs(vz_dat2) >= 0.)   # [bool]
            z_dat2  = z_dat2[sel]               # [pc]
            vz_dat2 = vz_dat2[sel]              # [km/s]

            # determine sigma_v
            sig_dat_bin2 = np.zeros(gp.nipol)
            sig_dat_err_bin2 = np.zeros(gp.nipol)
            for i in range(gp.nipol):
                sel = (z_dat2 > zbinmin[i])*(z_dat2 < zbinmax[i])   # select bin
                vtemp = np.array(vz_dat2[sel])                     # [km/s]
                sig_dat_bin2[i] = np.sqrt(np.mean(vtemp**2) - np.mean(vtemp)**2) # [km/s]
                sig_dat_err_bin2[i] = sig_dat_bin2[i]/(1.*np.sum(sel))   # [km/s]


            nu_dat_bin2 = np.zeros(gp.nipol)
            nu_dat_err_bin2 = np.zeros(gp.nipol)
            for i in range(gp.nipol):
                sel = (z_dat2 > zbinmin[i])*(z_dat2 < zbinmax[i])   # select bin
                nu_dat_bin2[i] = 1.*np.sum(sel)/(1.*(zbinmax[i]-zbinmin[i])) # [1/tot. area/pc]
                nu_dat_err_bin2[i] = nu_dat_bin2[i] / np.sqrt(np.sum(sel)) # [1/tot. area/pc], Poisson distributed

            renorm = 1.*max(nu_dat_bin2)       # [1/tot.area/pc]
            nu_dat_bin2 /= renorm              # [1]
            nu_dat_err_bin2 /= renorm          # [1]


        # if gp.bprior:
        # Load the baryonic model:
        if gp.baryonmodel == 'silvia':
            zvis,sigexpvis,sigexpviserr,sigsecvis,sigsecviserr = gh.readcoln('/home/ast/user/jread/Data/Local_dm/Vis/Sigma_MM.txt')
            # [kpc, Munit/pc^2, Msun/pc^2, Msun/pc^2, Msun/pc^2]
            sigusevis    = sigsecvis      # [Munit/pc^2]
            siguseviserr = sigsecviserr   # [Munit/pc^2]
        elif gp.baryonmodel == 'sim':
            zvis, sigusevis, siguseviserr = gh.readcol3(gp.files.surfdenfiles[0])
            # [kpc, Munit/pc^2, Munit/pc^2]
            zvis *= 1000.                     # [pc]
            sigusevis    = gh.ipol(zvis, sigusevis, gp.xipol)   # interpolate to xipol radius array
            siguseviserr = gh.ipol(zvis, siguseviserr, gp.xipol)
            zvis = gp.xipol                          # [pc]

            # read in DM surface density
            zdm, sigusedm, sigusedmerr = gh.readcol3(gp.files.surfdenfiles[1])
            # [kpc, Munit/pc^2, Munit/pc^2]
            zdm *= 1000.                                # [pc]
            sigusedm = gh.ipol(zdm, sigusedm, gp.xipol) # interpolate to xipol radius array
            sigusedmerr = gh.ipol(zdm, sigusedmerr, gp.xipol)
            zdm = gp.xipol                    # [pc]
        elif gp.baryonmodel == 'simple':
            zvis = gp.xipol                   # [pc]
            D = 250.                          # [pc]
            K = 1.65
            sigusevis = K*zvis/sqrt(zvis**2.+D**2.) / (2.0*np.pi*G1)
            siguseviserr = sigusevis*0.01

        # baryonic surface density, really a Sig
        gp.dat.Mx   = gp.xipol                # [pc]
        gp.dat.Mrdat = sigusevis              # [Munit/pc^2]
        gp.dat.Mrerr = siguseviserr           # [Munit/pc^2]

        # total surface density (same z array as baryonic)
        gp.Mmodel = sigusevis + sigusedm         # [Munit/pc^2]
        Kz_zstar = -gp.Mmodel * (2.*np.pi*gu.G1__pcMsun_1km2s_2) # [1000/pc (km/s)^2]

        # should be kappa data (not sure whether this is necessary)
        gp.dat.densx     = gp.xipol                       # [pc]
        gp.dat.densdat   = phys.kappa(gp.dat.densx, Kz_zstar)
        gp.dat.denserr   = gp.dat.densdat/np.sqrt(len(Kz_zstar))

        gp.dat.nux1   = gp.xipol          # [pc]
        gp.dat.nu1 = nu_dat_bin           # [Munit/pc^3]
        gp.dat.nuerr1 = nu_dat_err_bin    # [Munit/pc^3]

        gp.dat.sigx1   = gp.xipol         # [pc]
        gp.dat.sig1 = sig_dat_bin         # [km/s]
        gp.dat.sigerr1 = sig_dat_err_bin  # [km/s]

        if gp.pops == 2:
            gp.dat.nux2 = gp.xipol               # [pc]
            gp.dat.nu2 = nu_dat_bin2             # [Munit/pc^3]
            gp.dat.nuerr2 = nu_dat_err_bin2      # [Munit/pc^3]

            gp.dat.sigx2 = gp.xipol              # [pc]
            gp.dat.sig2 = sig_dat_bin2           # [km/s]
            gp.dat.sigerr2 = sig_dat_err_bin2    # [km/s]

        gp.dat.output()
        return gp.dat
## \fn disc_sim(gp)
# read in disc simulation file
# @param gp
