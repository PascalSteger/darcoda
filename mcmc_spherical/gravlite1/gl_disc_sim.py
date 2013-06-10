#!/usr/bin/python
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
'''read data from simulation'''
import gl_params as gp
import gl_helper as gh
import numpy as np
import gl_plot as gpl
from binsmooth import *
from bincount import *
import scipy
import scipy.integrate
import scipy.special as ss
from scipy.integrate import simps,trapz
import pdb
from binsmooth import *

def disc_sim():
    if not gp.investigate == 'sim':
        print 'wrong file included'
        return

    gp.zpmin = -1
    gp.zpmax = -1

    #import all data from files
    if gp.importdata:
      z_nu1_raw,nu1_dat_raw,nu1_dat_err_raw = gh.readcoln(gp.files.nufiles[0])
      z_sig1_raw,sig1_dat_raw,sig1_dat_err_raw = gh.readcoln(gp.files.sigfiles[0])
      #z_surf_raw,surftot_dat_raw,surftot_dat_err_raw = gh.readcoln(gp.files.surfdenfiles[0])
      z_surf_raw,surfdm_dat_raw,surfdm_dat_err_raw = gh.readcoln(gp.files.surfdenfiles[1])
      z_surf_raw,surfbar_dat_raw,surfbar_dat_err_raw = gh.readcoln(gp.files.surfdenfiles[0])

      selnu1 = (z_nu1_raw > 0)
      selsig1 = (z_sig1_raw > 0)
      selsurf = (z_surf_raw > 0)

      if gp.pops == 2:
        z_nu2_raw,nu2_dat_raw,nu2_dat_err_raw = gh.readcoln(gp.files.nufiles[1])
        z_sig2_raw,sig2_dat_raw,sig2_dat_err_raw = gh.readcoln(gp.files.sigfiles[2])

        selnu2 = (z_nu2_raw > 0)
        selsig2 = (z_sig2_raw > 0)  

      #baryonic surface density
      gp.dat.Mx = z_surf_raw[selsurf]*1000.         # [pc]
      gp.dat.Mdat = surfbar_dat_raw[selsurf]          # [Msun/pc^2]
      gp.dat.Merr = surfbar_dat_err_raw[selsurf]  

    
      #total surface density
      gp.Mmodel = surftot_dat_raw[selsusrf]           # [Msun/pc^2]
  
      Kz_zstar = -gp.Mmodel * (2.*np.pi*gp.G1)  
 
      #should be kappa data (not sure whether this is necessary)      
      gp.dat.densx     = z_surf_raw[selsusrf]*1000.         # [pc]
      gp.dat.densdat   = phys.kappa(gp.dat.densx, Kz_zstar)   #should be the total kappa, not sure about x array though
      gp.dat.denserr   = gp.dat.densdat  #not necessary for a certainty

      gp.dat.nux1 = z_nu1_raw[selnu1]*1000.
      gp.dat.nudat1 = nu1_dat_raw[selnu1]
      gp.dat.nuerr1 = nu1_dat_err_raw[selnu1]
        
      gp.dat.sigx1 = z_sig1_raw[selsig1]*1000.
      gp.dat.sigdat1 = sig1_dat_raw[selsig1]
      gp.dat.sigerr1 = sig1_dat_err_raw[selsig1]

      if gp.pops == 2:
        gp.dat.nux2 = z_nu2_raw[selnu2]*1000.
        gp.dat.nudat2 = nu2_dat_raw[selnu2]
        gp.dat.nuerr2 = nu2_dat_err_raw[selnu2]

        gp.dat.sigx2 = z_sig2_raw[selsig2]*1000.
        gp.dat.sigdat2 = z_sig2_raw[selsig2]
        gp.dat.sigerr2 = z_sig2_raw[selsig2]

      gp.dat.output()
      gp.dat.save(gp.files.dir+'pp') # pickle
      return gp.dat  
  

    #import simulation datapoints and calculate nu,sig
    else:
      zmin = 100. ; zmax = 1300.

      # Read in the data:
      mass, x_dat,y_dat,z_dat, vx_dat,vy_dat,vz_dat, pot_dat = gh.readcoln(gp.files.posvelfiles[0])
      if max(mass) != min(mass):
        print '**** Multimass data not yet supported ****'
        exit(1)

      #change to [pc]
      x_dat*=1000.
      y_dat*=1000.
      z_dat*=1000.
      #v is in units [100 km/s]

      z_mean  = np.sum(z_dat)/(1.*len(z_dat))
      vz_mean = np.sum(vz_dat)/(1.*len(z_dat))

      # center on coordinate, if also negative z values read in
      if min(z_dat)<0:
        z_dat = z_dat - z_mean
        vz_dat = vz_dat - vz_mean
    
      # Add errors:
      if gp.adderrors:
        # Assume normal errors for now: 
        xerrfac  = 10.0;  yerrfac = 10.0;  zerrfac = 10.0
        vxerrfac = 10.0; vyerrfac = 10.0; vzerrfac = 10.0
        x_dat_err = abs(x_dat/xerrfac);    y_dat_err = abs(y_dat/yerrfac);    z_dat_err = abs(z_dat/zerrfac)
        vx_dat_err = abs(vx_dat/vxerrfac); vy_dat_err = abs(vy_dat/vyerrfac); vz_dat_err = abs(vz_dat/vzerrfac)
        
        x_dat = x_dat + npr.normal(-1.,1.,len(z_dat)) * x_dat_err
        y_dat = y_dat + npr.normal(-1.,1.,len(z_dat)) * y_dat_err
        z_dat = z_dat + npr.normal(-1.,1.,len(z_dat)) * z_dat_err
        vx_dat = vx_dat + npr.normal(-1.,1.,len(z_dat)) * vx_dat_err
        vy_dat = vy_dat + npr.normal(-1.,1.,len(z_dat)) * vy_dat_err
        vz_dat = vz_dat + npr.normal(-1.,1.,len(z_dat)) * vz_dat_err

      #Cut on zmax
      sel = (z_dat < zmax)
      z_dat = zstar[sel];   vz_dat = vzstar[sel]
    
      # Cut zero velocities:
      sel = (abs(vz_dat) > 0)
      z_dat = z_dat[sel];   vz_dat = vz_dat[sel] 
        
      # Calculate and output binned data:
      zbin = gp.nipol
      index = np.argsort(z_dat)

      # rout, arrayout, count_bin = binsmoo(r, array, low, high, bin, nanflag):

      z_dat_bin,sig_dat_bin,count_bin = binsmoo(z_dat[index],vz_dat[index],zmin,zmax,zbin,0.)
      sig_dat_bin = np.sqrt(sig_dat_bin)
      sig_dat_err_bin = sig_dat_bin / np.sqrt(count_bin)

      # rout, arrrayout, count_bin = bincou(r, low, high, bin):

      z_dat_bin, nu_dat_bin, count_bin = bincou(z_dat[index],zmin,zmax,zbin)
      nu_dat_err_bin = nu_dat_bin / np.sqrt(count_bin)
      renorm = 1.*max(nu_dat_bin)
      nu_dat_bin = nu_dat_bin / renorm
      nu_dat_err_bin = nu_dat_err_bin / renorm

      #[TODO]: downsampling if needed 
    
      if gp.pops == 2:
        # Read in the data:
        mass2, x_dat2,y_dat2,z_dat2, vx_dat2,vy_dat2,vz_dat2, pot_dat2 = gh.readcoln(gp.files.posvelfiles[1])
        if max(mass2) != min(mass2):
          print '**** Multimass data not yet supported ****'
          exit(1)

        #change to [pc]
        x_dat2*=1000.
        y_dat2*=1000.
        z_dat2*=1000.
        #v is in units [100 km/s]

        z_mean2  = np.sum(z_dat2)/(1.*len(z_dat2))
        vz_mean2 = np.sum(vz_dat2)/(1.*len(z_dat2))

        # center on coordinate, if also negative z values read in
        if min(z_dat2)<0:
          z_dat2 = z_dat2 - z_mean2
          vz_dat2 = vz_dat2 - vz_mean2

        # Add errors:
        if gp.adderrors:
          # Assume normal errors for now: 
          xerrfac  = 10.0;  yerrfac = 10.0;  zerrfac = 10.0
          vxerrfac = 10.0; vyerrfac = 10.0; vzerrfac = 10.0
          x_dat_err2 = abs(x_dat2/xerrfac);    y_dat_err2 = abs(y_dat2/yerrfac);    z_dat_err2 = abs(z_dat2/zerrfac)
          vx_dat_err2 = abs(vx_dat2/vxerrfac); vy_dat_err2 = abs(vy_dat2/vyerrfac); vz_dat_err2 = abs(vz_dat2/vzerrfac)
        
          x_dat2 = x_dat2 + npr.normal(-1.,1.,len(z_dat2)) * x_dat_err2
          y_dat2 = y_dat2 + npr.normal(-1.,1.,len(z_dat2)) * y_dat_err2
          z_dat2 = z_dat2 + npr.normal(-1.,1.,len(z_dat2)) * z_dat_err2
          vx_dat2 = vx_dat2 + npr.normal(-1.,1.,len(z_dat2)) * vx_dat_err2
          vy_dat2 = vy_dat2 + npr.normal(-1.,1.,len(z_dat2)) * vy_dat_err2
          vz_dat2 = vz_dat2 + npr.normal(-1.,1.,len(z_dat2)) * vz_dat_err2

        #Cut on zmax
        sel = (z_dat2 < zmax)
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


    
      # if gp.bprior:
      # Load the baryonic model:
      if gp.baryonmodel == 'silvia':
        zvis,sigexpvis,sigexpviserr,sigsecvis,sigsecviserr = gh.readcoln('/home/ast/user/jread/Data/Local_dm/Vis/Sigma_MM.txt')
        sigusevis = sigsecvis
        siguseviserr = sigsecviserr
      elif gp.baryonmodel == 'sim':
        # read in baryonic surface density
        zvis,sigusevis,siguseviserr = gh.readcol(gp.files.surfdenfiles[0])
        # read in DM surface density
        zdm,sigusedm,sigusedmerr = gh.readcol(gp.files.surfdenfiles[1])
      elif gp.baryonmodel == 'simple':
        zth = np.arange(zpnts) * (zmax-zmin)/(zpnts-1.) + zmin
        zvis = zth
        D = 0.25
        sigusevis = K*zvis/sqrt(zvis**2.+D**2.) / (2.0*np.pi*G1) / 1000**2.
        siguseviserr = sigusevis*0.01

      sel = (z_dat_bin>0)
      xip = z_dat_bin[sel]  

      #baryonic surface density
      gp.dat.Mx = zvis*1000.         # [pc]
      gp.dat.Mdat = sigusevis          # [Msun/pc^2]
      gp.dat.Merr = siguseviserr

      #total surface density (same z array as baryonic)
      gp.Mmodel = sigusevis + sigusedm        # [Msun/pc^2]
  
      Kz_zstar = -gp.Mmodel * (2.*np.pi*gp.G1)  
 
      #should be kappa data (not sure whether this is necessary)      
      gp.dat.densx     = zvis*1000.         # [pc]
      gp.dat.densdat   = phys.kappa(gp.dat.densx, Kz_zstar)   #should be the total kappa, not sure about x array though
      gp.dat.denserr   = gp.dat.densdat  #not necessary for a certainty

      gp.dat.nux1 = xip
      gp.dat.nudat1 = nu_dat_bin
      gp.dat.nuerr1 = nu_dat_bin_err
        
      gp.dat.sigx1 = xip
      gp.dat.sigdat1 = sig_dat_bin
      gp.dat.sigerr1 = sig_dat_bin_err

      if gp.pops == 2:
        gp.dat.nux2 = xip
        gp.dat.nudat2 = nu_dat_bin2
        gp.dat.nuerr2 = nu_dat_bin_err2

        gp.dat.sigx2 = xip
        gp.dat.sigdat2 = sig_dat_bin2
        gp.dat.sigerr2 = sig_dat_bin_err2

      gp.dat.output()
      gp.dat.save(gp.files.dir+'pp') # pickle
      return gp.dat 


