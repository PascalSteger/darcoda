#!/usr/bin/python
'''read data from simulation'''
import gl_params as gp
import gl_helper as gh
import numpy as np
import gl_plot as gplot
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

    # Read in the data:
    mass, x_dat,y_dat,z_dat, vx_dat,vy_dat,vz_dat, pot_dat = gh.readcoln(gp.files.posvelfiles[0])
    if max(mass) != min(mass):
        print '**** Multimass data not yet supported ****'
        exit(1)
    
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

    # Slice and dice the data:
    if gp.slicedata:
        # Cut in angle:
        R_dat = np.sqrt(x_dat**2. + y_dat**2.)
        theta = np.atan(y_dat,x_dat)
        theta = theta - np.sum(theta)/(1.*len(theta))
        angle = max(theta) / 4.0
        sel = (abs(theta) < angle)
        x_cut = x_dat[sel];   y_cut = y_dat[sel];   z_cut = z_dat[sel]
        vx_cut = vx_dat[sel]; vy_cut = vy_dat[sel]; vz_cut = vz_dat[sel]
        
        # Cut in z:
        # zmin = min(z_cut)
        zmin = 0.
        sel = (z_cut > zmin)
        x_cut = x_cut[sel];   y_cut = y_cut[sel];   z_cut = z_cut[sel]
        vx_cut = vx_cut[sel]; vy_cut = vy_cut[sel]; vz_cut = vz_cut[sel]
        
        # Cut in R: 
        R_cut = np.sqrt(x_cut**2. + y_cut**2.)
        Rmin = min(R_cut); Rmax = max(R_cut)
        sel = (R_cut > Rmin and R_cut < Rmax)
        x_dat = x_cut[sel];   y_dat = y_cut[sel];   z_dat = z_cut[sel]
        vx_dat = vx_cut[sel]; vy_dat = vy_cut[sel]; vz_dat = vz_cut[sel]
        
        # Downsample:
        nfrac = 1.0
        step = 1.0/nfrac
        nn = 0
        # TODO: random downsampling with npr.sample() possible?
        for jj in range(0,len(z_dat),step):
            x_cut[nn] = x_dat[jj];   y_cut[nn] = y_dat[jj];   z_cut[nn] = z_dat[jj]
            vx_cut[nn] = vx_dat[jj]; vy_cut[nn] = vy_dat[jj]; vz_cut[nn] = vz_dat[jj]
            nn = nn + 1
        
        x_dat = x_cut[:nn-1];   y_dat = y_cut[:nn-1];   z_dat = z_cut[:nn-1]
        vx_dat = vx_cut[:nn-1]; vy_dat = vy_cut[:nn-1]; vz_dat = vz_cut[:nn-1]

    R_dat = np.sqrt(x_dat**2.+y_dat**2.)
    vR_dat = (vx_dat*x_dat+vy_dat*y_dat)/R_dat
    
    if False: #gp.testplot: #TODO: enable if really wanted to see output
        # Test whether dist. are really Gaussian or not:
        import pylab as pl
        # pl.hist(x, bins, normed=1, histtype='bar', rwidth=0.8)
        pl.clf(); histvz,binsvz,pa = pl.hist(vz_dat,  rwidth=2.5)
        pl.clf(); histz,binsz,pa   = pl.hist(z_dat,   rwidth=0.05)
        pl.clf(); histvzvR,binsvzvR,pa = pl.hist(vz_dat*vR_dat,rwidth=25.0)
        
        # Cut a zbin 
        zmin = 0.4;  zmax = 0.6
        sel = (z_dat > zmin) * (z_dat < zmax)
        zcut = z_dat[sel]; vzcut = vz_dat[sel]; vRcut = vR_dat[sel]
        
        pl.clf(); histvzbin,binsvzcut,pa  = pl.hist(vzcut,rwidth=2.5)
        pl.clf(); histzbin,binszcut,pa    = pl.hist(zcut, rwidth=0.05)
        pl.clf(); histvzvRbin,binsvzvRcut,pa = pl.hist(vzcut*vRcut, rwidth=100.0)

        sigvz = 20.
        pl.clf(); gplot.plot(binsvzcut,np.exp(-binsvzcut**2./2./sigvz**2.))
        
        # Assume "normal product" distribution [e.g. Wolfram Mathworld]: 
        sigRz = 2500. 
        meanbfunc = 100.
        tmin = -20001.;        tmax = 20001. ;        tpnts = 10000 
        test = np.arange(tpnts)*(tmax-tmin)/(1.*tpnts) + tmin
        bfunc = ss.kv(0,(abs(test-meanbfunc))/sigRz)
        pl.clf(); gplot.plot(test,bfunc/max(bfunc))

        mean = simps(bfunc*test,test) / simps(bfunc,test)
        
    # Calculate and output binned data:
    zbin = gp.nipol
    index = np.argsort(z_dat)

    # rout, arrayout, count_bin = binsmoo(r, array, low, high, bin, nanflag):

    z_dat_bin,sig_dat_bin,count_bin = binsmoo(z_dat[index],vz_dat[index],-1.2,1.2,zbin,0.)
    sig_dat_bin = np.sqrt(sig_dat_bin)
    sig_dat_err_bin = sig_dat_bin / np.sqrt(count_bin)
    z_dat_bin,sigRz_dat_bin,count_bin = binsmoo(z_dat[index],vz_dat[index]*vR_dat[index],-1.2,1.2,zbin,0.)
    sigRz_dat_err_bin = sigRz_dat_bin / np.sqrt(count_bin)

    # rout, arrrayout, count_bin = bincou(r, low, high, bin):

    z_dat_bin, nu_dat_bin, count_bin = bincou(z_dat[index],-1.2,1.2,zbin)
    nu_dat_err_bin = nu_dat_bin / np.sqrt(count_bin)
    renorm = 1.*max(nu_dat_bin)
    nu_dat_bin = nu_dat_bin / renorm
    nu_dat_err_bin = nu_dat_err_bin / renorm
    
    # z_dat_bin,nu_dat_bin,nu_dat_err_bin = readcol(nufile)
    # z_dat_bin,sig_dat_bin,sig_dat_err_bin = readcol(sigfile)
    sel = (z_dat_bin > 0)
    z_dat = z_dat_bin[sel]

    if not gp.uselike:
        gp.dat.nur1 = z_dat_bin[sel]
        gp.dat.nudat1 = nu_dat_bin[sel]
        gp.dat.nuerr1 = nu_dat_err_bin[sel]
        
        gp.dat.sigr1 = z_dat_bin[sel]
        gp.dat.sigdat1 = sig_dat_bin[sel]
        gp.dat.sigerr1 = sig_dat_err_bin[sel]

        sigRz_dat = sigRz_dat_bin[sel]
        sigRz_dat_err = sigRz_dat_err_bin[sel]

    # Output binned data (now that basename is fully specified): 
    # arraydump(basename+'_nu_dat_bin.txt',[[z_dat_bin],[nu_dat_bin],[nu_dat_err_bin]],3)
    # arraydump(basename+'_sigz_dat_bin.txt',[[z_dat_bin],[sig_dat_bin],[sig_dat_err_bin]],3)
    
    # if gp.bprior:
    # Load the baryonic model:
    if gp.baryonmodel == 'silvia':
        zvis,sigexpvis,sigexpviserr,sigsecvis,sigsecviserr = gh.readcoln('/home/ast/user/jread/Data/Local_dm/Vis/Sigma_MM.txt')
        sigusevis = sigsecvis
        siguseviserr = sigsecviserr
    elif gp.baryonmodel == 'sim':
        zvis,sigusevis,siguseviserr = gh.readcol(gp.files.massfile)
        gp.dat.Mx = zvis
        gp.dat.Mdat = sigusevis
        gp.dat.Merr = siguseviserr
    elif gp.baryonmodel == 'simple':
        zth = np.arange(zpnts) * (zmax-zmin)/(zpnts-1.) + zmin
        zvis = zth
        D = 0.25
        sigusevis = K*zvis/sqrt(zvis**2.+D**2.) / (2.0*np.pi*G1) / 1000**2.
        siguseviserr = sigusevis*0.01
        
    sigvismin = sigusevis - siguseviserr

    # Binning in z:
    if gp.xpmin < 0: gp.xpmin = 1.*min(z_dat)
    if gp.xpmax < 0: gp.xpmax = 1.*max(z_dat)
    if gp.nipol > 0:
        zp_kz = np.arange(gp.nipol) * (gp.xpmax - gp.xpmin) / (gp.nipol-1.) + gp.xpmin
    else:
        # Adaptive auto-binning here:
        zabs = abs(z_dat)
        index = np.sort(zabs)
        zsort = zabs(index)
        numperbin = abs(nbin) 
        ztmp = np.zeros(len(zsort))
        jl = 0.
        jr = 0.
        nbin = 0
        while jr < len(zsort)-1:
            jr = jl + numperbin
            if jr > len(zsort)-1: jr = len(zsort)-1 
            ztmp[nbin] = np.sum(zsort[jl:jr-1])/(1.*jr-jl)
            nbin = nbin + 1
            jl = jl + numperbin
            
        zp_kz = ztmp[:nbin-1]


    # Set up bprior:
    if gp.bprior : 
        blow = gh.ipol(zvis,sigvismin,zp_kz) 
        gp.blow = bsurf * (2.*np.pi*G1) * 1000.**2.
    else:
        gp.blow = np.zeros(len(zp_kz)) + 1e-5

