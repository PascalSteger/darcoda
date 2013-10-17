#!/usr/bin/env ipython-python3.2
# (c) 2013 Pascal S.P. Steger
'''calculate surface mass density falloff of circular rings around center of mass'''
'''calculate velocity dispersion of 2D rings from a Walker dataset'''
'''calculate 4th order velocity moment (kurtosis) of 2D rings from a Walker dataset'''
# do this consistently with always the same sets of particles per bin
# 3D version, see grw_MCMCbin for 2D version only

import sys
import pdb
import math
import numpy as np
from scipy.stats import kurtosis
from pylab import *

import gl_params as gp
import gr_params as gpr
import gl_file as gfile
from gl_helper import expDtofloat, bin_r_linear, bin_r_log, bin_r_const_tracers
from gl_class_files import *
from BiWeight import meanbiweight



def run():
    xall,yall = np.loadtxt(gpr.get_com_file(0),skiprows=1,usecols=(0,1),unpack=True) # 2*[Rcore]
    # calculate 2D radius on the skyplane
    R = np.sqrt(xall**2+yall**2) #[Rcore]
    # set number and size of (linearly spaced) bins
    Rmin = 0. # [Rcore]
    Rmax = max(r) if gpr.rprior<0 else 1.0*gpr.rprior # [Rcore]
    print('Rmax [Rcore] = ', Rmax)
    R = R[(R<Rmax)]

    # determine radius once and for all
    # this must not be changed between readout and gravlite run
    # if you wish to change: set gp.getnewdata = True in gl_params.py
    if gp.lograd:
        print(gpr.nbins,' bins in log spacings')
        Binmin, Binmax, Rbin = bin_r_log(Rmax/gpr.nbins, Rmax, gpr.nbins)
    elif gp.consttr:
        print(len(R)/gpr.nbins,' particles per bin')
        Binmin, Binmax, Rbin = bin_r_const_tracers(R, len(R)/gpr.nbins)
    else:
        print(gpr.nbins, ' bins in linear spacings')
        Binmin, Binmax, Rbin = bin_r_linear(Rmin, Rmax, gpr.nbins)


    # volume of a circular ring from binmin to binmax
    vol = np.zeros(gpr.nbins)
    for k in range(gpr.nbins):
        vol[k] = 4.*np.pi/3.*(Binmax[k]**3 - Binmin[k]**3) # [Rcore^3]


    for comp in range(gpr.ncomp):
        print('#######  working on component ',comp)
        print('input: ',gpr.get_com_file(comp)+'_3D')
        # start from data centered on COM already:
        if gfile.bufcount(gpr.get_com_file(comp)+'_3D')<2: continue
        x,y,z,v = np.loadtxt(gpr.get_com_file(comp)+'_3D',\
                           skiprows=1,usecols=(0,1,2,3),unpack=True) # 3*[Rcore], [km/s]

        # calculate 2D radius on the skyplane
        r = np.sqrt(x**2+y**2) # [Rcore]
        
        # set maximum radius (if gpr.rprior is set)
        rmax = max(r) if gpr.rprior<0 else 1.0*gpr.rprior # [Rcore]
        print('rmax [Rcore] = ', rmax)
        sel = (r<=rmax)
        x = x[sel]; y = y[sel]; z = z[sel]; v = v[sel]; r = r[sel] # [Rcore]
        totmass = 1.*len(x) # [munit], munit = 1/star
            
        rs = r                   # + possible starting offset, [Rcore]
        vlos = v                 # + possible starting offset, [km/s]
        
        print('output density: ')
        print(gpr.get_ntracer_file(comp)+'_3D')
        tr = open(gpr.get_ntracer_file(comp)+'_3D','w')
        print(totmass, file=tr)
        tr.close()

        print(gpr.get_dens_file(comp)+'_3D')
        de = open(gpr.get_dens_file(comp)+'_3D','w')
        print('rbin','binmin','binmax','nu(r)/nu(0)','error', file=de)

        print(gpr.get_enc_mass_file(comp)+'_3D')
        em = open(gpr.get_enc_mass_file(comp)+'_3D','w')
        print('rbin','binmin','binmax','M(<r)','error', file=em)


        # gpr.n=30 iterations for getting random picked radius values
        density = np.zeros((gpr.nbins,gpr.n))
        a       = np.zeros((gpr.nbins,gpr.n)) # shared by density, siglos, kappa calcs
        for k in range(gpr.n):
            rsi = gpr.Rerror * np.random.randn(len(rs)) + rs # [Rcore]
            vlosi = gpr.vrerror*np.random.randn(len(vlos)) + vlos # [km/s]
            for i in range(gpr.nbins):
                ind1 = np.argwhere(np.logical_and(rsi>=Binmin[i], rsi<Binmax[i])).flatten() # [1]
                density[i][k] = (1.*len(ind1))/vol[i]*totmass # [munit/Rcore^2]
                vlos1 = vlosi[ind1]                           # [km/s]
                a[i][k] = 1.*len(ind1)                        # [1]

        # output density
        dens0 = np.sum(density[0])/(1.*gpr.n) # [munit/Rcore^3]
        print('dens0 = ',dens0,' [munit/Rcore^3]')
        crcore = open(gpr.get_params_file(comp)+'_3D','r')
        Rcore = np.loadtxt(crcore, comments='#', skiprows=1, unpack=False)
        crcore.close()

        cdens = open(gpr.get_params_file(comp)+'_3D','a')
        print(dens0, file=cdens)               # [munit/Rcore^3]
        print(dens0/Rcore**3, file=cdens)      # [munit/pc^3]
        print(totmass, file=cdens)             # [munit]
        cdens.close()

        ab0   = np.sum(a[0])/(1.*gpr.n)     # [1]
        denserr0 = dens0/np.sqrt(ab0)       # [munit/Rcore^3]
        p_dens  = np.zeros(gpr.nbins);  p_edens = np.zeros(gpr.nbins)
        for b in range(gpr.nbins):
            dens = np.sum(density[b])/(1.*gpr.n) # [munit/Rcore^3]
            ab   = np.sum(a[b])/(1.*gpr.n)       # [1]
            denserr = dens/np.sqrt(ab)       # [munit/Rcore^3]
            denserror = np.sqrt((denserr/dens0)**2+(dens*denserr0/(dens0**2))**2) #[1]
            if(math.isnan(denserror)):
                denserror = 0. # [1]
                p_dens[b] = p_dens[b-1]  # [1]
                p_edens[b]= p_edens[b-1] # [1]
            else:
                p_dens[b] = dens/dens0   # [1]
                p_edens[b]= denserror    # [1] #100/rbin would be artificial guess

            print(Rbin[b],Binmin[b],Binmax[b],p_dens[b],p_edens[b], file=de) # [Rcore], 2*[dens0]
            indr = (r<Binmax[b])
            menclosed = 1.0*np.sum(indr)/totmass # for normalization to 1  # [totmass]
            merror = menclosed/np.sqrt(ab) # artificial menclosed/10 # [totmass]
            print(Rbin[b],Binmin[b],Binmax[b],menclosed,merror, file=em) # [rcore], 2*[totmass]
            # TODO: check: take rbinmax for MCMC?
        de.close()
        em.close()


        if not gpr.showplots: continue
        # plot density
        ion(); subplot(111)
        print('rbin = ',Rbin)
        print('p_dens = ',p_dens)
        print('p_edens = ',p_edens)

        plot(Rbin,p_dens,'b',lw=1)
        lbound = p_dens-p_edens; lbound[lbound<1e-6] = 1e-6
        ubound = p_dens+p_edens; 
        fill_between(Rbin,lbound,ubound,alpha=0.5,color='r')
        yscale('log')
        xlim([0,3.])
        ylim([np.min(lbound),np.max(ubound)])
        xlabel(r'$r [r_c]$')
        ylabel(r'$\nu(r)/\nu(0)$')
        savefig(gpr.get_dens_png(i)+'_3D.png')
        if gpr.showplots:
            ioff(); show(); clf()

if __name__ == '__main__':
    gpr.showplots = True
    run()

