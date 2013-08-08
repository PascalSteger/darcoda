#!/usr/bin/python
# (c) 2013 Pascal S.P. Steger
'''calculate surface mass density falloff of circular rings around center of mass'''
'''calculate velocity dispersion of 2D rings from a Walker dataset'''
'''calculate 4th order velocity moment (kurtosis) of 2D rings from a Walker dataset'''
# do this consistently with always the same sets of particles per bin

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
    xall,yall = np.loadtxt(gpr.get_com_file(0),skiprows=1,usecols=(0,1),unpack=True) # 2*[rcore]
    # calculate 2D radius on the skyplane
    r = np.sqrt(xall**2+yall**2) #[rcore]
    # set number and size of (linearly spaced) bins
    rmin = 0. #[rcore]
    rmax = max(r) if gpr.rprior<0 else 1.0*gpr.rprior #[rcore]
    print 'rmax [rcore] = ', rmax
    r = r[(r<rmax)]

    # determine radius once and for all
    # this must not be changed between readout and gravlite run
    # if you wish to change: set gp.getnewdata = True in gl_params.py
    if gp.lograd:
        print gpr.nbins,' bins in log spacings'
        binmin, binmax, rbin = bin_r_log(rmax/gpr.nbins, rmax, gpr.nbins)
    elif gp.consttr:
        print len(r)/gpr.nbins,' particles per bin'
        binmin, binmax, rbin = bin_r_const_tracers(r, len(r)/gpr.nbins)
    else:
        print gpr.nbins, ' bins in linear spacings'
        binmin, binmax, rbin = bin_r_linear(rmin, rmax, gpr.nbins)


    # volume of a circular ring from binmin to binmax
    vol = np.zeros(gpr.nbins)
    for k in range(gpr.nbins):
        vol[k] = np.pi*(binmax[k]**2-binmin[k]**2) # [rcore^2]


    for comp in range(gpr.ncomp):
        print '#######  working on component ',comp
        print 'input: ',gpr.get_com_file(comp)
        # start from data centered on COM already:
        if gfile.bufcount(gpr.get_com_file(comp))<2: continue
        x,y,v = np.loadtxt(gpr.get_com_file(comp),\
                           skiprows=1,usecols=(0,1,2),unpack=True) #[rcore], [rcore], [km/s]

        # calculate 2D radius on the skyplane
        r = np.sqrt(x**2+y**2) #[rcore]
        
        # set maximum radius (if gpr.rprior is set)
        rmax = max(r) if gpr.rprior<0 else 1.0*gpr.rprior #[rcore]
        print 'rmax [rcore] = ', rmax
        sel = (r<=rmax)
        x = x[sel]; y = y[sel]; v = v[sel]; r = r[sel] #[rcore]
        totmass = 1.*len(x) # [munit], munit = 1/star
            
        rs = r                   # + possible starting offset, [rcore]
        vlos = v                 # + possible starting offset, [km/s]
        
        print 'output density: '
        print gpr.get_ntracer_file(comp)
        tr = open(gpr.get_ntracer_file(comp),'w')
        print >> tr,totmass
        tr.close()

        print gpr.get_dens_file(comp)
        de = open(gpr.get_dens_file(comp),'w')
        print >> de,'r','nu(r)/nu(0)','error'

        print gpr.get_enc_mass_file(comp)
        em = open(gpr.get_enc_mass_file(comp),'w')
        print >> em,'r','M(<r)','error'


        print 'output siglos: ',gpr.get_siglos_file(comp)
        sigfil = open(gpr.get_siglos_file(comp),'w')
        print >> sigfil,'r','sigma_r(r)','error'


        print 'output kurtosis: ',gpr.get_kurtosis_file(comp)
        kappafil = open(gpr.get_kurtosis_file(comp),'w')
        print >> kappafil,'r','kappa_los(r)','error'


        # gpr.n=30 iterations for getting random picked radius values
        density = np.zeros((gpr.nbins,gpr.n))
        dispvelocity = np.zeros((gpr.nbins,gpr.n))
        mom4         = np.zeros((gpr.nbins,gpr.n))
        a            = np.zeros((gpr.nbins,gpr.n)) # shared by density, siglos, kappa calcs
        for k in range(gpr.n):
            rsi = gpr.rerror * np.random.randn(len(rs)) + rs # [rcore]
            vlosi = gpr.vrerror*np.random.randn(len(vlos)) + vlos #[km/s]
            for i in range(gpr.nbins):
                ind1 = np.argwhere(np.logical_and(rsi>=binmin[i],rsi<binmax[i])).flatten() # [1]
                density[i][k] = (1.*len(ind1))/vol[i]*totmass # [munit/rcore**2]
                vlos1 = vlosi[ind1] #[km/s]

                if(len(ind1)<=1):
                    dispvelocity[i][k] = dispvelocity[i-1][k]
                    mom4[i][k] = mom4[i-1][k]
                    # attention! should be 0, uses last value
                else:
                    dispvelocity[i][k] = meanbiweight(vlos1,ci_perc=68.4,ci_mean=True,ci_std=True)[1]
                                        # [km/s], see BiWeight.py
                    mom4[i][k] = kurtosis(vlos1, axis=0, fisher=False, bias=False) # [1]

                a[i][k] = 1.*len(ind1) #[1]

        # output density
        dens0 = np.sum(density[0])/(1.*gpr.n) # [munit/rcore**2]
        print 'dens0 = ',dens0,'[munit/rcore**2]'
        crcore = open(gpr.get_params_file(comp),'r')
        rcore = np.loadtxt(crcore, comments='#', skiprows=1, unpack=False)
        crcore.close()

        cdens = open(gpr.get_params_file(comp),'a')
        print >> cdens, dens0               # [munit/rcore**2]
        print >> cdens, dens0/rcore**2      # [munit/pc**2]
        print >> cdens, totmass             # [munit]
        cdens.close()

        ab0   = np.sum(a[0])/(1.*gpr.n)     # [1]
        denserr0 = dens0/np.sqrt(ab0)       # [munit/rcore**2]
        p_dens  = np.zeros(gpr.nbins);  p_edens = np.zeros(gpr.nbins)
        for b in range(gpr.nbins):
            dens = np.sum(density[b])/(1.*gpr.n) # [munit/rcore**2]
            ab   = np.sum(a[b])/(1.*gpr.n)       # [1]
            denserr = dens/np.sqrt(ab)       # [munit/rcore**2]
            denserror = np.sqrt((denserr/dens0)**2+(dens*denserr0/(dens0**2))**2) #[1]
            if(math.isnan(denserror)):
                denserror = 0. # [1]
                p_dens[b] = p_dens[b-1]  # [1]
                p_edens[b]= p_edens[b-1] # [1]
            else:
                p_dens[b] = dens/dens0   # [1]
                p_edens[b]= denserror    # [1] #100/rbin would be artificial guess

            print >> de,rbin[b],p_dens[b],p_edens[b] # [rcore], [dens0], [dens0]
            indr = (r<binmax[b])
            menclosed = 1.0*np.sum(indr)/totmass # for normalization to 1  #[totmass]
            merror = menclosed/np.sqrt(ab) # artificial menclosed/10 #[totmass]
            print >> em,rbin[b],menclosed,merror # [rcore], [totmass], [totmass]
            # TODO: check: take rbinmax for MCMC?
        de.close()
        em.close()


        # output siglos
        p_dvlos = np.zeros(gpr.nbins);        p_edvlos = np.zeros(gpr.nbins)
        for b in range(gpr.nbins):
            dispvel = np.sum(dispvelocity[b])/gpr.n #[km/s]
            ab = np.sum(a[b])/(1.*gpr.n) #[1]
            if ab == 0:
                dispvelerror = p_edvlos[b-1] #[km/s]
                # attention! uses last error
            else:
                dispvelerror = dispvel/np.sqrt(ab) #[km/s]
            p_dvlos[b] = dispvel      #[km/s]
            p_edvlos[b]= dispvelerror #[km/s]

        maxvlos = max(p_dvlos) #[km/s]
        print 'maxvlos = ',maxvlos,'[km/s]'
        fpars = open(gpr.get_params_file(comp),'a')
        print >> fpars,maxvlos          #[km/s]
        fpars.close()
        
        for b in range(gpr.nbins):
            #             [rcore]  [maxvlos]                  [maxvlos]
            print >> sigfil,rbin[b], np.abs(p_dvlos[b]/maxvlos),np.abs(p_edvlos[b]/maxvlos)
            # TODO: check uncommented /np.sqrt(n))
        sigfil.close()



        # output kurtosis kappa
        p_kappa = np.zeros(gpr.nbins) # needed for plotting later
        p_ekappa = np.zeros(gpr.nbins)
        for b in range(gpr.nbins):
            kappavel = np.sum(mom4[b])/gpr.n #[1]
            ab = np.sum(a[b])/(1.*gpr.n) #[1]
            if ab == 0:
                kappavelerror = p_edvlos[b-1] #[1]
                # attention! uses last error
            else:
                kappavelerror = np.abs(kappavel/np.sqrt(ab)) #[1]
            p_kappa[b] = kappavel
            p_ekappa[b] = kappavelerror
            
            print >> kappafil,rbin[b], kappavel, kappavelerror # [rcore], 2*[1]
            # TODO: /np.sqrt(n))
        kappafil.close()


    


        if not gp.showplot_readout: continue
        # plot density
        ion(); subplot(111)
        print 'rbin = ',rbin
        print 'p_dens = ',p_dens
        print 'p_edens = ',p_edens

        plot(rbin,p_dens,'b',lw=1)
        lbound = p_dens-p_edens; lbound[lbound<1e-6] = 1e-6
        ubound = p_dens+p_edens; 
        fill_between(rbin,lbound,ubound,alpha=0.5,color='r')
        yscale('log')
        xlim([0,3.])
        ylim([np.min(lbound),np.max(ubound)])
        xlabel(r'$r [r_c]$')
        ylabel(r'$\nu(r)/\nu(0)$')
        savefig(gpr.get_dens_png(i))
        if gpr.showplots:
            ioff(); show(); clf()

        # plot siglos
        ion(); subplot(111)
        print 'rbin = ',rbin,' rcore'
        print 'p_dvlos = ',p_dvlos,' km/s'
        print 'p_edvlos = ',p_edvlos, 'km/s'
        plot(rbin,p_dvlos,'b',lw=1)
        fill_between(rbin,p_dvlos-p_edvlos,p_dvlos+p_edvlos,alpha=0.5,color='r') #[rcore],2*[km/s]

        xlabel(r'$r [\mathrm{rcore}]$')
        ylabel(r'$\langle\sigma_{\mathrm{LOS}}\rangle [\mathrm{km/s}]$')
        ylim([-5,30])
        xlim([0,3])
        savefig(gpr.get_siglos_png(comp))
        if gpr.showplots:
            ioff();show();clf()


        # plot kappa
        ion(); subplot(111)
        print 'rbin = ',rbin,' rcore'
        print 'p_kappa = ',p_kappa
        print 'p_ekappa = ',p_ekappa
        plot(rbin,p_kappa,'b',lw=1)
        fill_between(rbin,p_kappa-p_ekappa,p_kappa+p_ekappa,alpha=0.5,color='r') #[rcore],2*[1]
        xlabel(r'$r [\mathrm{rcore}]$')
        ylabel(r'$\langle\kappa_{\mathrm{LOS}}\rangle [1]$')
        ylim([-3,3])
        xlim([0,3])
        savefig(gpr.get_kurtosis_png(comp))
        if gpr.showplots:
            ioff();show();clf()


if __name__ == '__main__':
    gp.showplot_readout = True
    run()

