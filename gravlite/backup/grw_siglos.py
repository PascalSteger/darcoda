#!/usr/bin/python
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
'''calculate velocity dispersion of 2D rings from a Walker dataset'''

from pylab import *
import numpy as np
import sys
import math
from BiWeight import meanbiweight
import gl_params as gp
import gr_params as gpr
import gl_file as gfile
from gl_class_files import *
from gl_helper import bin_r_linear, bin_r_log, bin_r_const_tracers



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
        print 'comp = ',comp
        print 'input: ',gpr.get_com_file(comp)
        # start from data centered on COM already:
        if gfile.bufcount(gpr.get_com_file(comp))<2: continue
        x,y,vlos = np.loadtxt(gpr.get_com_file(comp),\
                              skiprows=1,usecols=(0,1,2),unpack=True) #[rcore], [rcore], [km/s]

        # calculate 2D radius on the skyplane
        r = np.sqrt(x**2+y**2) #[rcore]
        
        # set maximum radius (if gpr.rprior is set)
        rmax = max(r) if gpr.rprior<0 else 1.0*gpr.rprior #[rcore]
        print 'rmax [rcore] = ', rmax
        sel = (r<=rmax)
        x = x[sel]; y = y[sel]; vlos = vlos[sel]; r = r[sel] #[rcore]
        totmass = 1.*len(x) #[munit], munit = 1/star

        rs = r
        # no offset from the start!
        # rs = gpr.rerror*np.random.randn(len(r))+r #[rcore]
        # vlos = gpr.vrerror*np.random.randn(len(vlos))+vlos #[km/s]



        print 'output: ',gpr.get_siglos_file(comp)
        vfil = open(gpr.get_siglos_file(comp),'w')
        print >> vfil,'r','sigma_r(r)','error'

        # 30 iterations for drawing a given radius in bin
        dispvelocity = np.zeros((gpr.nbins,gpr.n))
        a = np.zeros((gpr.nbins,gpr.n))
        p_dvlos = np.zeros(gpr.nbins)
        p_edvlos = np.zeros(gpr.nbins)

        for k in range(gpr.n):
            rsi = gpr.rerror*np.random.randn(len(rs))+rs #[rcore]
            vlosi = gpr.vrerror*np.random.randn(len(vlos))+vlos #[km/s]
            for i in range(gpr.nbins):
                ind1 = np.argwhere(np.logical_and(rsi>binmin[i],rsi<binmax[i])).flatten()
                a[i][k] = len(ind1) #[1]
                vlos1 = vlosi[ind1] #[km/s]
                if(len(ind1)<=1):
                    dispvelocity[i][k] = dispvelocity[i-1][k]
                    # attention! should be 0, uses last value
                else:
                    dispvelocity[i][k] = meanbiweight(vlos1,ci_perc=68.4,\
                                                      ci_mean=True,ci_std=True)[1]
                    # [km/s], see BiWeight.py

        for i in range(gpr.nbins):
            dispvel = np.sum(dispvelocity[i])/gpr.n #[km/s]
            ab = np.sum(a[i])/(1.*gpr.n) #[1]
            if ab == 0:
                dispvelerror = p_edvlos[i-1] #[km/s]
                # attention! uses last error
            else:
                dispvelerror = dispvel/np.sqrt(ab) #[km/s]
            p_dvlos[i] = dispvel      #[km/s]
            p_edvlos[i]= dispvelerror #[km/s]

        maxvlos = max(p_dvlos) #[km/s]
        print 'maxvlos = ',maxvlos,'[km/s]'
        fpars = open(gpr.get_params_file(comp),'a')
        print >> fpars,maxvlos          #[km/s]
        fpars.close()

        for i in range(gpr.nbins):
            #             [rcore]  [maxvlos]                  [maxvlos]
            print >> vfil,rbin[i], np.abs(p_dvlos[i]/maxvlos),np.abs(p_edvlos[i]/maxvlos) #/np.sqrt(n))
        vfil.close()

        if not gp.testplot_read: continue

        ion(); subplot(111)
        print 'rbin = ',rbin,' rcore'
        print 'p_dvlos = ',p_dvlos,' km/s'
        print 'p_edvlos = ',p_edvlos, 'km/s'
        plot(rbin,p_dvlos,'b',lw=1)
        fill_between(rbin,p_dvlos-p_edvlos,p_dvlos+p_edvlos,alpha=0.5,color='r') #[rcore],[km/s],[km/s]

        xlabel(r'$r [\mathrm{rcore}]$')
        ylabel(r'$\langle\sigma_{\mathrm{LOS}}\rangle [\mathrm{km/s}]$')
        ylim([-5,30])
        # xscale('log')
        xlim([np.min(rbin),np.max(rbin)])
        #plt.legend(['\rho','\rho'],'lower left'); #title(dwarf)
        savefig(gpr.get_siglos_png(comp))
        if gpr.showplots:
            ioff();show();clf()


if __name__ == '__main__':
    gp.testplot_read = True
    run()

