#!/usr/bin/python
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
'''calculate 4th order velocity moment (kurtosis) of 2D rings from a Walker dataset'''

from pylab import *
import numpy as np
from scipy.stats import kurtosis

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
        print >> vfil,'r','kappa_los(r)','error'

        # 30 iterations for drawing a given radius in bin
        mom4         = np.zeros((gpr.nbins,gpr.n))
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
                    mom4[i][k] = mom4[i-1][k]
                else:
                    mom4[i][k] = kurtosis(vlos1, axis=0, fisher=True, bias=False)

        for i in range(gpr.nbins):
            kappavel = np.sum(mom4[i])/gpr.n #[1]
            ab = np.sum(a[i])/(1.*gpr.n) #[1]
            if ab == 0:
                kappavelerror = p_edvlos[i-1] #[1]
                # attention! uses last error
            else:
                kappavelerror = kappavel/np.sqrt(ab) #[1]
            p_dvlos[i] = kappavel      #[1]
            p_edvlos[i]= kappavelerror #[1]


        for i in range(gpr.nbins):
            #             [rcore]  [1]               [1]
            print >> vfil,rbin[i], np.abs(p_dvlos[i]),np.abs(p_edvlos[i]) #/np.sqrt(n))
        vfil.close()

        if not gp.showplot_readout: continue

        ion(); subplot(111)
        print 'rbin = ',rbin,' rcore'
        print 'p_dvlos = ',p_dvlos
        print 'p_edvlos = ',p_edvlos
        plot(rbin,p_dvlos,'b',lw=1)
        fill_between(rbin,p_dvlos-p_edvlos,p_dvlos+p_edvlos,alpha=0.5,color='r') #[rcore],2*[1]

        xlabel(r'$r [\mathrm{rcore}]$')
        ylabel(r'$\langle\kappa_{\mathrm{LOS}}\rangle [1]$')
        ylim([-3,3])
        # xscale('log')
        xlim([0,3])
        #plt.legend(['\rho','\rho'],'lower left'); #title(dwarf)
        savefig(gpr.get_siglos_png(comp))
        if gpr.showplots:
            ioff();show();clf()


if __name__ == '__main__':
    gp.showplot_readout = True
    run()

