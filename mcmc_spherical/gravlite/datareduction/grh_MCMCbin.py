#!/usr/bin/env python3
# (c) 2013 Pascal S.P. Steger
'''calculate surface mass density falloff of circular rings around center of mass'''
'''calculate velocity dispersion of 2D rings from a Hernquist dataset'''
'''calculate 4th order velocity moment (kurtosis) of 2D rings from a Hernquistdataset'''
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
    # determine radius once and for all from all tracers
    R, Phi, vzall = np.loadtxt(gpr.fileposspherical[0],
                               comments='#',unpack=True) # 2*[Rcore], [km/s]
    # set number and size of (linearly spaced) bins
    Rmin = 0. # [rcore]
    Rmax = max(R) if gpr.rprior<0 else 1.0*gpr.rprior # [Rcore]
    print('Rmax [Rcore] = ', Rmax)
    R = R[(R<=Rmax)]


    # this must not be changed between readout and gravlite run
    # if you wish to change: set gp.getnewdata = True in gl_params.py
    if gp.lograd:
        Binmin, Binmax, Rbin = bin_r_log(Rmax/gpr.nbins, Rmax, gpr.nbins)
        print(gpr.nbins,' bins in log spacings')
    elif gp.consttr:
        Binmin, Binmax, Rbin = bin_r_const_tracers(R, len(R)/gpr.nbins)
        print(len(R)/gpr.nbins,' particles per bin')
    else:
        Binmin, Binmax, Rbin = bin_r_linear(Rmin, Rmax, gpr.nbins)
        print(gpr.nbins, ' bins in linear spacings')

    # volume of a circular ring from binmin to binmax
    Vol = np.zeros(gpr.nbins)
    for k in range(gpr.nbins):
        Vol[k] = np.pi*(Binmax[k]**2-Binmin[k]**2) # [Rcore^2]


    for comp in range(gpr.ncomp):
        print('#######  working on component ',comp)
        print('grh_MCMCbin: input: ',gpr.fileposspherical[comp])
        # start from data centered on COM already:
        if gfile.bufcount(gpr.fileposspherical[comp])<2: continue
        R,Phi,v = np.loadtxt(gpr.fileposspherical[comp],\
                             comments='#',unpack=True)
                             # 2*[rcore], [km/s]
        
        # set maximum radius (if gpr.rprior is set)
        Rmax = max(R) if gpr.rprior<0 else 1.0*gpr.rprior # [Rcore]
        print('Rmax [Rcore] = ', Rmax)
        sel = (R<=Rmax)
        R = R[sel]; v = v[sel] # [Rcore], [km/s]
        totmass = 1.*len(R) # [munit], munit = 1/star
            
        Rs = R                   # + possible starting offset, [Rcore]
        vlos = v                 # + possible starting offset, [km/s]
        
        print('grh_MCMCbin: output density: ')
        print(gpr.get_ntracer_file(comp))
        tr = open(gpr.get_ntracer_file(comp),'w')
        print(totmass, file=tr)
        tr.close()

        print(gpr.filedenfalloff[comp])
        de = open(gpr.filedenfalloff[comp],'w')
        print('Rbin [Rcore]','Binmin [Rcore]','Binmax [Rcore]',\
              'Nu(R)/Nu(0) [1]','error [1]', file=de)

        print(gpr.filemass[comp])
        em = open(gpr.filemass[comp],'w')
        print('R [Rcore]','Binmin [Rcore]','Binmax [Rcore]',\
              'M(<Binmax) [Msun]','error [Msun]', file=em)


        print('grh_MCMCbin: output siglos: ',gpr.filesig[comp])
        sigfil = open(gpr.filesig[comp],'w')
        print('R [Rcore]','Binmin [Rcore]','Binmax [Rcore]',\
              'sigma_r(R) [km/s]','error [km/s]', file=sigfil)


        print('grh_MCMCbin: output kurtosis: ',gpr.filekappa[comp])
        kappafil = open(gpr.filekappa[comp],'w')
        print('R [Rcore]','Binmin [Rcore]','Binmax [Rcore]',\
              'kappa_los(R) [1]','error [1]', file=kappafil)


        # gpr.n=30 iterations for getting random picked radius values
        Density = np.zeros((gpr.nbins,gpr.n))
        dispvelocity = np.zeros((gpr.nbins,gpr.n))
        mom4         = np.zeros((gpr.nbins,gpr.n))
        a            = np.zeros((gpr.nbins,gpr.n))
        # 'a' shared by density, siglos, kappa calcs
        for k in range(gpr.n):
            Rsi = gpr.Rerror * np.random.randn(len(Rs)) + Rs # [Rcore]
            vlosi = gpr.vrerror * np.random.randn(len(vlos)) + vlos # [km/s]
            for i in range(gpr.nbins):
                ind1 = np.argwhere(np.logical_and(Rsi >= Binmin[i],Rsi<Binmax[i])).flatten() # [1]
                Density[i][k] = (1.*len(ind1))/Vol[i]*totmass # [munit/rcore**2]
                vlos1 = vlosi[ind1] # [km/s]

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
        Dens0 = np.sum(Density[0])/(1.*gpr.n) # [munit/Rcore^2]
        print('Dens0 = ',Dens0,'[munit/Rcore^2]')
        crcore = open(gpr.get_params_file(comp),'r')
        Rcore = np.loadtxt(crcore, comments='#', unpack=False)
        crcore.close()

        cdens = open(gpr.get_params_file(comp),'a')
        print(Dens0, file=cdens)               # [munit/Rcore^2]
        print(Dens0/Rcore**2, file=cdens)      # [munit/pc^2]
        print(totmass, file=cdens)             # [munit]
        cdens.close()

        ab0   = np.sum(a[0])/(1.*gpr.n)     # [1]
        Denserr0 = Dens0/np.sqrt(ab0)       # [munit/Rcore^2]
        P_dens  = np.zeros(gpr.nbins);  P_edens = np.zeros(gpr.nbins)
        for b in range(gpr.nbins):
            Dens = np.sum(Density[b])/(1.*gpr.n) # [munit/Rcore^2]
            ab   = np.sum(a[b])/(1.*gpr.n)       # [1]
            Denserr = Dens/np.sqrt(ab)       # [munit/Rcore^2]
            Denserror = np.sqrt((Denserr/Dens0)**2+(Dens*Denserr0/(Dens0**2))**2) # [1]
            if(math.isnan(Denserror)):
                Denserror = 0. # [1]
                P_dens[b] = P_dens[b-1]  # [1]
                P_edens[b]= P_edens[b-1] # [1]
            else:
                P_dens[b] = Dens/Dens0   # [1]
                P_edens[b]= Denserror    # [1] #100/rbin would be artificial guess

            print(Rbin[b],Binmin[b],Binmax[b],P_dens[b],P_edens[b], file=de)
            # [rcore], [dens0], [dens0]
            indr = (R<Binmax[b])
            Menclosed = 1.0*np.sum(indr)/totmass # for normalization to 1  #[totmass]
            Merror = Menclosed/np.sqrt(ab) # or artificial Menclosed/10 #[totmass]
            print(Rbin[b], Binmin[b], Binmax[b], Menclosed, Merror, file=em)
            # [Rcore], 2* [totmass]
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
        print('maxvlos = ',maxvlos,'[km/s]')
        fpars = open(gpr.get_params_file(comp),'a')
        print(maxvlos, file=fpars)          #[km/s]
        fpars.close()
        
        for b in range(gpr.nbins):
            #             [rcore]  [maxvlos]                  [maxvlos]
            print(Rbin[b],Binmin[b],Binmax[b], np.abs(p_dvlos[b]/maxvlos),np.abs(p_edvlos[b]/maxvlos), file=sigfil)
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
            
            print(Rbin[b],Binmin[b],Binmax[b],\
                  kappavel, kappavelerror, file=kappafil) # [rcore], 2*[1]
            # TODO: /np.sqrt(n))
        kappafil.close()


    


        if not gpr.showplots: continue
        # plot density
        ion(); subplot(111)
        print('Rbin = ',Rbin)
        print('P_dens = ',P_dens)
        print('P_edens = ',P_edens)

        plot(Rbin,P_dens,'b',lw=1)
        lbound = P_dens-P_edens; lbound[lbound<1e-6] = 1e-6
        ubound = P_dens+P_edens; 
        fill_between(Rbin,lbound,ubound,alpha=0.5,color='r')
        yscale('log')
        xlim([0,3.])
        ylim([np.min(lbound),np.max(ubound)])
        xlabel(r'$R [R_c]$')
        ylabel(r'$\nu_{2D}(R)/\nu_{2D}(0)$')
        savefig(gpr.get_dens_png(comp))
        from gl_analytic import Sigma_anf
        plot(Rbin, Sigma_anf(Rbin*Rcore))  # argument in [pc] !

        if gpr.showplots:
            ioff(); show(); clf()

        # plot siglos
        ion(); subplot(111)
        print('Rbin = ',Rbin,' Rcore')
        print('p_dvlos = ',p_dvlos,' km/s')
        print('p_edvlos = ',p_edvlos, 'km/s')
        plot(Rbin,p_dvlos,'b',lw=1)
        fill_between(Rbin,p_dvlos-p_edvlos,p_dvlos+p_edvlos,alpha=0.5,color='r') #[rcore],2*[km/s]

        xlabel(r'$R [\mathrm{Rcore}]$')
        ylabel(r'$\langle\sigma_{\mathrm{LOS}}\rangle [\mathrm{km/s}]$')
        ylim([0.,1.5*max(p_dvlos)])
        xlim([0,3])
        savefig(gpr.get_siglos_png(comp))
        from gl_analytic import sig_los_anf
        plot(Rbin,sig_los_anf(Rbin*Rcore)) # argument must be [pc] !

        if gpr.showplots:
            ioff();show();clf()


        # plot kappa
        ion(); subplot(111)
        print('Rbin = ',Rbin,' Rcore')
        print('p_kappa = ',p_kappa)
        print('p_ekappa = ',p_ekappa)
        plot(Rbin,p_kappa,'b',lw=1)
        fill_between(Rbin,p_kappa-p_ekappa,p_kappa+p_ekappa,alpha=0.5,color='r') #[rcore],2*[1]
        xlabel(r'$R [\mathrm{Rcore}]$')
        ylabel(r'$\langle\kappa_{\mathrm{LOS}}\rangle [1]$')
        ylim([0,5])
        xlim([0,3])
        savefig(gpr.get_kurtosis_png(comp))
        if gpr.showplots:
            ioff();show();clf()


if __name__ == '__main__':
    gpr.showplots = True
    run()

