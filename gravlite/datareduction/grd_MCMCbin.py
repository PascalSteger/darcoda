#!/usr/bin/env python3

##
# @file
# calculate surface mass density falloff of circular rings around center of mass
# calculate velocity dispersion of 2D rings from an observed dwarf dataset
# calculate 4th order velocity moment (kurtosis) of 2D rings
# do this consistently with always the same sets of particles per bin

# (c) 2013 Pascal S.P. Steger

import sys
import pdb
import math
import numpy as np
from scipy.stats import kurtosis
from pylab import *

import gr_params as gpr
import gl_file as gfile
from gl_helper import expDtofloat, bin_r_linear, bin_r_log, bin_r_const_tracers
from gl_class_files import *
from BiWeight import meanbiweight



def run(gp):
    xall,yall = np.loadtxt(gpr.get_com_file(0), skiprows=1, usecols=(0,1), unpack=True) # 2*[Rscale]
    # calculate 2D radius on the skyplane
    R = np.sqrt(xall**2+yall**2) # [Rscale]
    # set number and size of (linearly spaced) bins
    Rmin = 0. #[rscale]
    Rmax = max(R) if gp.maxR<0 else 1.0*gp.maxR # [Rscale]
    print('Rmax [Rscale] = ', Rmax)
    R = R[(R<Rmax)]

    # determine radius once and for all
    # this must not be changed between readout and gravlite run
    # if you wish to change: set gp.getnewdata = True in gl_params.py
    if gp.lograd:
        print(gp.nipol, ' bins in log spacings')
        Binmin, Binmax, Rbin = bin_r_log(Rmax/gp.nipol, Rmax, gp.nipol)
    elif gp.consttr:
        print(len(R)/gp.nipol,' particles per bin')
        Binmin, Binmax, Rbin = bin_r_const_tracers(R, len(R)/gp.nipol)
    else:
        print(gp.nipol, ' bins in linear spacings')
        Binmin, Binmax, Rbin = bin_r_linear(Rmin, Rmax, gp.nipol)


    # volume of a circular ring from binmin to binmax
    Vol = np.zeros(gp.nipol)
    for k in range(gp.nipol):
        Vol[k] = np.pi*(Binmax[k]**2-Binmin[k]**2) # [Rscale^2]


    for comp in range(gpr.ncomp):
        print('#######  working on component ',comp)
        print('input: ',gpr.get_com_file(comp))
        # start from data centered on COM already:
        if gfile.bufcount(gpr.get_com_file(comp))<2: continue
        x,y,v = np.loadtxt(gpr.get_com_file(comp),\
                           skiprows=1,usecols=(0,1,2),unpack=True)
                           # [rscale], [rscale], [km/s]

        # calculate 2D radius on the skyplane
        R = np.sqrt(x**2+y**2) #[rscale]
        
        # set maximum radius (if gp.maxR is set)
        Rmax = max(R) if gp.maxR<0 else 1.0*gp.maxR # [Rscale]
        print('Rmax [Rscale] = ', Rmax)
        sel = (R<=Rmax)
        x = x[sel]; y = y[sel]; v = v[sel]; R = R[sel] # [Rscale]
        totmass = 1.*len(x) # [munit], munit = 1/star
            
        Rs = R                   # + possible starting offset, [Rscale]
        vlos = v                 # + possible starting offset, [km/s]
        
        tr = open(gp.files.get_ntracer_file(comp),'w')
        print(totmass, file=tr)
        tr.close()

        de, em, sigfil, kappafil = gfile.write_headers(gp, comp)

        # gpr.n=30 iterations for getting random picked radius values
        Density = np.zeros((gp.nipol,gpr.n))
        dispvelocity = np.zeros((gp.nipol,gpr.n))
        mom4         = np.zeros((gp.nipol,gpr.n))
        a            = np.zeros((gp.nipol,gpr.n)) # shared by density, siglos, kappa calcs
        for k in range(gpr.n):
            Rsi = gpr.Rerror * np.random.randn(len(Rs)) + Rs # [Rscale]
            vlosi = gpr.vrerror * np.random.randn(len(vlos)) + vlos # [km/s]
            for i in range(gp.nipol):
                ind1 = np.argwhere(np.logical_and(Rsi >= Binmin[i],Rsi<Binmax[i])).flatten() # [1]
                Density[i][k] = (1.*len(ind1))/Vol[i]*totmass # [munit/rscale**2]
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
        Dens0 = np.sum(Density[0])/(1.*gpr.n) # [munit/Rscale^2]
        print('Dens0 = ', Dens0, '[munit/Rscale^2]')
        crscale = open(gp.files.get_scale_file(comp),'r')
        Rscale = np.loadtxt(crscale, comments='#', skiprows=1, unpack=False)
        crscale.close()

        cdens = open(gp.files.get_scale_file(comp),'a')
        print(Dens0, file=cdens)               # [munit/Rscale^2]
        Dens0pc = Dens0/Rscale**2              # [munis/pc^2]
        print(Dens0pc, file=cdens)             # [munit/pc^2]
        print(totmass, file=cdens)             # [munit]
        cdens.close()

        ab0   = np.sum(a[0])/(1.*gpr.n)     # [1]
        Denserr0 = Dens0/np.sqrt(ab0)       # [munit/Rscale^2]
        P_dens  = np.zeros(gp.nipol);  P_edens = np.zeros(gp.nipol)
        for b in range(gp.nipol):
            Dens = np.sum(Density[b])/(1.*gpr.n) # [munit/Rscale^2]
            ab   = np.sum(a[b])/(1.*gpr.n)       # [1]
            Denserr = Dens/np.sqrt(ab)       # [munit/Rscale^2]
            # compare data and analytic profile <=> get stellar density or mass ratio from Matt Walker
            Denserror = np.sqrt((Denserr/Dens0)**2+(Dens*Denserr0/(Dens0**2))**2) # [1]
            if(math.isnan(Denserror)):
                Denserror = 0. # [1]
                P_dens[b] = P_dens[b-1]  # [1]
                P_edens[b]= P_edens[b-1] # [1]
            else:
                P_dens[b] = Dens/Dens0   # [1]
                P_edens[b]= Denserror    # [1] #100/rbin would be artificial guess

            print(Rbin[b], Binmin[b], Binmax[b], P_dens[b], P_edens[b], file=de)
            # 3*[rscale], [dens0], [dens0]
            indr = (R<Binmax[b])
            Menclosed = 1.0*np.sum(indr)/totmass # for normalization to 1  #[totmass]
            Merror = Menclosed/np.sqrt(ab) # or artificial Menclosed/10 #[totmass]
            print(Rbin[b], Binmin[b], Binmax[b], Menclosed, Merror, file=em) # [Rscale], 2* [totmass]
        de.close()
        em.close()


        # output siglos
        p_dvlos = np.zeros(gp.nipol);        p_edvlos = np.zeros(gp.nipol)
        for b in range(gp.nipol):
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
        print('maxvlos = ', maxvlos, '[km/s]')
        fpars = open(gp.files.get_scale_file(comp),'a')
        print(maxvlos, file=fpars)          #[km/s]
        fpars.close()
        
        for b in range(gp.nipol):
            print(Rbin[b], Binmin[b], Binmax[b], np.abs(p_dvlos[b]/maxvlos),np.abs(p_edvlos[b]/maxvlos), file=sigfil)
            # 3*[rscale], 2*[maxvlos]
            # TODO: check uncommented /np.sqrt(n))
        sigfil.close()


        # output kurtosis kappa
        p_kappa = np.zeros(gp.nipol) # needed for plotting later
        p_ekappa = np.zeros(gp.nipol)
        for b in range(gp.nipol):
            kappavel = np.sum(mom4[b])/gpr.n #[1]
            ab = np.sum(a[b])/(1.*gpr.n) #[1]
            if ab == 0:
                kappavelerror = p_edvlos[b-1] #[1]
                # attention! uses last error
            else:
                kappavelerror = np.abs(kappavel/np.sqrt(ab)) #[1]
            p_kappa[b] = kappavel
            p_ekappa[b] = kappavelerror
            
            print(Rbin[b],Binmin[b],Binmax[b], kappavel, kappavelerror, file=kappafil) # [rscale], 2*[1]
            # TODO: /np.sqrt(n))
        kappafil.close()


    
        if not gpr.showplots: continue
        # plot density
        ion(); subplot(111)
        print('Rbin = ', Rbin)
        print('P_dens = ', P_dens)
        print('P_edens = ', P_edens)

        plot(Rbin, P_dens*Dens0pc, 'b', lw=1)
        lbound = (P_dens-P_edens)*Dens0pc; lbound[lbound<1e-6] = 1e-6
        ubound = (P_dens+P_edens)*Dens0pc
        fill_between(Rbin, lbound, ubound, alpha=0.5, color='r')
        yscale('log')
        # xlim([0, gp.maxR])
        # ylim([np.min(lbound),np.max(ubound)])
        xlabel(r'$R [R_c]$')
        ylabel(r'$\nu_{2D}(R) [\mathrm{Msun/pc/pc}]$')
        savefig(gpr.get_dens_png(i))
        ioff(); show(); clf()

        # plot siglos
        ion(); subplot(111)
        print('Rbin = ',Rbin,' Rscale')
        print('p_dvlos = ',p_dvlos,' km/s')
        print('p_edvlos = ',p_edvlos, 'km/s')
        plot(Rbin, p_dvlos, 'b', lw=1)
        fill_between(Rbin, p_dvlos-p_edvlos, p_dvlos+p_edvlos, alpha=0.5, color='r')
        # [rscale],2*[km/s]

        xlabel(r'$R [\mathrm{Rscale}]$')
        ylabel(r'$\langle\sigma_{\mathrm{LOS}}\rangle [\mathrm{km/s}]$')
        ylim([-1, 30])
        # xlim([0, 3])
        savefig(gpr.get_siglos_png(comp))
        ioff(); show(); clf()


        # plot kappa
        ion(); subplot(111)
        print('Rbin = ', Rbin, ' Rscale')
        print('p_kappa = ', p_kappa)
        print('p_ekappa = ', p_ekappa)
        plot(Rbin, p_kappa, 'b', lw=1)
        fill_between(Rbin, p_kappa-p_ekappa, p_kappa+p_ekappa, alpha=0.5, color='r')
        # [rscale], 2*[1]
        xlabel(r'$R [\mathrm{Rscale}]$')
        ylabel(r'$\langle\kappa_{\mathrm{LOS}}\rangle [1]$')
        ylim([0, 5.])
        # xlim([0, gp.maxR])
        savefig(gpr.get_kurtosis_png(comp))
        ioff(); show(); clf()


if __name__ == '__main__':
    gpr.showplots = True
    import gl_params
    gp = gl_params.Params()
    run(gp)

