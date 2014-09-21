#!/usr/bin/env ipython3

##
# @file
# calculate surface mass density falloff of circular rings around center of mass
# calculate velocity dispersion of 2D rings from a Hernquist dataset
# calculate 4th order velocity moment (kurtosis) of 2D rings from a Hernquistdataset
# do this consistently with always the same sets of particles per bin

# (c) 2013 Pascal S.P. Steger

import sys, math, pdb
import numpy as np
from scipy.stats import kurtosis

import gr_params as gpr
import gl_file as gfile
from gl_helper import expDtofloat, bin_r_linear, bin_r_log, bin_r_const_tracers
from gl_class_files import *
from BiWeight import meanbiweight


def run(gp):
    # determine radius once and for all from all tracers
    R, Phi, vzall = np.loadtxt(gpr.fileposspherical[0],
                               comments='#',unpack=True) # 2*[Rscale], [km/s]
    Rscale0 = gfile.read_Rscale(gp.files.get_scale_file(0)) # [pc]

    # set number and size of (linearly spaced) bins
    Rmin = 0.                                         # [Rscale]
    Rmax = max(R) if gp.maxR<0 else 1.0*gp.maxR       # [Rscale]
    print('Rmax [Rscale] = ', Rmax)                   # [Rscale]
    R = R[(R<=Rmax)]                                  # [Rscale]

    # this must not be changed between readout and gravlite run
    # if you wish to change: set gp.getnewdata = True in gl_params.py
    if gpr.lograd:
        Binmin, Binmax, Rbin = bin_r_log(Rmax/gp.nipol, Rmax, gp.nipol)
        print(gp.nipol,' bins in log spacings')
    elif gp.consttr:
        Binmin, Binmax, Rbin = bin_r_const_tracers(R, len(R)/gp.nipol)
        print(len(R)/gp.nipol,' particles per bin')
    else:
        Binmin, Binmax, Rbin = bin_r_linear(Rmin, Rmax, gp.nipol)
        print(gp.nipol, ' bins in linear spacings')


    # volume of a circular ring from binmin to binmax
    Vol = np.zeros(gp.nipol)
    for k in range(gp.nipol):
        Vol[k] = np.pi*(Binmax[k]**2-Binmin[k]**2) # [Rscale^2]


    for pop in range(gpr.pops):
        print('#######  working on component ',pop)
        print('grh_MCMCbin: input: ',gpr.fileposspherical[pop])
        # start from data centered on COM already:
        if gfile.bufcount(gpr.fileposspherical[pop])<2: continue
        R, Phi, v = np.loadtxt(gpr.fileposspherical[pop],\
                               comments='#',unpack=True)
                               # [Rscale], [1], [km/s]
        
        # set maximum radius (if gp.maxR is set)
        Rmax = max(R) if gp.maxR<0 else 1.0*gp.maxR         # [Rscale]
        print('Rmax [Rscale] = ', Rmax)
        Rscalei = gfile.read_Rscale(gp.files.get_scale_file(pop)) # [pc]
        sel = (R<=Rmax)
        R = R[sel]; v = v[sel] # [Rscale], [km/s]
        totmass = 1.*len(R) # [Munit], Munit = 1/star
            
        Rs = R                   # + possible starting offset, [Rscale]
        vlos = v                 # + possible starting offset, [km/s]
        
        print('grh_MCMCbin: output density: ')
        tr = open(gp.files.get_ntracer_file(pop),'w')
        print(totmass, file=tr)
        tr.close()

        de, em, sigfil, kappafil = gfile.write_headers_2D(gp, pop)

        # gpr.n=30 iterations for getting random picked radius values
        Density = np.zeros((gp.nipol,gpr.n))
        dispvelocity = np.zeros((gp.nipol, gpr.n))
        mom4         = np.zeros((gp.nipol, gpr.n))
        tpb          = np.zeros((gp.nipol, gpr.n))
        zetaa        = np.zeros((gp.nipol, gpr.n))
        zetab        = np.zeros((gp.nipol, gpr.n))
        # 'a' shared by density, siglos, kappa calcs
        for k in range(gpr.n):
            Rsi = gpr.Rerr * np.random.randn(len(Rs)) + Rs # [Rscale]
            vlosi = gpr.vrerr * np.random.randn(len(vlos)) + vlos # [km/s]
            for i in range(gp.nipol):
                sel = np.argwhere(np.logical_and(Rsi * Rscalei >= Binmin[i] * Rscale0,\
                                                 Rsi * Rscalei <  Binmax[i] * Rscale0)).flatten() # [1]

                # gives positions of particles in bin i, in an array
                lsel = float(len(sel)) # is total number of particles in this bin
                tpb[i][k] = lsel
                Density[i][k] = (1.*len(sel))/Vol[i]*totmass # [Munit/Rscale0^2]
                vlos1 = vlosi[sel] # [km/s]

                if(len(sel)<=1):
                    dispvelocity[i][k] = dispvelocity[i-1][k]
                    mom4[i][k] = mom4[i-1][k]
                    zetaa[i][k] = zetaa[i-1][k]
                    zetab[i][k] = zetab[i-1][k]
                    # attention! should be 0, uses last value
                else:
                    dispvelocity[i][k] = meanbiweight(vlos1,ci_perc=68.4,ci_mean=True,ci_std=True)[1]
                                        # [km/s], see BiWeight.py
                    mom4[i][k] = kurtosis(vlos1, axis=0, fisher=False, bias=False) # [1]
                    # zetaa = TODO
                    # zetab = TODO

        ### TODO: use a method for this part
        # output density
        Dens0 = np.sum(Density[0])/(1.*gpr.n) # [Munit/Rscale^2]
        print('Dens0 = ', Dens0, '[Munit/Rscale^2]')
        crscale = open(gp.files.get_scale_file(pop),'r')
        Rscale = np.loadtxt(crscale, comments='#', unpack=False) # [pc]
        crscale.close()

        cdens = open(gp.files.get_scale_file(pop),'a')
        print(Dens0, file=cdens)               # [Munit/Rscale^2]
        Dens0pc = Dens0*Rscale**2
        print(Dens0pc, file=cdens)     # [Munit/pc^2]
        print(totmass, file=cdens)             # [Munit]
        cdens.close()

        tpbb0   = np.sum(tpb[0])/(1.*gpr.n)     # [1]
        Denserr0 = Dens0/np.sqrt(tpbb0)       # [Munit/Rscale^2]
        P_dens  = np.zeros(gp.nipol)
        P_edens = np.zeros(gp.nipol)
        for b in range(gp.nipol):
            Dens = np.sum(Density[b])/float(gpr.n) # [Munit/Rscale^2]
            tpbb   = np.sum(tpb[b])/float(gpr.n)       # [1]
            Denserr = Dens/np.sqrt(tpbb)       # [Munit/Rscale^2]
            if(math.isnan(Denserr)):
                P_dens[b] = P_dens[b-1]  # [1]
                P_edens[b]= P_edens[b-1] # [1]
            else:
                P_dens[b] = Dens/Dens0   # [1]
                P_edens[b]= Denserr/Dens0    # [1] #100/rbin would be artificial guess

            print(Rbin[b], Binmin[b], Binmax[b], P_dens[b], P_edens[b], file=de)
            # [Rscale], [dens0], [dens0]
            indr = (R<Binmax[b])
            Menclosed = 1.0*np.sum(indr)/totmass # for normalization to 1  # [totmass]
            Merr = Menclosed/np.sqrt(tpbb) # or artificial Menclosed/10 # [totmass]
            print(Rbin[b], Binmin[b], Binmax[b], Menclosed, Merr, file=em)
            # [Rscale], 2* [totmass]
        de.close()
        em.close()


        # output siglos
        p_dvlos = np.zeros(gp.nipol);        p_edvlos = np.zeros(gp.nipol)
        for b in range(gp.nipol):
            dispvel = np.sum(dispvelocity[b])/gpr.n # [km/s]
            ab = np.sum(tpb[b])/(1.*gpr.n)            # [1]
            if ab == 0:
                dispvelerr = p_edvlos[b-1] # [km/s]
                # attention! uses last error
            else:
                dispvelerr = dispvel/np.sqrt(ab) # [km/s]
            p_dvlos[b] = dispvel                 # [km/s]
            p_edvlos[b]= dispvelerr              # [km/s]

        maxsiglos = max(p_dvlos)          # [km/s]
        print('maxsiglos = ',maxsiglos,'[km/s]')
        fpars = open(gp.files.get_scale_file(pop),'a')
        print(maxsiglos, file=fpars)      # [km/s]
        fpars.close()
        
        for b in range(gp.nipol):
            #     [Rscale] [Rscale]  [Rscale]   [maxsiglos]                  [maxsiglos]
            print(Rbin[b], Binmin[b], Binmax[b], np.abs(p_dvlos[b]/maxsiglos), np.abs(p_edvlos[b]/maxsiglos), file=sigfil)
        sigfil.close()

        # output kurtosis kappa
        p_kappa = np.zeros(gp.nipol) # needed for plotting later
        p_ekappa = np.zeros(gp.nipol)
        for b in range(gp.nipol):
            kappavel = np.sum(mom4[b])/gpr.n # [1]
            ab = np.sum(tpb[b])/(1.*gpr.n)     # [1]
            if ab == 0:
                kappavelerr = p_edvlos[b-1] # [1]
                # attention! uses last error
            else:
                kappavelerr = np.abs(kappavel/np.sqrt(ab))
            p_kappa[b] = kappavel
            p_ekappa[b] = kappavelerr
            
            print(Rbin[b], Binmin[b], Binmax[b],\
                  kappavel, kappavelerr, file=kappafil) # 3*[Rscale], 2*[1]
        kappafil.close()

        if gpr.showplots:
            gpr.show_plots_dens_2D(pop, Rbin, P_dens, P_edens, Dens0pc)
            gpr.show_plots_sigma(pop, Rbin, p_dvlos, p_edvlos)
            gpr.show_plots_kappa(pop, Rbin, p_kappa, p_ekappa)

            
if __name__ == '__main__':
    gpr.showplots = True
    run(gp)

