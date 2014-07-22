#!/usr/bin/env ipython3

##
# @file
# calculate density falloff for Hernquist datasets

# (c) 2013 ETHZ, Pascal S.P. Steger, psteger@phys.ethz.ch

import numpy as np
import pdb
import multiprocessing as mp
from pylab import *


import gl_params
gp = gl_params.Params()
import gr_params as gpr
from gl_helper import expDtofloat, bin_r_linear, bin_r_log, bin_r_const_tracers

def run():
    for comp in range(gpr.ncomp):
        print('input:',gpr.fileposspherical[comp])
        R, Phi, vlos = np.loadtxt(gpr.fileposspherical[comp],\
                                  comments='#', unpack=True)
        Totmass = 1.*len(R) # [Munit], Munit = 1/star
        # Rs=gpr.Rerr*np.random.randn(len(r))+r
        # Rs not changed later on, ever: misplacement, for all realizations. wanted?
        # yes, we scatter radii in foo_pool
        Rs = R;    Rmin = min(Rs);     Rmax = max(Rs)
        if gpr.lograd:
            print(gp.nipol,' bins in log spacings')
            Binmin, Binmax, Rbin = bin_r_log(Rmax/gp.nipol, Rmax, gp.nipol)
        elif gp.consttr:
            Binmin, Binmax, Rbin = bin_r_const_tracers(Rs, len(Rs)/gp.nipol)
            print(len(R)/gp.nipol,' particles per bin')
        else:
            print(gp.nipol, ' bins in linear spacings')
            Binmin, Binmax, Rbin = bin_r_linear(Rmin, Rmax, gp.nipol)

        # volume of a bin with height binlength, 2D
        Vol = np.zeros(gpr.bins)
        for i in range(gpr.bins):
            Vol[i] = np.pi*(Binmax[i]**2-Binmin[i]**2)

        # gpr.n=30 iterations for getting random picked radius values
        Density = np.zeros((gp.nipol,gpr.n))
        A       = np.zeros((gp.nipol,gpr.n)) # shared by density, siglos, kappa calcs
        for k in range(gpr.n):
            Rsi = gpr.Rerr * np.random.randn(len(Rs)) + Rs # [Rscale]
            vlosi = gpr.vrerr * np.random.randn(len(vlos)) + vlos
            for i in range(gp.nipol):
                ind1 = np.argwhere(np.logical_and(Rsi>=Binmin[i], Rsi<Binmax[i])).flatten() # [1]
                Density[i][k] = (1.*len(ind1))/Vol[i]*Totmass # [Munit/Rscale^2]
                vlos1 = vlosi[ind1]                           # [km/s]
                A[i][k] = 1.*len(ind1)                        # [1]

        # output density
        Dens0 = np.sum(Density[0])/(1.*gpr.n) # [Munit/Rscale^3]
        print('Dens0 = ',Dens0,' [Munit/Rscale^2]')
        crscale = open(gp.files.get_scale_file(comp),'r')
        Rscale = np.loadtxt(cscale, comments='#', unpack=False)
        crscale.close()

        cdens = open(gp.files.get_scale_file(comp),'a')
        print(Dens0, file=cdens)               # [Munit/Rscale^2]
        print(Dens0/Rscale**2, file=cdens)      # [Munit/pc^2]
        print(Totmass, file=cdens)             # [Munit]
        cdens.close()

        de = open(gp.files.Sigfiles[comp], 'w')
        print('# Rbin [Rscale]','Binmin [Rscale]','Binmax [Rscale]',
              'Sig(R)/Sig(0) [1]','error [1]', file=de)
    
        em = open(gp.files.massfiles[comp], 'w')
        print('# R [Rscale]','Binmin [Rscale]','Binmax [Rscale]',\
              'M(<Binmax) [Munit]','error [Munit]', file=em)
    
        AB0   = np.sum(A[0])/(1.*gpr.n)     # [1]
        Denserr0 = Dens0/np.sqrt(AB0)       # [Munit/Rscale^3]
        P_dens  = np.zeros(gp.nipol);  P_edens = np.zeros(gp.nipol)
        for b in range(gp.nipol):
            Dens = np.sum(Density[b])/(1.*gpr.n) # [Munit/Rscale^3]
            AB   = np.sum(A[b])/(1.*gpr.n)       # [1]
            Denserr = Dens/np.sqrt(AB)       # [Munit/Rscale^3]

            if(np.isnan(Denserr):
                P_dens[b] = P_dens[b-1]  # [1]
                P_edens[b]= P_edens[b-1] # [1]
            else:
                P_dens[b] = Dens/Dens0   # [1]
                P_edens[b]= Denserr/Dens0    # [1] #100/rbin would be artificial guess

            print(Rbin[b], Binmin[b], Binmax[b], P_dens[b], P_edens[b], file=de)
            # [Rscale], 2*[dens0]

            # normalization to 1:
            indr = (R<Binmax[b])
            menclosed = 1.0*np.sum(indr)/Totmass  # [Totmass]
            
            merr = menclosed/np.sqrt(AB) # artificial menclosed/10 # [Totmass]
            print(Rbin[b], Binmin[b], Binmax[b], menclosed, merr, file=em)
            # [Rscale], 2*[Totmass]
        de.close()
        em.close()
        
        if gpr.showplots:
            gpr.show_plots_dens_2D(Rbin, P_dens, P_edens, gp)
    return


if __name__=="__main__":
    gpr.showplots = True
    run()
