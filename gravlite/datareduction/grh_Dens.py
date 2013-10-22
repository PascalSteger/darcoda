#!/usr/bin/env python3

##
# @file
# calculate density falloff for Hernquist datasets

# (c) 2013 ETHZ, Pascal S.P. Steger, psteger@phys.ethz.ch

import numpy as np
import math
import pdb
import multiprocessing as mp
from pylab import *


import gl_params as gp
import gr_params as gpr
from gl_helper import expDtofloat, bin_r_linear, bin_r_log, bin_r_const_tracers

def run():
    for comp in range(gpr.ncomp):
        print('input:',gpr.fileposspherical[comp])
        R, Phi, vlos = np.loadtxt(gpr.fileposspherical[comp],\
                                  comments='#', unpack=True)
        Totmass = 1.*len(R) # [munit], munit = 1/star
        # Rs=gpr.Rerror*np.random.randn(len(r))+r
        # Rs not changed later on, ever: misplacement, for all realizations. wanted?
        # yes, we scatter radii in foo_pool
        Rs = R;    Rmin = min(Rs);     Rmax = max(Rs)
        if gp.lograd:
            print(gpr.nbins,' bins in log spacings')
            Binmin, Binmax, Rbin = bin_r_log(Rmax/gpr.nbins, Rmax, gpr.nbins)
        elif gp.consttr:
            Binmin, Binmax, Rbin = bin_r_const_tracers(Rs, len(Rs)/gpr.nbins)
            print(len(R)/gpr.nbins,' particles per bin')
        else:
            print(gpr.nbins, ' bins in linear spacings')
            Binmin, Binmax, Rbin = bin_r_linear(Rmin, Rmax, gpr.nbins)

        # volume of a bin with height binlength, 2D
        Vol = np.zeros(gpr.bins)
        for i in range(gpr.bins):
            Vol[i] = np.pi*(Binmax[i]**2-Binmin[i]**2)

        # gpr.n=30 iterations for getting random picked radius values
        Density = np.zeros((gpr.nbins,gpr.n))
        A       = np.zeros((gpr.nbins,gpr.n)) # shared by density, siglos, kappa calcs
        for k in range(gpr.n):
            Rsi = gpr.Rerror * np.random.randn(len(Rs)) + Rs # [Rcore]
            vlosi = gpr.vrerror * np.random.randn(len(vlos)) + vlos
            for i in range(gpr.nbins):
                ind1 = np.argwhere(np.logical_and(Rsi>=Binmin[i], Rsi<Binmax[i])).flatten() # [1]
                Density[i][k] = (1.*len(ind1))/Vol[i]*Totmass # [munit/Rcore^2]
                vlos1 = vlosi[ind1]                           # [km/s]
                A[i][k] = 1.*len(ind1)                        # [1]

        # output density
        Dens0 = np.sum(Density[0])/(1.*gpr.n) # [munit/Rcore^3]
        print('Dens0 = ',Dens0,' [munit/Rcore^2]')
        crcore = open(gpr.get_params_file(comp),'r')
        Rcore = np.loadtxt(crcore, comments='#', unpack=False)
        crcore.close()

        cdens = open(gpr.get_params_file(comp),'a')
        print(Dens0,file=cdens)               # [munit/Rcore^2]
        print(Dens0/Rcore**2,file=cdens)      # [munit/pc^2]
        print(Totmass,file=cdens)             # [munit]
        cdens.close()

        print(gpr.get_dens_file(comp))
        de = open(gpr.get_dens_file(comp),'w')
        print('# Rbin [Rcore]','Binmin [Rcore]','Binmax [Rcore]',
              'Nu(R)/Nu(0) [1]','error [1]', file=de)
    
        print(gpr.get_enc_mass_file(comp))
        em = open(gpr.get_enc_mass_file(comp),'w')
        print('# R [Rcore]','Binmin [Rcore]','Binmax [Rcore]',\
              'M(<Binmax) [Msun]','error [Msun]', file=em)

    
        AB0   = np.sum(A[0])/(1.*gpr.n)     # [1]
        Denserr0 = Dens0/np.sqrt(AB0)       # [munit/Rcore^3]
        P_dens  = np.zeros(gpr.nbins);  P_edens = np.zeros(gpr.nbins)
        for b in range(gpr.nbins):
            Dens = np.sum(Density[b])/(1.*gpr.n) # [munit/Rcore^3]
            AB   = np.sum(A[b])/(1.*gpr.n)       # [1]
            Denserr = Dens/np.sqrt(AB)       # [munit/Rcore^3]
            Denserror = np.sqrt((Denserr/Dens0)**2+(Dens*Denserr0/(Dens0**2))**2) #[1]
            if(math.isnan(Denserror)):
                Denserror = 0. # [1]
                P_dens[b] = P_dens[b-1]  # [1]
                P_edens[b]= P_edens[b-1] # [1]
            else:
                P_dens[b] = Dens/Dens0   # [1]
                P_edens[b]= Denserror    # [1] #100/rbin would be artificial guess

            print(Rbin[b], Binmin[b], Binmax[b], P_dens[b], P_edens[b], file=de)
            # [Rcore], 2*[dens0]

            # normalization to 1:
            indr = (R<Binmax[b])
            menclosed = 1.0*np.sum(indr)/Totmass  # [Totmass]
            
            merror = menclosed/np.sqrt(AB) # artificial menclosed/10 # [Totmass]
            print(Rbin[b], Binmin[b], Binmax[b], menclosed, merror, file=em)
            # [rcore], 2*[Totmass]
            # TODO: check: take rbinmax for MCMC?
        de.close()
        em.close()
        
        if not gpr.showplots:
            continue
        
        # plot density
        ion(); subplot(111)
        print('Rbin = ',Rbin)
        print('P_dens = ',P_dens)
        print('P_edens = ',P_edens)
        
        plot(Rbin,P_dens,'b',lw=1)
        lbound = P_dens-P_edens; lbound[lbound<1e-6] = 1e-6
        ubound = P_dens+P_edens; 
        fill_between(Rbin, lbound, ubound, alpha=0.5, color='r')
        yscale('log')
        # xlim([0,3.])
        # ylim([np.min(lbound),np.max(ubound)])
        xlabel(r'$R [R_c]$')
        ylabel(r'$\nu_{2D}(R)/\nu_{2D}(0)$')
        savefig(gpr.get_dens_png(i))
        if gpr.showplots:
            ioff(); show(); clf()
    return
            
if __name__=="__main__":
    run()
