#!/usr/bin/env ipython3

##
# @file
# calculate surface mass density falloff of circular rings
# around center of mass calculate velocity dispersion of 2D rings from
# a Walker dataset calculate 4th order velocity moment (kurtosis) of
# 2D rings from a Walker dataset do this consistently with always the
# same sets of particles per bin

# (c) 2012-2014 ETHZ, Pascal S.P. Steger

import sys, pdb
import numpy as np
from scipy.stats import kurtosis
from pylab import *

import gr_params as gpr
import gl_file as gfile
from gl_helper import expDtofloat, add_errors
from BiWeight import meanbiweight

def set_bndry(R, gp):
    Rmin = 0. #[Rscale]
    if gp.maxR < 0:
        Rmax = max(R)
    else:
        Rmax = 1.0*gp.maxR # [Rscale]
    print('Rmax [Rscale] = ', Rmax)
    return Rmin, Rmax
## \fn set_bndry(R, gp)
# set number and size of (linearly spaced) bins
# @param R
# @param gp global parameters


def run(gp):
    xall,yall = np.loadtxt(gpr.get_com_file(0), skiprows=1, usecols=(0,1), unpack=True)
    # 2*[Rscale_0]

    Rscale0 = gfile.read_Rscale(gp.files.get_scale_file(0)) # [pc]
    
    # calculate 2D radius on the skyplane
    R = np.sqrt(xall**2+yall**2) # [Rscale0]
    Rmin, Rmax = set_bndry(R, gp) # [Rscale0]
    R = R[(R<Rmax)] # [Rscale_0]
    Binmin, Binmax, Rbin = gpr.determine_radius(R, Rmin, Rmax, gp) # [Rscale0]
    Vol = gpr.volume_circular_ring(Binmin, Binmax, gp) # [Rscale0^2]

    for comp in range(gpr.ncomp):
        if gfile.empty(gpr.get_com_file(comp)): continue
        print('####### working on component ',comp)

        # start from data centered on COM already:
        print('input: ', gpr.get_com_file(comp))
        x,y,v = np.loadtxt(gpr.get_com_file(comp),\
                           skiprows=1,usecols=(0,1,2),unpack=True) #[Rscale_i], [Rscale_i], [km/s]

        R = np.sqrt(x**2+y**2) #[Rscale_i]
        
        # set maximum radius (if gp.maxR is set)
        Rmin, Rmax = set_bndry(R, gp) # [Rscale_i]
        Rscalei = gfile.read_Rscale(gp.files.get_scale_file(comp)) # [pc]
        sel = (R * Rscalei <= Rmax * Rscale0) # [pc]
        x = x[sel]; y = y[sel]; v = v[sel]; R = R[sel] # [Rscale_i], [km/s]
        totmass = float(len(x)) # [Munit] We assume Mstar = 1 Munit for all stars
        Rs = R                  #+possible starting offset, [Rscale_i]
        vlos = v                #+possible starting offset, [km/s]

        gfile.write_tracer_file(gp.files.get_ntracer_file(comp), totmass)
        de, em, sigfil, kappafil = gfile.write_headers_2D(gp, comp)

        # tracers per bin, shared by density, siglos, kappa
        # calculations
        tpb     = np.zeros((gp.nipol, gpr.n))
        Density = np.zeros((gp.nipol, gpr.n))
        sigma   = np.zeros((gp.nipol, gpr.n))
        kappa   = np.zeros((gp.nipol, gpr.n))
        zetaa   = np.zeros((gp.nipol, gpr.n))
        zetab   = np.zeros((gp.nipol, gpr.n))
        # gpr.n=30 iterations for getting randomly picked radius and
        # v_LOS values
        # this gives us a handle on the stochastic error, which has a gaussian / normal distribution
        for k in range(gpr.n):
            Rsi   = add_errors(Rs,   gpr.Rerr) # [Rscale_i]
            vlosi = add_errors(vlos, gpr.vrerr) # [km/s]
            for i in range(gp.nipol):
                sel = np.argwhere(np.logical_and(Rsi*Rscalei >= Binmin[i]*Rscale0,\
                                                 Rsi*Rscalei <  Binmax[i]*Rscale0)).flatten() # [pc]
                # gives positions of particles in bin i, in an array
                lsel = float(len(sel)) # is total number of particles in this bin
                tpb[i][k] = lsel
                Density[i][k] = lsel/Vol[i] # [Munit/Rscale_0^2]

                if(lsel<=1): # missing values => use last value
                    sigma[i][k] = sigma[i-1][k]
                    kappa[i][k] = kappa[i-1][k]
                    zetaa[i][k] = zetaa[i-1][k]
                    zetab[i][k] = zetab[i-1][k]
                else:
                    sigma[i][k] = meanbiweight(vlosi[sel], ci_perc=68.4, ci_mean=True, ci_std=True)[1]
                                        # [km/s], see BiWeight.py
                    kappa[i][k] = kurtosis(vlosi[sel], axis=0, fisher=False, bias=False) # [1]
                    # zetaa[i][k] = TODO
                    # zetab[i][k] = TODO

        # central density averaged over gpr.n iterations
        Dens0 = np.sum(Density[0])/float(gpr.n) # [Munit/Rscale0^2]
        Dens0pc = Dens0/Rscale0**2              # [Munit/pc^2]
        gfile.write_density(gp.files.get_scale_file(comp), Dens0pc, totmass)

        # number of tracers in central bin, it. averaged
        tpb0   = np.sum(tpb[0])/float(gpr.n)     # [1]

        # error on the central density (error of the mean)
        Denserr0 = Dens0/np.sqrt(tpb0)       # [Munit/Rscale^2]
        P_dens  = np.zeros(gp.nipol);  P_edens = np.zeros(gp.nipol)
        for b in range(gp.nipol):
            Dens = np.sum(Density[b])/float(gpr.n) # [Munit/Rscale^2]
            tpbb   = np.sum(tpb[b])/float(gpr.n)       # [1]
            Denserr = Dens/np.sqrt(tpbb)       # [Munit/Rscale^2]
            
            if(np.isnan(Denserr)):
                P_dens[b] = P_dens[b-1]  # [1]
                P_edens[b]= P_edens[b-1] # [1]
            else:
                P_dens[b] = Dens/Dens0   # [1]
                P_edens[b]= Denserr/Dens0    # [1]

            print(Rbin[b], Binmin[b], Binmax[b], P_dens[b], P_edens[b], file=de)
            # 3*[Rscale0], [dens0], [dens0]
            indr = (R<Binmax[b])
            Menclosed = 1.0*np.sum(indr)/totmass # for normalization
                                                 # to 1 #[totmass]
            Merr = Menclosed/np.sqrt(tpbb) # or artificial
                                           # Menclosed/10 #[totmass]
            print(Rbin[b], Binmin[b], Binmax[b], Menclosed, Merr, file=em)
            # [Rscale], 2* [totmass]
        de.close()
        em.close()

        # output siglos
        p_dvlos = np.zeros(gp.nipol);        p_edvlos = np.zeros(gp.nipol)
        for b in range(gp.nipol):
            dispvel = np.sum(sigma[b])/gpr.n #[km/s] # mean sigma
            tpbb = np.sum(tpb[b])/(float(gpr.n)) #[1]
            if tpbb == 0:
                dispvelerr = p_edvlos[b-1] #[km/s] attention! uses
                # last error
            else:
                dispvelerr = dispvel/np.sqrt(tpbb) #[km/s]
            p_dvlos[b] = dispvel      #[km/s]
            p_edvlos[b]= dispvelerr #[km/s]

        maxsiglos = max(p_dvlos) #[km/s]
        print('maxsiglos = ', maxsiglos, '[km/s]')
        fpars = open(gp.files.get_scale_file(comp),'a')
        print(maxsiglos, file=fpars)          #[km/s]
        fpars.close()
        
        for b in range(gp.nipol):
            print(Rbin[b], Binmin[b], Binmax[b], \
                  np.abs(p_dvlos[b]/maxsiglos),np.abs(p_edvlos[b]/maxsiglos), \
                  file=sigfil)
            # 3*[rscale], 2*[maxsiglos]
        sigfil.close()


        # output kurtosis kappa
        p_kappa = np.zeros(gp.nipol) # needed for plotting later
        p_ekappa = np.zeros(gp.nipol)
        for b in range(gp.nipol):
            kappavel = np.sum(kappa[b])/gpr.n #[1]
            tpbb = np.sum(tpb[b])/float(gpr.n) #[1]
            if tpbb == 0:
                kappavelerr = p_edvlos[b-1] #[1] attention! uses
                # last error
            else:
                kappavelerr = np.abs(kappavel/np.sqrt(tpbb)) #[1]
            p_kappa[b] = kappavel
            p_ekappa[b] = kappavelerr
            
            print(Rbin[b], Binmin[b], Binmax[b],\
                  kappavel, kappavelerr, \
                  file=kappafil)
            # [rscale], 2*[1]
        kappafil.close()
    
        if gpr.showplots:
            gpr.show_plots_dens_2D(comp, Rbin, P_dens, P_edens, Dens0pc)
            gpr.show_plots_sigma(comp, Rbin, p_dvlos, p_edvlos)
            gpr.show_plots_kappa(comp, Rbin, p_kappa, p_ekappa)
            
## \fn run(gp)
# main functionality
# @param gp global parameters
            

if __name__ == '__main__':
    gpr.showplots = True
    import gl_params
    gp = gl_params.Params()
    run(gp)

