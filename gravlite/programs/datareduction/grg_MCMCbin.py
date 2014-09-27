#!/usr/bin/env ipython3

##
# @file
# calculate surface mass density falloff of circular rings around center of mass
# calculate velocity dispersion of 2D rings from a Walker dataset
# calculate 4th order velocity moment (kurtosis) of 2D rings from a Walker dataset
# do this consistently with always the same sets of particles per bin

# (c) 2013 Pascal S.P. Steger

import sys, pdb
import numpy as np
from scipy.stats import kurtosis
from pylab import *

import gr_params as gpr
import gl_file as gfile
from gl_helper import expDtofloat, add_errors
from gl_class_files import *
from BiWeight import meanbiweight
from gl_project import Rho_INT_rho



def run(gp):
    xall, yall = np.loadtxt(gpr.get_com_file(0), skiprows=1, usecols=(0,1), unpack=True)
    # 2*[Rscale]
    # calculate 2D radius on the skyplane
    R = np.sqrt(xall**2+yall**2) # [Rscale0]
    # set number and size of (linearly spaced) bins
    Rmin = 0. # [Rscale0]
    Rmax = max(R) if gp.maxR<0 else 1.0*gp.maxR # [Rscale0]
    R = R[(R<Rmax)]             # [Rscale0] exclude any NaNs

    Binmin, Binmax, Rbin = gpr.determine_radius(R, Rmin, Rmax, gp) # [Rscale0]
    Vol = gpr.volume_circular_ring(Binmin, Binmax, gp) # [Rscale0^2]
    Rscale0 = gfile.read_Xscale(gp.files.get_scale_file(0)) # [pc]

    for pop in range(gp.pops+1):
        print('#######  working on component ',pop)
        print('input: ',gpr.get_com_file(pop))
        # start from data centered on COM already:
        if gfile.bufcount(gpr.get_com_file(pop))<2: continue
        x,y,v = np.loadtxt(gpr.get_com_file(pop),\
                           skiprows=1,usecols=(0,1,2),unpack=True)
                           # [Rscalei], [Rscalei], [km/s]

        # calculate 2D radius on the skyplane
        R = np.sqrt(x**2+y**2) # [Rscalei]
        Rscalei = gfile.read_Xscale(gp.files.get_scale_file(pop)) # [pc]

        # set maximum radius (if gp.maxR is set)
        Rmax = max(R) if gp.maxR < 0 else 1.0*gp.maxR # [Rscale0]
        sel = (R * Rscalei <= Rmax * Rscale0) # [pc]
        x = x[sel]; y = y[sel]; v = v[sel]; R = R[sel] # [Rscalei]
        totmass = float(len(x)) # [Munit], Munit = 1/star

        Rs = R                   # + possible starting offset, [Rscalei]
        vlos = v                 # + possible starting offset, [km/s]

        gfile.write_tracer_file(gp.files.get_ntracer_file(pop), totmass)
        f_Sig, f_nu, f_mass, f_sig, f_kap = gfile.write_headers_2D(gp, pop)

        tpb     = np.zeros((gp.nipol, gpr.n)) # tracers per bin
        Density = np.zeros((gp.nipol, gpr.n))
        sigma   = np.zeros((gp.nipol, gpr.n))
        kappa   = np.zeros((gp.nipol, gpr.n))
        zetaa   = np.zeros((gp.nipol, gpr.n))
        zetab   = np.zeros((gp.nipol, gpr.n))

        # gpr.n=30 iterations for getting random picked radius values
        for k in range(gpr.n):
            Rsi   = add_errors(Rs,   gpr.Rerr)  # [Rscalei]
            vlosi = add_errors(vlos, gpr.vrerr) # [km/s]
            for i in range(gp.nipol):
                sel = np.argwhere(np.logical_and(Rsi*Rscalei >= Binmin[i]*Rscale0,\
                                                 Rsi*Rscalei <  Binmax[i]*Rscale0)).flatten() # [1]
                tpb[i][k] = float(len(sel)) #[1]
                Density[i][k] = float(len(sel))*totmass/Vol[i] # [Munit/Rscale0^2]

                if(len(sel)<=1):
                    sigma[i][k] = sigma[i-1][k]
                    kappa[i][k] = kappa[i-1][k]
                    zetaa[i][k] = zetaa[i-1][k]
                    zetab[i][k] = zetab[i-1][k]
                    # attention! should be 0, uses last value
                else:
                    sigma[i][k] = meanbiweight(vlosi[sel], ci_perc=68.4,\
                                               ci_mean=True, ci_std=True)[1]
                                        # [km/s], see BiWeight.py
                    kappa[i][k] = kurtosis(vlosi[sel], axis=0, fisher=False, bias=False) # [1]
                    # zetaa[i][k] = TODO
                    # zetab[i][k] = TODO


        # output density
        Dens0 = np.sum(Density[0])/(1.*gpr.n) # [Munit/Rscale^2]
        Dens0pc = Dens0/Rscale0**2              # [munis/pc^2]
        gfile.write_Sig_scale(gp.files.get_scale_file(pop), Dens0pc, totmass)

        tpb0   = np.sum(tpb[0])/float(gpr.n)     # [1]
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
                P_edens[b]= Denserr/Dens0    # [1] #100/rbin would be artificial guess

            print(Rbin[b], Binmin[b], Binmax[b], \
                  P_dens[b], P_edens[b], \
                  file=f_Sig)
            indr = (R<Binmax[b])
            Menclosed = float(np.sum(indr))/totmass # for normalization to 1#[totmass]
            Merr = Menclosed/np.sqrt(tpbb) # or artificial Menclosed/10 #[totmass]
            print(Rbin[b], Binmin[b], Binmax[b], \
                  Menclosed, Merr,\
                  file=f_mass)
        f_Sig.close()
        f_mass.close()


        # deproject Sig to get nu
        numedi = Rho_INT_rho(Rbin*Rscalei, Dens0pc*P_dens, gp)
        numin  = Rho_INT_rho(Rbin*Rscalei, Dens0pc*(P_dens-P_edens), gp)
        numax  = Rho_INT_rho(Rbin*Rscalei, Dens0pc*(P_dens+P_edens), gp)

        nu0pc  = numedi[0]
        gfile.write_nu_scale(gp.files.get_scale_file(pop), nu0pc)

        nuerr  = numax-numedi
        for b in range(gp.nipol):
            print(Rbin[b], Binmin[b], Binmax[b],\
                  numedi[b]/nu0pc, nuerr[b]/nu0pc, \
                  file = f_nu)
        f_nu.close()

        # output siglos
        p_dvlos = np.zeros(gp.nipol)
        p_edvlos = np.zeros(gp.nipol)
        for b in range(gp.nipol):
            dispvel = np.sum(sigma[b])/gpr.n #[km/s]
            tpbb = np.sum(tpb[b])/(1.*gpr.n) #[1]
            if tpbb == 0:
                dispvelerr = p_edvlos[b-1] #[km/s]
                # attention! uses last error
            else:
                dispvelerr = dispvel/np.sqrt(tpbb) #[km/s]
            p_dvlos[b] = dispvel      #[km/s]
            p_edvlos[b]= dispvelerr #[km/s]

        maxsiglos = max(p_dvlos) #[km/s]
        print('maxsiglos = ',maxsiglos,'[km/s]')
        fpars = open(gp.files.get_scale_file(pop),'a')
        print(maxsiglos, file=fpars)          #[km/s]
        fpars.close()

        for b in range(gp.nipol):
            print(Rbin[b], Binmin[b], Binmax[b], \
                  p_dvlos[b]/maxsiglos, p_edvlos[b]/maxsiglos, \
                  file=f_sig)
        f_sig.close()


        # output kurtosis kappa
        p_kappa = np.zeros(gp.nipol) # needed for plotting later
        p_ekappa = np.zeros(gp.nipol)
        for b in range(gp.nipol):
            kappavel = np.sum(kappa[b])/gpr.n #[1]
            ab = np.sum(tpb[b])/(1.*gpr.n) #[1]
            if ab == 0:
                kappavelerr = p_edvlos[b-1] #[1]  # TODO: /np.sqrt(n))
                # attention! uses last error
            else:
                kappavelerr = np.abs(kappavel/np.sqrt(ab)) #[1]
            p_kappa[b] = kappavel
            p_ekappa[b] = kappavelerr
            print(Rbin[b], Binmin[b], Binmax[b],\
                  kappavel, kappavelerr, file=f_kap) # [rscale], 2*[1]
        f_kap.close()


        if gpr.showplots:
            import gl_analytic as ga

            Ri = Rbin * Rscalei
            Sig_dat = Dens0pc*P_dens
            Sig_err = Dens0pc*P_edens

            #fill_between(Ri, Sig_dat-Sig_err, Sig_dat+Sig_err,\
            #             color='red', alpha=0.5)
            loglog(Ri, Sig_dat, 'r.-')

            nu_dat = Rho_INT_rho(Ri, Dens0pc*P_dens, gp)
            #loglog(Ri, nu_dat, 'r.-')

            Sig_gaia = ga.Sig_gaia(Ri, gp)
            loglog(Ri, Sig_gaia, 'b.-')

            nu_gaia = ga.rho_gaia(Ri, gp)[1]
            #loglog(Ri, nu_gaia, 'b.-')

            pdb.set_trace()

if __name__ == '__main__':
    gpr.showplots = True
    import gl_params
    gp = gl_params.Params()
    import gl_analytic as ga
    run(gp)
