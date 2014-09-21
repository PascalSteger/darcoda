#!/usr/bin/env ipython3

##
# @file
# calculate surface mass density falloff of circular rings around center of mass
# calculate velocity dispersion of 2D rings from an observed dwarf dataset
# calculate 4th order velocity moment (kurtosis) of 2D rings
# do this consistently with always the same sets of particles per bin

# (c) 2013 Pascal S.P. Steger

import sys, pdb
import numpy as np
from scipy.stats import kurtosis
from pylab import *

import gr_params as gpr
import gl_file as gfile
import gl_helper as gh
from BiWeight import meanbiweight


def run(gp):
    xall,yall = np.loadtxt(gpr.get_com_file(0), skiprows=1, \
                           usecols=(0,1), unpack=True)
    # 2*[Rscale0]

    R = np.sqrt(xall**2+yall**2) # [Rscale0]
    # set number and size of (linearly spaced) bins
    Rmin = 0. #[Rscale0]
    Rmax = max(R) if gp.maxR < 0 else 1.0*gp.maxR # [Rscale0]
    R = R[(R<Rmax)] # [Rscale0]

    Binmin, Binmax, Rbin = gpr.determine_radius(R, Rmin, Rmax, gp) # [Rscale0]
    Vol = gpr.volume_circular_ring(Binmin, Binmax, gp) # [Rscale0^2]

    Rscale0 = gfile.read_Xscale(gp.files.get_scale_file(0)) # [pc]

    for pop in range(gpr.pops):
        print('#######  working on component ',pop)
        print('input: ', gpr.get_com_file(pop))
        # start from data centered on COM already:
        if gfile.bufcount(gpr.get_com_file(pop))<2: continue
        x,y,v = np.loadtxt(gpr.get_com_file(pop),\
                           skiprows=1,usecols=(0,1,2),unpack=True)
                           # [Rscalei], [Rscalei], [km/s]

        # calculate 2D radius on the skyplane
        R = np.sqrt(x**2+y**2) #[Rscalei]
        Rscalei = gfile.read_Xscale(gp.files.get_scale_file(pop)) # [pc]

        # set maximum radius (if gp.maxR is set)
        Rmax = max(R) if gp.maxR<0 else 1.0*gp.maxR # [Rscale0]
        print('Rmax [Rscale0] = ', Rmax)
        sel = (R * Rscalei <= Rmax * Rscale0)
        x = x[sel]; y = y[sel]; v = v[sel]; R = R[sel] # [Rscalei]
        totmass = float(len(x)) # [Munit], Munit = 1/star
        
        Rs = R                   # + possible starting offset, [Rscalei]
        vlos = v                 # + possible starting offset, [km/s]
        
        tr = open(gp.files.get_ntracer_file(pop),'w')
        print(totmass, file=tr)
        tr.close()

        de, em, sigfil, kappafil = gfile.write_headers_2D(gp, pop)

        Density_kin   = np.zeros((gp.nipol, gpr.n))
        sigma     = np.zeros((gp.nipol, gpr.n))
        kappa     = np.zeros((gp.nipol, gpr.n))
        zetaa     = np.zeros((gp.nipol, gpr.n))
        zetab     = np.zeros((gp.nipol, gpr.n))
        # shared by density, siglos, kappa c alcs
        tpb       = np.zeros((gp.nipol,gpr.n))
        for k in range(gpr.n):
            Rsi   = gh.add_errors(Rs,   gpr.Rerr)   # [Rscalei]
            vlosi = gh.add_errors(vlos, gpr.vrerr)   # [km/s]
            for i in range(gp.nipol):
                ind1 = np.argwhere(np.logical_and(Rsi * Rscalei >= Binmin[i] * Rscale0, \
                                                  Rsi * Rscalei <  Binmax[i] * Rscale0)).flatten() # [1]
                tpb[i][k] = float(len(ind1)) #[1]
                Density_kin[i][k] = float(len(ind1))*totmass/Vol[i]
                # [Munit/rscale**2]

                if(len(ind1)<=1):
                    sigma[i][k] = sigma[i-1][k]
                    print('### using last value, missing data')
                    kappa[i][k] = kappa[i-1][k]
                    zetaa[i][k] = zetaa[i-1][k]
                    zetab[i][k] = zetab[i-1][k]
                    # attention! should be 0, uses last value
                else:
                    sigma[i][k] = meanbiweight(vlosi[ind1], ci_perc=68.4, \
                                               ci_mean=True, ci_std=True)[1]
                                        # [km/s], see BiWeight.py
                    kappa[i][k] = kurtosis(vlosi[ind1], axis=0, \
                                           fisher=False, bias=False) # [1]
                    # zetaa[i][k] = TODO
                    # zetab[i][k] = TODO


        # output density
        for kkk in range(gp.nipol):
            Rw = (Binmax[kkk]-Binmin[kkk])*Rscale0 # [pc]
            kpc = 1000 # [pc]
            DL = {0: lambda x: x * (138),#+/- 8 for Fornax
                  1: lambda x: x * (101),#+/- 5 for Carina
                  2: lambda x: x * (79),  #+/- 4 for Sculptor
                  3: lambda x: x * (86) #+/- 4 for Sextans
              }[gp.case](kpc)
            arcmin = 2.*np.pi* DL / (360 * 60) # [pc] at distance 138 kpc for Fornax
            width = (Binmax[kkk]-Binmin[kkk])*Rscale0/arcmin # [arcmin]
            # TODO change scale according to width
            if width < 1.: # [arcmin]
                wsize_ipol = '0.5'
            elif width < 3.5: # [arcmin]
                wsize_ipol = '2.0'
            elif width < 7.5: # [arcmin]
                wsize_ipol = '5.0'
            else:          # [arcmin]
                wsize_ipol = '10.0'

            A = np.loadtxt(gp.files.dir+'w_'+wsize_ipol+'.dat')

            Rpt, wpt = A.T # [arcmin], [1]
            Rpt *= arcmin # [pc]
            w_ipol = wpt[np.where(abs(Rw-Rpt) == min(abs(Rw-Rpt)))]
            # collapsing to value at nearest radius
            # this is not exact, but gives much faster code than
            # gh.ipol(Rpt, wpt, Rw) # all radii in [pc]
            Density_phot = Density_kin / w_ipol

        Dens0 = np.sum(Density_phot[0])/float(gpr.n) # [Munit/Rscale^2]
        Dens0pc = Dens0/Rscale0**2              # [munis/pc^2]
        gfile.write_density(gp.files.get_scale_file(pop), Dens0pc, totmass)

        tpb0   = np.sum(tpb[0])/float(gpr.n)     # [1]
        Denserr0 = Dens0/np.sqrt(tpb0)       # [Munit/Rscale^2]
        P_dens  = np.zeros(gp.nipol)
        P_edens = np.zeros(gp.nipol)
        for b in range(gp.nipol):
            Dens = np.sum(Density_kin[b])/(1.*gpr.n) # [Munit/Rscale^2]
            tpbb   = np.sum(tpb[b])/float(gpr.n)       # [1], mean number of tracers in bin
            Denserr = Dens/np.sqrt(tpbb)       # [Munit/Rscale^2], Poissonian error
            # compare data and analytic profile <=> get stellar
            # density or mass ratio from Matt Walker
            if(np.isnan(Denserr)):
                P_dens[b] = P_dens[b-1]  # [1]
                P_edens[b]= P_edens[b-1] # [1]
            else:
                P_dens[b] = Dens/Dens0   # [1]
                P_edens[b]= Denserr/Dens0 # [1]

            print(Rbin[b], Binmin[b], Binmax[b], P_dens[b], P_edens[b], file=de)
            # 3*[rscale], [dens0], [dens0]
            indr = (R<Binmax[b])
            Menclosed = float(np.sum(indr))/totmass # for normalization to 1#[totmass]
            Merr = Menclosed/np.sqrt(tpbb) # or artificial Menclosed/10 #[totmass]
            print(Rbin[b], Binmin[b], Binmax[b], Menclosed, Merr, file=em) # [Rscale0], 2* [totmass]
        de.close()
        em.close()


        # output siglos
        p_dvlos = np.zeros(gp.nipol);        p_edvlos = np.zeros(gp.nipol)
        for b in range(gp.nipol):
            dispvel = np.sum(sigma[b])/gpr.n #[km/s]
            tpbb = np.sum(tpb[b])/float(gpr.n) #[1]
            if tpbb == 0:
                dispvelerr = p_edvlos[b-1] #[km/s]
                # attention! uses last error
            else:
                dispvelerr = dispvel/np.sqrt(tpbb) #[km/s]
            p_dvlos[b] = dispvel    #[km/s]
            p_edvlos[b]= dispvelerr #[km/s]

        maxsiglos = max(p_dvlos) #[km/s]
        print('maxsiglos = ', maxsiglos, '[km/s]')
        fpars = open(gp.files.get_scale_file(pop),'a')
        print(maxsiglos, file=fpars)          #[km/s]
        fpars.close()

        for b in range(gp.nipol):
            print(Rbin[b], Binmin[b], Binmax[b], np.abs(p_dvlos[b]/maxsiglos),\
                  np.abs(p_edvlos[b]/maxsiglos), file=sigfil)
            # 3*[rscale], 2*[maxsiglos]
        sigfil.close()

        # output kurtosis kappa
        p_kappa = np.zeros(gp.nipol) # needed for plotting later
        p_ekappa = np.zeros(gp.nipol)
        for b in range(gp.nipol):
            kappavel = np.sum(kappa[b])/gpr.n #[1]
            tpbb = np.sum(tpb[b])/float(gpr.n) #[1]
            if tpbb == 0:
                kappavelerr = p_edvlos[b-1] #[1]
                # attention! uses last error
            else:
                kappavelerr = np.abs(kappavel/np.sqrt(tpbb)) #[1]
            p_kappa[b] = kappavel
            p_ekappa[b] = kappavelerr
            
            print(Rbin[b], Binmin[b], Binmax[b], \
                  kappavel, kappavelerr, file=kappafil)
            # [rscale], 2*[1]
        kappafil.close()
    
        if gpr.showplots:
            gpr.show_plots_dens_2D(pop, Rbin, P_dens, P_edens, Dens0pc)
            gpr.show_plots_sigma(pop, Rbin, p_dvlos, p_edvlos)
            gpr.show_plots_kappa(pop, Rbin, p_kappa, p_ekappa)



if __name__ == '__main__':
    gpr.showplots = True
    import gl_params
    gp = gl_params.Params()
    run(gp)

