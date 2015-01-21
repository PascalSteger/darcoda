#!/usr/bin/env ipython3

##
# @file
# calculate surface mass density falloff of circular rings around center of mass
# calculate velocity dispersion of 2D rings from an observed dwarf dataset
# calculate 4th order velocity moment (kurtosis) of 2D rings
# do this consistently with always the same sets of particles per bin

# (c) GPL v3 2014 Pascal S.P. Steger

import pdb
import numpy as np
from scipy.stats import kurtosis

import gi_file as gf
import gi_helper as gh
import gi_project as gip
from BiWeight import meanbiweight

def run(gp):
    import gr_params
    gpr = gr_params.grParams(gp)
    xall,yall = np.loadtxt(gp.files.get_com_file(0), skiprows=1, \
                           usecols=(0,1), unpack=True)
    # 2*[Rscale0]

    R = np.sqrt(xall**2+yall**2) # [Rscale0]
    # set number and size of (linearly spaced) bins
    Rmin = 0. #[Rscale0]
    Rmax = max(R) if gp.maxR < 0 else 1.0*gp.maxR # [Rscale0]
    R = R[(R<Rmax)] # [Rscale0]

    Binmin, Binmax, Rbin = gh.determine_radius(R, Rmin, Rmax, gp) # [Rscale0]
    gp.xipol = Rbin
    minr = min(Rbin)                           # [pc]
    maxr = max(Rbin)                           # [pc]
    gp.xepol = np.hstack([minr/8., minr/4., minr/2.,\
                          Rbin, \
                          2*maxr, 4*maxr, 8*maxr]) # [pc]
    Vol = gh.volume_circular_ring(Binmin, Binmax, gp) # [Rscale0^2]

    Rscale0 = gf.read_Xscale(gp.files.get_scale_file(0)) # [pc]

    pop=0
    print('#######  working on component ',pop)
    print('input: ', gp.files.get_com_file(pop))
    # start from data centered on COM already:
    if gf.bufcount(gp.files.get_com_file(pop))<2:
        return

    # only read in data if needed: pops = 1: reuse data from pop=0 part
    x,y = np.loadtxt(gp.files.get_com_file(pop),\
                     skiprows=1,usecols=(0,1),unpack=True)
        # [Rscalei], [Rscalei]

    # calculate 2D radius on the skyplane
    R = np.sqrt(x**2+y**2) #[Rscalei]
    Rscalei = gf.read_Xscale(gp.files.get_scale_file(pop)) # [pc]

    # set maximum radius (if gp.maxR is set)
    Rmax = max(R) if gp.maxR<0 else 1.0*gp.maxR # [Rscale0]
    print('Rmax [Rscale0] = ', Rmax)
    sel = (R * Rscalei <= Rmax * Rscale0)
    x = x[sel] # [Rscalei]
    y = y[sel] # [Rscalei]
    R = R[sel] # [Rscalei]
    totmass_tracers = float(len(x)) # [Munit], Munit = 1/star

    Rs = R                   # + possible starting offset, [Rscalei]

    tr = open(gp.files.get_ntracer_file(pop),'w')
    print(totmass_tracers, file=tr)
    tr.close()

    f_Sig, f_nu, f_mass, f_sig, f_kap, f_zeta = gf.write_headers_2D(gp, pop)

    Sig_phot   = np.zeros((gp.nipol, gpr.n))
    # particle selections, shared by density, siglos, kappa and zeta calculations
    tpb       = np.zeros((gp.nipol,gpr.n))
    for k in range(gpr.n):
        Rsi   = gh.add_errors(Rs,   gpr.Rerr)   # [Rscalei]
        for i in range(gp.nipol):
            ind1 = np.argwhere(np.logical_and(Rsi * Rscalei >= Binmin[i] * Rscale0, \
                                          Rsi * Rscalei <  Binmax[i] * Rscale0)).flatten() # [1]
            tpb[i][k] = float(len(ind1)) #[1]
            Sig_phot[i][k] = float(len(ind1))*totmass_tracers/Vol[i] # [Munit/rscale^2]

    # do the following for all populations
    Sig0 = np.sum(Sig_phot[0])/float(gpr.n) # [Munit/Rscale^2]
    Sig0pc = Sig0/Rscale0**2              # [munis/pc^2]
    gf.write_Sig_scale(gp.files.get_scale_file(pop), Sig0pc, totmass_tracers)



    # calculate density and mass profile, store it
    # ----------------------------------------------------------------------
    P_dens  = np.zeros(gp.nipol)
    P_edens = np.zeros(gp.nipol)
    for b in range(gp.nipol):
        Sig = np.sum(Sig_phot[b])/(1.*gpr.n) # [Munit/Rscale^2]
        tpbb   = np.sum(tpb[b])/float(gpr.n)       # [1], mean number of tracers in bin
        Sigerr = Sig/np.sqrt(tpbb)       # [Munit/Rscale^2], Poissonian error
        # compare data and analytic profile <=> get stellar
        # density or mass ratio from Matt Walker
        if(np.isnan(Sigerr)):
            P_dens[b] = P_dens[b-1]  # [1]
            P_edens[b]= P_edens[b-1] # [1]
        else:
            P_dens[b] = Sig/Sig0   # [1]
            P_edens[b]= Sigerr/Sig0 # [1]

        print(Rbin[b], Binmin[b], Binmax[b], P_dens[b], P_edens[b], file=f_Sig)
        # 3*[rscale], [dens0], [dens0]
        indr = (R<Binmax[b])
        Menclosed = float(np.sum(indr))/totmass_tracers # for normalization to 1#[totmass_tracers]
        Merr = Menclosed/np.sqrt(tpbb) # or artificial Menclosed/10 #[totmass_tracers]
        print(Rbin[b], Binmin[b], Binmax[b], Menclosed, Merr, file=f_mass) # [Rscale0], 2* [totmass_tracers]
    f_Sig.close()
    f_mass.close()


    # deproject Sig to get nu
    numedi = gip.Sig_INT_rho(Rbin*Rscalei, Sig0pc*P_dens, gp)
    #numin  = gip.Sig_INT_rho(Rbin*Rscalei, Sig0pc*(P_dens-P_edens), gp)
    numax  = gip.Sig_INT_rho(Rbin*Rscalei, Sig0pc*(P_dens+P_edens), gp)

    nu0pc  = numedi[0]
    gf.write_nu_scale(gp.files.get_scale_file(pop), nu0pc)

    nuerr  = numax-numedi
    for b in range(gp.nipol):
        print(Rbin[b], Binmin[b], Binmax[b],\
              numedi[b]/nu0pc, nuerr[b]/nu0pc, \
              file = f_nu)
    f_nu.close()

    # write dummy sig scale, not to be used later on
    maxsiglos = -1. #[km/s]
    fpars = open(gp.files.get_scale_file(pop),'a')
    print(maxsiglos, file=fpars)          #[km/s]
    fpars.close()

if __name__ == '__main__':
    import gi_params
    gp = gi_params.Params()
    run(gp)
