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

def obs_Sig_phot(Binmin, Binmax, Rscale0, Sig_kin, gp, gpr):
    Sig_phot   = np.zeros((gp.nipol, gpr.n))
    for kbin in range(gp.nipol):
        Rw = (Binmax[kbin]-Binmin[kbin])*Rscale0 # [pc]
        kpc = 1000 # [pc]
        DL = {1: lambda x: x * (138),#+/- 8 for Fornax
              2: lambda x: x * (101),#+/- 5 for Carina
              3: lambda x: x * (79), #+/- 4 for Sculptor
              4: lambda x: x * (86), #+/- 4 for Sextans
              5: lambda x: x * (80)  #+/- 10 for Draco
        }[gp.case](kpc)
        arcmin = 2.*np.pi* DL / (360 * 60) # [pc] at distance 138 kpc for Fornax
        width = (Binmax[kbin]-Binmin[kbin])*Rscale0/arcmin # [arcmin]
        # change scale according to width
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
        Sig_phot[kbin] = Sig_kin[kbin] / w_ipol
    return Sig_phot
## \fn obs_Sig_phot(Binmin, Binmax, Rscale0, Sig_kin, gp)
# return photometric surface density for observations with selection function w
# only in 'obs' investigation
# @param Binmin [pc]
# @param Binmax [pc]
# @param Rscale0 [pc]
# @param Sig_kin [Msun/pc^2]
# @param gp global parameters

def run(gp):
    import gr_params
    gpr = gr_params.grParams(gp)
    xall,yall = np.loadtxt(gp.files.get_com_file(0), skiprows=1, usecols=(0,1), unpack=True)
    # 2*[Rscale0]
    R = np.sqrt(xall**2+yall**2) # [Rscale0]
    # set number and size of (linearly spaced) bins
    Rmin = 0. #[Rscale0]
    Rmax = max(R) if gp.maxR < 0 else 1.0*gp.maxR # [Rscale0]
    R = R[(R<Rmax)] # [Rscale0]
    Binmin, Binmax, Rbin = gh.determine_radius(R, Rmin, Rmax, gp) # [Rscale0]
    gp.xipol = Rbin
    minr = min(Rbin) # [pc]
    maxr = max(Rbin) # [pc]
    gp.xepol = np.hstack([minr/8., minr/4., minr/2., Rbin, 2*maxr, 4*maxr, 8*maxr]) # [pc]
    Vol = gh.volume_circular_ring(Binmin, Binmax, gp) # [Rscale0^2]
    Rscale0 = gf.read_Xscale(gp.files.get_scale_file(0)) # [pc]
    for pop in range(gp.pops+1):
        print('#######  working on component ',pop)
        print('input: ', gp.files.get_com_file(pop))
        # exclude second condition if self-consistent approach wished
        if gp.investigate == "obs" and gp.case==1 and pop==0:
            # for Fornax, overwrite first Sigma with deBoer data
            import gr_MCMCbin_for
            gr_MCMCbin_for.run(gp)
            continue
        # start from data centered on COM already:
        if gf.bufcount(gp.files.get_com_file(pop))<2:
            continue
        # only read in data if needed: pops = 1: reuse data from pop=0 part
        if (gp.pops == 1 and pop < 1 or gp.pops == 2) or gp.investigate == 'obs':
            x,y,v = np.loadtxt(gp.files.get_com_file(pop), skiprows=1,usecols=(0,1,2),unpack=True)
            # [Rscalei], [Rscalei], [km/s]
            # calculate 2D radius on the skyplane
            R = np.sqrt(x**2+y**2) #[Rscalei]
            Rscalei = gf.read_Xscale(gp.files.get_scale_file(pop)) # [pc]
            # set maximum radius (if gp.maxR is set)
            Rmax = max(R) if gp.maxR<0 else 1.0*gp.maxR # [Rscale0]
            print('Rmax [Rscale0] = ', Rmax)
            sel = (R * Rscalei <= Rmax * Rscale0)
            x = x[sel]
            y = y[sel]
            v = v[sel]
            R = R[sel] # [Rscalei]
            totmass_tracers = float(len(x)) # [Munit], Munit = 1/star
            Rs = R                   # + possible starting offset, [Rscalei]
            vlos = v                 # + possible starting offset, [km/s]
        tr = open(gp.files.get_ntracer_file(pop),'w')
        print(totmass_tracers, file=tr)
        tr.close()
        f_Sig, f_nu, f_mass, f_sig, f_kap, f_zeta = gf.write_headers_2D(gp, pop)
        if (gp.pops == 1 and pop < 1) or gp.pops == 2 or gp.investigate == 'obs':
            Sig_kin   = np.zeros((gp.nipol, gpr.n))
            siglos    = np.zeros((gp.nipol, gpr.n))
            if gp.usekappa:
                kappa     = np.zeros((gp.nipol, gpr.n))
            if gp.usezeta:
                v2        = np.zeros((gp.nipol, gpr.n))
                v4        = np.zeros((gp.nipol, gpr.n))
                Ntot      = np.zeros(gpr.n)
                zetaa     = np.zeros(gpr.n)
                zetab     = np.zeros(gpr.n)
            # particle selections, shared by density, siglos, kappa and zeta calculations
            tpb       = np.zeros((gp.nipol,gpr.n))
            for k in range(gpr.n):
                Rsi   = gh.add_errors(Rs,   gpr.Rerr)   # [Rscalei]
                vlosi = gh.add_errors(vlos, gpr.vrerr)   # [km/s]
                for i in range(gp.nipol):
                    ind1 = np.argwhere(np.logical_and(Rsi * Rscalei >= Binmin[i] * Rscale0, Rsi * Rscalei <  Binmax[i] * Rscale0)).flatten() # [1]
                    tpb[i][k] = float(len(ind1)) #[1]
                    Sig_kin[i][k] = float(len(ind1))*totmass_tracers/Vol[i] # [Munit/rscale**2]
                    if(len(ind1)<=1):
                        siglos[i][k] = siglos[i-1][k]
                        print('### using last value, missing data')
                        if gp.usekappa:
                            kappa[i][k] = kappa[i-1][k]
                            # attention! should be 0, uses last value
                        if gp.usezeta:
                            v2[i][k] = v2[i-1][k]
                            v4[i][k] = v4[i-1][k]
                    else:
                        siglos[i][k] = meanbiweight(vlosi[ind1], ci_perc=68.4, \
                                                    ci_mean=True, ci_std=True)[1]
                        # [km/s], see BiWeight.py
                        if gp.usekappa:
                            kappa[i][k] = kurtosis(vlosi[ind1], axis=0, \
                                                   fisher=False, bias=False) # [1]
                        if gp.usezeta:
                            ave, adev, sdev, var, skew, curt = gh.moments(vlosi[ind1])
                            v2[i][k] = var
                            v4[i][k] = (curt+3)*var**2
                Sigma = Sig_kin[:,k]
                if gp.usezeta:
                    pdb.set_trace()
                    Ntot[k] = gh.Ntot(Rbin, Sigma, gp)
                    zetaa[k] = gh.starred(Rbin, v4[:,k], Sigma, Ntot[k], gp)
                    v2denom = (gh.starred(Rbin, v2[:,k], Sigma, Ntot[k], gp))**2
                    zetaa[k] /= v2denom
                    zetab[k] = gh.starred(Rbin, v4[:,k]*Rbin**2, Sigma, Ntot[k], gp)
                    zetab[k] /= v2denom
                    zetab[k] /= (gh.starred(Rbin, Rbin, Sigma, Ntot[k], gp))**2
            if gp.investigate == 'obs' and gp.case < 5:
                Sig_phot = obs_Sig_phot(Binmin, Binmax, Rscale0, Sig_kin, gp, gpr)
            else:
                Sig_phot = Sig_kin
        # do the following for all populations
        Sig0 = np.sum(Sig_phot[0])/float(gpr.n) # [Munit/Rscale^2]
        Sig0pc = Sig0/Rscale0**2              # [munis/pc^2]
        gf.write_Sig_scale(gp.files.get_scale_file(pop), Sig0pc, totmass_tracers)
        # calculate density and mass profile, store it
        # ----------------------------------------------------------------------
        #tpb0   = np.sum(tpb[0])/float(gpr.n)     # [1]
        #Sigerr0 = Sig0/np.sqrt(tpb0)       # [Munit/Rscale^2]
        P_dens  = np.zeros(gp.nipol)
        P_edens = np.zeros(gp.nipol)
        for b in range(gp.nipol):
            Sig = np.sum(Sig_kin[b])/(1.*gpr.n) # [Munit/Rscale^2]
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
            print(Rbin[b], Binmin[b], Binmax[b], numedi[b]/nu0pc, nuerr[b]/nu0pc, file = f_nu)
        f_nu.close()
        # calculate and output siglos
        # --------------------------------------------
        p_dvlos = np.zeros(gp.nipol)
        p_edvlos = np.zeros(gp.nipol)
        for b in range(gp.nipol):
            sig = np.sum(siglos[b])/gpr.n #[km/s]
            tpbb = np.sum(tpb[b])/float(gpr.n) #[1]
            if tpbb == 0:
                sigerr = p_edvlos[b-1] #[km/s]
                # attention! uses last error
            else:
                sigerr = sig/np.sqrt(tpbb) #[km/s]
            p_dvlos[b] = sig    #[km/s]
            p_edvlos[b]= sigerr #[km/s]
        maxsiglos = max(p_dvlos) #[km/s]
        print('maxsiglos = ', maxsiglos, '[km/s]')
        fpars = open(gp.files.get_scale_file(pop),'a')
        print(maxsiglos, file=fpars)          #[km/s]
        fpars.close()
        for b in range(gp.nipol):
            print(Rbin[b], Binmin[b], Binmax[b], np.abs(p_dvlos[b]/maxsiglos),\
                  np.abs(p_edvlos[b]/maxsiglos), file=f_sig)
            # 3*[rscale], 2*[maxsiglos]
        f_sig.close()
        # calculate and output kurtosis kappa
        # --------------------------------------------
        if gp.usekappa:
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
                      kappavel, kappavelerr, file=f_kap)
                # [rscale], 2*[1]
            f_kap.close()
        # output zetas
        # -------------------------------------------------------------
        if gp.usezeta:
            print(np.median(zetaa), np.median(zetab), file=f_zeta)
            f_zeta.close()
        if gpr.showplots:
            gpr.show_plots_dens_2D(Rbin*Rscalei, P_dens, P_edens, Sig0pc)
            gpr.show_plots_sigma(Rbin*Rscalei, p_dvlos, p_edvlos)
            if gp.usekappa:
                gpr.show_plots_kappa(Rbin*Rscalei, p_kappa, p_ekappa)

        # overwrite Sig profile if photometric data is used
        if gp.investigate == 'obs' and gp.case==1 and pop==1 and not gp.selfconsistentnu:
            import os
            os.system('cp '+gp.files.get_scale_file(0)+' '+gp.files.get_scale_file(1))
            # replace last line with actual maxsiglos from tracer particles
            os.system("sed -i '$s/^.*/"+str(maxsiglos)+"/' "+gp.files.get_scale_file(1))
            os.system('cp '+gp.files.Sigfiles[0]+' '+gp.files.Sigfiles[1])
            continue

if __name__ == '__main__':
    import gi_params
    gp = gi_params.Params()
    run(gp)
