#!/usr/bin/env ipython3

##
# @file
# calculate surface mass density falloff of circular rings around center of mass
# for triaxial systems

# (c) 2015 GPL v.3 Pascal S.P. Steger, pascal@steger.aero

import pdb
import numpy as np
import gi_file as gf
import gi_helper as gh
import gi_project as gip

def run(gp):
    import gr_params
    gpr = gr_params.grParams(gp)
    print('scalefile: ', gp.files.get_scale_file(0))
    Rscale0 = gf.read_Xscale(gp.files.get_scale_file(0)) # [pc]
    print('input: ',gp.files.get_com_file(0))
    # start from data centered on COM already:
    x,y,v = np.loadtxt(gp.files.get_com_file(0),\
                       skiprows=1,usecols=(0,1,2),unpack=True) #[Rscalei], [Rscalei], [km/s]

    for pop in range(2):
        # calculate 2D radius on the skyplane
        R = np.sqrt(x**2+y**2) # [Rscalei]
        Rscalei = gf.read_Xscale(gp.files.get_scale_file(pop)) # [pc]
        # set number and size of bins
        Rmin = 0. # [rscale]
        Rmax = max(R) if gp.maxR < 0 else float(gp.maxR)   # [Rscale0]


        sel = (R * Rscalei < Rmax * Rscale0)
        x = x[sel]; y = y[sel]; v = v[sel] #[rscale]
        totmass_tracers = 1.*len(x) #[munit], munit = 1/star

        Binmin, Binmax, Rbin = gh.determine_radius(R, Rmin, Rmax, gp) # [Rscale0]
        gp.xipol = Rbin
        minr = min(Rbin)                           # [pc]
        maxr = max(Rbin)                           # [pc]
        gp.xepol = np.hstack([minr/8., minr/4., minr/2.,\
                              Rbin, \
                              2*maxr, 4*maxr, 8*maxr]) # [pc]
        Vol = gh.volume_circular_ring(Binmin, Binmax, gp) # [Rscale0^2]

        # rs = gpr.Rerr*np.random.randn(len(r))+r
        Rs = R  # [Rscale] # if no initial offset is whished

        tr = open(gp.files.get_ntracer_file(pop),'w')
        print(totmass_tracers, file=tr)
        tr.close()

        f_Sig, f_nu, f_mass, f_sig, f_kap, f_zeta = gf.write_headers_2D(gp, pop)

        # 30 iterations for getting random picked radius values
        Density = np.zeros((gp.nipol,gpr.n))
        tpb       = np.zeros((gp.nipol,gpr.n))
        for k in range(gpr.n):
            Rsi = gh.add_errors(Rs, gpr.Rerr) # [Rscalei]
            for j in range(gp.nipol):
                ind1 = np.argwhere(np.logical_and(Rsi * Rscalei >= Binmin[j] * Rscale0,\
                                                  Rsi * Rscalei <  Binmax[j] * Rscale0)).flatten() # [1]
                Density[j][k] = float(len(ind1))/Vol[j]*totmass_tracers # [munit/Rscale0^2]
                tpb[j][k] = float(len(ind1)) #[1]

        Dens0 = np.sum(Density[0])/float(gpr.n) # [Munit/Rscale0^2]
        Dens0pc = Dens0/Rscale0**2 # [Munit/pc^2]
        gf.write_Sig_scale(gp.files.get_scale_file(pop), Dens0pc, totmass_tracers)

        tpbb0   = np.sum(tpb[0])/float(gpr.n)     # [1]
        Denserr0 = Dens0/np.sqrt(tpbb0)       # [Munit/rscale^2]

        p_dens  = np.zeros(gp.nipol)
        p_edens = np.zeros(gp.nipol)

        for b in range(gp.nipol):
            Dens = np.sum(Density[b])/float(gpr.n) # [Munit/rscale^2]
            tpbb = np.sum(tpb[b])/float(gpr.n)       # [1]
            Denserr = Dens/np.sqrt(tpbb)       # [Munit/rscale^2]
            if(np.isnan(Denserr)):
                p_dens[b] = p_dens[b-1]  # [1]
                p_edens[b]= p_edens[b-1] # [1]
            else:
                p_dens[b] = Dens/Dens0   # [1]
                p_edens[b]= Denserr/Dens0    # [1] #100/rbin would be artificial guess

        for b in range(gp.nipol):
            print(Rbin[b], Binmin[b], Binmax[b], p_dens[b], p_edens[b], file=f_Sig)
            # [rscale], [dens0], [dens0]
            indr = (R < Binmax[b])
            menclosed = float(np.sum(indr))/totmass_tracers
             # /totmass_tracers for normalization to 1 at last bin #[totmass_tracers]
            merr = menclosed/np.sqrt(tpbb) # artificial menclosed/10 gives good approximation #[totmass_tracers]
            print(Rbin[b], Binmin[b], Binmax[b], menclosed, merr, file=f_mass)
            # [rscale], [totmass_tracers], [totmass_tracers]
        f_Sig.close()
        f_mass.close()

        # deproject Sig to get nu
        numedi = gip.Sig_INT_rho(Rbin*Rscalei, Dens0pc*p_dens, gp)
        numin  = gip.Sig_INT_rho(Rbin*Rscalei, Dens0pc*(p_dens-p_edens), gp)
        numax  = gip.Sig_INT_rho(Rbin*Rscalei, Dens0pc*(p_dens+p_edens), gp)

        nu0pc  = numedi[0]
        gf.write_nu_scale(gp.files.get_scale_file(pop), nu0pc)

        nuerr  = numax-numedi
        for b in range(gp.nipol):
            print(Rbin[b], Binmin[b], Binmax[b],\
                  numedi[b]/nu0pc, nuerr[b]/nu0pc, \
                  file = f_nu)
        f_nu.close()

        if gpr.showplots:
            gpr.show_plots_dens_2D(Rbin*Rscalei, p_dens, p_edens, Dens0pc)


if __name__ == '__main__':
    gpr.showplots = True
    import gi_params
    gp = gi_params.Params()
    run(gp)
