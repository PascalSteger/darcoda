#!/usr/bin/env ipython3

##
# @file
# calculate surface mass density falloff of circular rings around center of mass
# for triaxial systems
# TODO: enable ellipsoidal bins

# (c) 2013 Pascal S.P. Steger


import sys, ipdb
import numpy as np
import gr_params as gpr
import gl_file as gf
import gl_helper as gh
import gl_project as glp

def run(gp):
    Rscale0 = gf.read_Xscale(gp.files.get_scale_file(0)) # [pc]
    print('input: ',gpr.get_com_file(0))
    # start from data centered on COM already:
    x,y,v = np.loadtxt(gpr.get_com_file(0),\
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
        totmass = 1.*len(x) #[munit], munit = 1/star

        Binmin, Binmax, Rbin = gpr.determine_radius(R, Rmin, Rmax, gp) # [Rscale0]
        Vol = gpr.volume_circular_ring(Binmin, Binmax, gp)

        # rs = gpr.Rerr*np.random.randn(len(r))+r
        Rs = R  # [Rscale] # if no initial offset is whished

        tr = open(gp.files.get_ntracer_file(0),'w')
        print(totmass, file=tr)
        tr.close()

        f_Sig, f_nu, f_mass, f_sig, f_kap = gf.write_headers_2D(gp, 0)

        # 30 iterations for getting random picked radius values
        Density = np.zeros((gp.nipol,gpr.n))
        tpb       = np.zeros((gp.nipol,gpr.n))
        for k in range(gpr.n):
            Rsi = gh.add_errors(Rs, gpr.Rerr) # [Rscalei]
            for j in range(gp.nipol):
                ind1 = np.argwhere(np.logical_and(Rsi * Rscalei >= Binmin[j] * Rscale0,\
                                                  Rsi * Rscalei <  Binmax[j] * Rscale0)).flatten() # [1]
                Density[j][k] = float(len(ind1))/Vol[j]*totmass # [munit/Rscale0^2]
                tpb[j][k] = float(len(ind1)) #[1]

        Dens0 = np.sum(Density[0])/float(gpr.n) # [Munit/Rscale0^2]
        Dens0pc = Dens0/Rscale0**2 # [Munit/pc^2]
        gf.write_Sig_scale(gp.files.get_scale_file(0), Dens0pc, totmass)

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
            menclosed = float(np.sum(indr))/totmass
             # /totmass for normalization to 1 at last bin #[totmass]
            merr = menclosed/np.sqrt(tpbb) # artificial menclosed/10 gives good approximation #[totmass]
            print(Rbin[b], Binmin[b], Binmax[b], menclosed, merr, file=f_mass)
            # [rscale], [totmass], [totmass]
        de.close()
        em.close()

        # deproject Sig to get nu
        numedi = glp.Sig_INT_rho(Rbin*Rscalei, Dens0pc*P_dens, gp)
        numin  = glp.Sig_INT_rho(Rbin*Rscalei, Dens0pc*(P_dens-P_edens), gp)
        numax  = glp.Sig_INT_rho(Rbin*Rscalei, Dens0pc*(P_dens+P_edens), gp)

        nu0pc  = numedi[0]
        gf.write_nu_scale(gp.files.get_scale_file(pop), nu0pc)

        nuerr  = numax-numedi
        for b in range(gp.nipol):
            print(Rbin[b], Binmin[b], Binmax[b],\
                  numedi[b]/nu0pc, nuerr[b]/nu0pc, \
                  file = f_nu)
        f_nu.close()


        if gpr.showplots:
            gpr.show_plots_dens_2D(0, Rbin, p_dens, p_edens, Dens0pc)


if __name__ == '__main__':
    gpr.showplots = True
    import gl_params
    gp = gl_params.Params()
    run(gp)
