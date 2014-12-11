#!/usr/bin/env python3

##
# @file
# calculate surface mass density falloff of circular rings around center of mass
# calculate velocity dispersion of 2D rings from a Walker dataset
# calculate 4th order velocity moment (kurtosis) of 2D rings from a Walker dataset
# do this consistently with always the same sets of particles per bin
# 3D version, see grw_MCMCbin for 2D version only

# (c) GPL v3 2014 Pascal S.P. Steger

import pdb
import numpy as np
#from scipy.stats import kurtosis
#from pylab import *
#ion()

import gl_file as gf
#from gl_helper import bin_r_linear, bin_r_log
#from gl_helper import bin_r_const_tracers
# TODO determine missing modules
import gl_class_files as gcf

# from BiWeight import meanbiweight
import gl_analytic as ga


def volume_spherical_shell(binmin, binmax, gp):
    vol = np.zeros(gp.nipol)
    for k in range(gp.nipol):
        vol[k] = 4.*np.pi/3.*(binmax[k]**3 - binmin[k]**3)
        # [rscale^3]
    return vol
## \fn volume_spherical_shell(binmin, binmax, gp)
# volume of a circular ring from binmin to binmax


def run(gp):
    import gr_params
    gpr = gr_params.grParams(gp)
    xall,yall,zall = np.loadtxt(gp.files.get_com_file(0),skiprows=1,\
                                usecols=(0,1,2),unpack=True) # 2*[rscale0]
    rscale0 = gf.read_Xscale(gp.files.get_scale_file(0)+'_3D')
    xall *= rscale0
    yall *= rscale0
    zall *= rscale0

    # calculate 3D
    r = np.sqrt(xall**2+yall**2+zall**2) #[pc]

    # set number and size of (linearly spaced) bins
    rmin = 0. # [pc]
    rmax = max(r) if gp.maxR < 0 else 1.0*gp.maxR           # [pc]
    print('rmax [rscale] = ', rmax)
    r = r[(r<rmax)] # [pc]

    binmin, binmax, rbin = gh.determine_radius(r, rmin, rmax, gp) # [pc]
    vol = volume_spherical_shell(binmin, binmax, gp) # [pc^3]

    for pop in range(gp.pops+1):
        print('#######  working on component ',pop)
        print('input: ',gp.files.get_com_file(pop)+'_3D')
        # start from data centered on COM already:
        if gf.bufcount(gp.files.get_com_file(pop)+'_3D')<2: continue
        x,y,z,v = np.loadtxt(gp.files.get_com_file(pop)+'_3D',\
                           skiprows=1,usecols=(0,1,2,3),unpack=True)
        # 3*[rscale], [km/s]
        rscalei = gf.read_Xscale(gp.files.get_scale_file(pop)) # [pc]
        x *= rscalei
        y *= rscalei
        z *= rscalei
        # calculate 2D radius on the skyplane
        r = np.sqrt(x**2+y**2+z**2) # [pc]

        # set maximum radius (if gp.maxR is set)
        rmax = max(r) if gp.maxR<0 else 1.0*gp.maxR # [pc]
        print('rmax [pc] = ', rmax)
        sel = (r<=rmax)
        x = x[sel]; y = y[sel]; z = z[sel]; v = v[sel]; r = r[sel] # [rscale]
        totmass_tracers = 1.*len(x) # [Munit], Munit = 1/star

        rs = r                   # + possible starting offset, [rscale]
        vlos = v                 # + possible starting offset, [km/s]

        gf.write_tracer_file(gp.files.get_ntracer_file(pop)+'_3D', totmass_tracers)
        de, em = gf.write_headers_3D(gp, pop)

        # gpr.n=30 iterations for getting random picked radius values
        density = np.zeros((gp.nipol,gpr.n))
        a       = np.zeros((gp.nipol,gpr.n)) # shared by density, siglos, kappa calcs
        for k in range(gpr.n):
            rsi = gpr.Rerr * np.random.randn(len(rs)) + rs # [pc]
            vlosi = gpr.vrerr*np.random.randn(len(vlos)) + vlos # [km/s]
            for i in range(gp.nipol):
                ind1 = np.argwhere(np.logical_and(rsi>=binmin[i], rsi<binmax[i])).flatten() # [1]
                density[i][k] = (1.*len(ind1))/vol[i]*totmass_tracers # [Munit/rscale^2]
                vlos1 = vlosi[ind1]                           # [km/s]
                a[i][k] = 1.*len(ind1)                        # [1]


        dens0 = np.sum(density[0])/(1.*gpr.n) # [Munit/rscale^3]
        print('dens0 = ',dens0,' [Munit/rscale^3]')

        dens0pc = dens0/rscale0**3
        gf.write_Sig_scale(gp.files.get_scale_file(pop)+'_3D', dens0pc, totmass_tracers)

        tpb0   = np.sum(a[0])/float(gpr.n)     # [1] tracers per bin
        denserr0 = dens0/np.sqrt(tpb0)       # [Munit/rscale^3]
        p_dens  = np.zeros(gp.nipol)
        p_edens = np.zeros(gp.nipol)
        for b in range(gp.nipol):
            dens = np.sum(density[b])/float(gpr.n) # [Munit/rscale^3]
            tpb  = np.sum(a[b])/float(gpr.n)       # [1]
            denserr = dens/np.sqrt(tpb)            # [Munit/rscale^3]

            if(np.isnan(denserr)):
                p_dens[b] = p_dens[b-1]                          # [1]
                p_edens[b]= p_edens[b-1]                         # [1]
            else:
                p_dens[b] = dens/dens0                           # [1]
                p_edens[b]= denserr/dens0 # [1] #100/rbin would be artificial guess

            print(rbin[b], binmin[b], binmax[b], p_dens[b], p_edens[b], file=de)
            # [rscale], 2*[dens0]
            indr = (r<binmax[b])
            menclosed = float(np.sum(indr))/totmass_tracers # for normalization to 1
            # [totmass_tracers]
            merr = menclosed/np.sqrt(tpb) # artificial menclosed/10 # [totmass_tracers]
            print(rbin[b], binmin[b], binmax[b], menclosed, merr, file=em)
            # [rscale], 2*[totmass_tracers]
        de.close()
        em.close()

        if gpr.showplots:
            print('plotting for pop ', pop)
            #show_plots_dens(rbin, p_dens, p_edens, gp)
            mf1 = 0.02 #1/totmass_tracers
            mf2 = 0.02
            rho_dm, rho_star1, rho_star2 = ga.rho_walk(rbin*rscale0, gp, mf1, mf2)

            if pop == 0:
                loglog(rbin*rscale0, rho_star1+rho_star2, 'k.-', lw=0.5)
            elif pop == 1:
                loglog(rbin*rscale0, rho_star1, 'b.-', lw = 0.5)
            elif pop == 2:
                loglog(rbin*rscale0, rho_star2, 'g.-', lw = 0.5)

            loglog(rbin*rscale0, dens0pc*p_dens, 'r.-')
            pdb.set_trace()
            clf()

if __name__ == '__main__':
    gpr.showplots = True
    import gl_params
    gp = gl_params.Params()
    run(gp)
