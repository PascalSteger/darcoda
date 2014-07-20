#!/usr/bin/env ipython3

##
# @file
# calculate surface mass density falloff of circular rings around center of mass
# calculate velocity dispersion of 2D rings from a Walker dataset
# calculate 4th order velocity moment (kurtosis) of 2D rings from a Walker dataset
# do this consistently with always the same sets of particles per bin
# 3D version, see grw_MCMCbin for 2D version only

# (c) 2013 Pascal S.P. Steger

import sys, pdb
import numpy as np
from scipy.stats import kurtosis
from pylab import *

import gr_params as gpr
import gl_file as gfile
from gl_helper import expDtofloat, bin_r_linear, bin_r_log, bin_r_const_tracers
from gl_class_files import *
from BiWeight import meanbiweight

def volume_spherical_shell(Binmin, Binmax, gp):
    vol = np.zeros(gp.nipol)
    for k in range(gp.nipol):
        vol[k] = 4.*np.pi/3.*(Binmax[k]**3 - Binmin[k]**3)
        # [Rscale^3]
## \fn volume_spherical_shell(Binmin, Binmax, gp)
# volume of a circular ring from binmin to binmax


def run(gp):
    xall,yall = np.loadtxt(gpr.get_com_file(0),skiprows=1,usecols=(0,1),unpack=True) # 2*[Rscale0]
    # calculate 2D radius on the skyplane
    R = np.sqrt(xall**2+yall**2) #[Rscale0]
    # set number and size of (linearly spaced) bins
    Rmin = 0. # [Rscale0]
    Rmax = max(r) if gp.maxR < 0 else 1.0*gp.maxR           # [Rscale0]
    print('Rmax [Rscale] = ', Rmax)
    R = R[(R<Rmax)] # [Rscale0]

    Binmin, Binmax, Rbin = gpr.determine_radius(R, Rmin, Rmax, gp) # [Rscale0]
    vol = volume_spherical_shell(Binmin, Binmax, gp) # [Rscale0^3]

    Rscale0 = gfile.read_Rscale(gp.files.get_scale_file(0)+'_3D')

    for comp in range(gpr.ncomp):
        print('#######  working on component ',comp)
        print('input: ',gpr.get_com_file(comp)+'_3D')
        # start from data centered on COM already:
        if gfile.bufcount(gpr.get_com_file(comp)+'_3D')<2: continue
        x,y,z,v = np.loadtxt(gpr.get_com_file(comp)+'_3D',\
                           skiprows=1,usecols=(0,1,2,3),unpack=True) # 3*[Rscale], [km/s]

        # calculate 2D radius on the skyplane
        r = np.sqrt(x**2+y**2) # [Rscale]
        
        # set maximum radius (if gp.maxR is set)
        rmax = max(r) if gp.maxR<0 else 1.0*gp.maxR # [Rscale]
        print('rmax [Rscale] = ', rmax)
        sel = (r<=rmax)
        x = x[sel]; y = y[sel]; z = z[sel]; v = v[sel]; r = r[sel] # [Rscale]
        totmass = 1.*len(x) # [Munit], Munit = 1/star
            
        rs = r                   # + possible starting offset, [Rscale]
        vlos = v                 # + possible starting offset, [km/s]
        
        gfile.write_tracer_file(gp.files.get_ntracer_file(comp)+'_3D', totmass)
        de, em = gfile.write_headers_3D(gp, comp)

        # gpr.n=30 iterations for getting random picked radius values
        density = np.zeros((gp.nipol,gpr.n))
        a       = np.zeros((gp.nipol,gpr.n)) # shared by density, siglos, kappa calcs
        for k in range(gpr.n):
            rsi = gpr.Rerr * np.random.randn(len(rs)) + rs # [Rscale]
            vlosi = gpr.vrerr*np.random.randn(len(vlos)) + vlos # [km/s]
            for i in range(gp.nipol):
                ind1 = np.argwhere(np.logical_and(rsi>=Binmin[i], rsi<Binmax[i])).flatten() # [1]
                density[i][k] = (1.*len(ind1))/vol[i]*totmass # [Munit/Rscale^2]
                vlos1 = vlosi[ind1]                           # [km/s]
                a[i][k] = 1.*len(ind1)                        # [1]

                
        dens0 = np.sum(density[0])/(1.*gpr.n) # [Munit/Rscale^3]
        print('dens0 = ',dens0,' [Munit/Rscale^3]')

        dens0pc = dens0/Rscale0**3
        gfile.write_density(gp.files.get_scale_file(comp)+'_3D', dens0pc, totmass)
        
        tpb0   = np.sum(a[0])/float(gpr.n)     # [1] tracers per bin
        denserr0 = dens0/np.sqrt(tpb0)       # [Munit/Rscale^3]
        p_dens  = np.zeros(gp.nipol);  p_edens = np.zeros(gp.nipol)
        for b in range(gp.nipol):
            dens = np.sum(density[b])/float(gpr.n) # [Munit/Rscale^3]
            tpb  = np.sum(a[b])/float(gpr.n)       # [1]
            denserr = dens/np.sqrt(tpb)            # [Munit/Rscale^3]

            if(np.isnan(denserr)):
                p_dens[b] = p_dens[b-1]                          # [1]
                p_edens[b]= p_edens[b-1]                         # [1]
            else:
                p_dens[b] = dens/dens0                           # [1]
                p_edens[b]= denserr/dens0 # [1] #100/rbin would be artificial guess

            print(Rbin[b], Binmin[b], Binmax[b], p_dens[b], p_edens[b], file=de) # [Rscale], 2*[dens0]
            indr = (r<Binmax[b])
            menclosed = float(np.sum(indr))/totmass # for normalization to 1  # [totmass]
            merr = menclosed/np.sqrt(tpb) # artificial menclosed/10 # [totmass]
            print(Rbin[b], Binmin[b], Binmax[b], menclosed, merr, file=em) # [rscale], 2*[totmass]
        de.close()
        em.close()

        if gpr.showplots:
            show_plots_dens(Rbin, p_dens, p_edens, gp)
            
if __name__ == '__main__':
    gpr.showplots = True
    import gl_params
    gp = gl_params.Params()
    run(gp)

