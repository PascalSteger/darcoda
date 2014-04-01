#!/usr/bin/env python3

##
# @file
# calculate approximative center of mass, assuming constant stellar mass
# 2D version, see grw_com for 3D

# (c) 2013 Pascal S.P. Steger

import numpy as np
import sys, pdb

from pylab import *
import gr_params as gpr
from gl_helper import expDtofloat
from gl_class_files import *
from gl_centering import *

def run(gp):
    print('input: ', gpr.fil)
    x0,y0,z0,vz0,vb0,Mg0,PM0,comp0 = np.genfromtxt(gpr.fil, skiprows = 0, unpack = True,\
                                                   usecols=(0, 1, 2, 5, 12, 13, 19, 20),\
                                                   dtype="d17",\
                                                   converters={0:expDtofloat,  # x0  in pc \
                                                               1:expDtofloat,  # y0  in pc \
                                                               2:expDtofloat,  # z0  in pc \
                                                               5:expDtofloat, # vz0 in km/s\
                                                               12:expDtofloat, # vb0(LOS due binary), km/s\
                                                               13:expDtofloat, # Mg0 in Angstrom\
                                                               19:expDtofloat, # PM0 [1]\
                                                               20:expDtofloat}) # comp0 1,2,3(background)
    # use component 12-1 instead of 6-1 for z velocity, to exclude observational errors

    # only use stars which are members of the dwarf: exclude pop3 by construction
    pm = (PM0 >= gpr.pmsplit) # exclude foreground contamination, outliers
    PM0 = PM0[pm]
    comp0 = comp0[pm]; x0=x0[pm]; y0=y0[pm]; z0=z0[pm]; vz0=vz0[pm]; vb0=vb0[pm]; Mg0=Mg0[pm]

    # assign population
    if gp.pops==2:
        pm1 = (comp0 == 1) # will be overwritten below if gp.metalpop is set
        pm2 = (comp0 == 2) # same same
    elif gp.pops==1:
        pm1 = (comp0 < 3)
        pm2 = (comp0 == -1) # assign none, but of same length as comp0

    if gp.metalpop:
        # drawing of populations based on metallicity
        # get parameters from function in pymcmetal.py
        import pymcmetal as pmc
        p, mu1, sig1, mu2, sig2, M = pmc.bimodal_gauss(Mg0)
        pm1, pm2 = pmc.assign_pop(Mg0, p, mu1, sig1, mu2, sig2)
        # output: changed pm1, pm2

    x1=x0[pm1]; y1=y0[pm1]; z1=z0[pm1]; comp1=comp0[pm1]; vz1=vz0[pm1]; vb1=vb0[pm1]; Mg1=Mg0[pm1]; PM1=PM0[pm1]
    x2=x0[pm2]; y2=y0[pm2]; z2=z0[pm2]; comp2=comp0[pm2]; vz2=vz0[pm2]; vb2=vb0[pm2]; Mg2=Mg0[pm2]; PM2=PM0[pm2]
    
    # cutting pm_i to a maximum of ntracers_i particles each:
    from random import shuffle
    
    ind1 = np.arange(len(x1))
    np.random.shuffle(ind1)     # random.shuffle already changes ind
    ind1 = ind1[:gp.ntracer[1-1]]

    ind2 = np.arange(len(x2))
    np.random.shuffle(ind2)     # random.shuffle already changes ind
    ind2 = ind2[:gp.ntracer[2-1]]

    x1=x1[ind1];y1=y1[ind1];z1=z1[ind1];comp1=comp1[ind1];vz1=vz1[ind1];vb1=vb1[ind1];Mg1=Mg1[ind1]; PMS1=PM1[ind1]
    x2=x2[ind2];y2=y2[ind2];z2=z2[ind2];comp2=comp2[ind2];vz2=vz2[ind2];vb2=vb2[ind2];Mg2=Mg2[ind2]; PMS2=PM2[ind2]

    x0 = np.hstack([x1, x2])
    y0 = np.hstack([y1, y2])
    z0 = np.hstack([z1, z2])
    vz0 = np.hstack([vz1, vz2])
    N1  = min(len(x1), gp.ntracer[1-1])
    N2  = min(len(x2), gp.ntracer[2-1])
    pm1 = np.hstack([np.ones(N1, dtype=bool), np.zeros(N2, dtype=bool)])
    pm2 = np.hstack([np.zeros(N1, dtype=bool), np.ones(N2,  dtype=bool)])
    pm = pm1 + pm2
    
    # get COM with shrinking sphere method
    com_x, com_y, com_z, com_vz = com_shrinkcircle_v(x0, y0, z0, vz0, pm) # [pc]
    print('COM [pc]: ', com_x, com_y, com_z)   # [pc]
    print('VOM [km/s]', com_vz)                # [pc]

    # from now on, work with 2D data only; z0 was only used to get center in (x,y) better
    x0 -= com_x; y0 -= com_y # [pc]
    vz0 -= com_vz #[km/s]
    
    R0 = np.sqrt(x0**2+y0**2) # [pc]
    Rc = R0                   # [pc]
    Rc.sort()                 # [pc]
    Rhalf = Rc[floor(len(Rc)/2)]        # [pc]
    Rscale = Rhalf                      # or gpr.r_DM # [pc]
    # Rscale = gpr.r_DM                 # we only work with data here, so gpr.r_DM is not known
    print('Rscale = ', Rscale, ' pc')
    print('max(R) = ', max(Rc), ' pc')
    print('last element of R : ', Rc[-1], ' pc')
    print('total number of stars: ', len(Rc))
    
    x0 = x0/Rscale; y0 = y0/Rscale           # [Rscale]
    
    i = -1
    for pmn in [pm, pm1, pm2]:
        pmr = ( R0 < (gp.maxR*Rscale) )  # read max extension for data (rprior*Rscale) from gl_params
        pmn = pmn*pmr                   # [1]
        print("fraction of members = ", 1.0*sum(pmn)/len(pmn))
        i = i+1
        x=x0[pmn]; y=y0[pmn]; vz=vz0[pmn]; vb=vb0[pmn];  # [1], [km/s]
        Mg=Mg0[pmn]; comp=comp0[pmn]; PMN=PM0[pmn]   # [ang], [1], [1]
        m = np.ones(len(pmn))
        
        R = np.sqrt(x*x+y*y)            # [Rscale]

        # store central values in scale_ file
        crscale = open(gp.files.get_scale_file(i),'w')
        print('# Rscale in [pc], surfdens_central (=dens0) in [munit/rscale**2], and in [munit/pc**2], and totmass [munit], and max(v_LOS) in [km/s]', file=crscale)
        print(np.median(R)*Rscale, file=crscale)
        crscale.close()

        print('output: ',gpr.get_com_file(i))
        c = open(gpr.get_com_file(i),'w')
        print('# x [Rscale],','y [Rscale],','vLOS [km/s],','Rscale = ',Rscale,' pc', file=c)
        for k in range(len(x)):
            print(x[k],y[k],vz[k], file=c)      # [Rscale], [Rscale], [km/s]
        c.close()
        
        
        if not gpr.showplots: continue
        
        ion(); subplot(111)
        res = (abs(x)<3)*(abs(y)<3)
        x = x[res]; y = y[res]           # [Rscale]
        en = len(x)
        if en == 0: continue
        scatter(x[:en], y[:en], c=pmn[:en], s=35, vmin=0.95, vmax=1.0, lw=0.0, alpha=0.2)
        # xscale('log'); yscale('log')
        if i == 0: colorbar()
        circ_HL=Circle((0,0), radius=Rscale/Rscale, fc='None', ec='g', lw=1)
        circ_DM=Circle((0,0), radius=gpr.r_DM/Rscale, fc='None', ec='r', lw=1)
        gca().add_patch(circ_HL); gca().add_patch(circ_DM)
        
        # visible region
        maxr = max(np.abs(x));  mayr = max(np.abs(y)) #[rscale]
        width2 = max([maxr,mayr]) #[rscale]
        xlim([-width2,width2]); ylim([-width2,width2])
        axes().set_aspect('equal')
    
        xlabel(r'$x [R_s]$'); ylabel(r'$y [R_s]$')
        # title(gpr.fil)
        savefig(gpr.get_com_png(i))
        if gpr.showplots:
            ioff();show();clf()
    
if __name__=='__main__':
    gpr.showplots = True
    import gl_params
    gp = gl_params.Params()
    run(gp)
