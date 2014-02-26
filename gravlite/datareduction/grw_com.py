#!/usr/bin/env python3

##
# @file
# calculate approximative center of mass, assuming constant stellar mass
# 3D version, see grw_COM for 2D

# (c) 2013 Pascal S.P. Steger

import numpy as np
import sys
import pdb

from pylab import *
import gr_params as gpr
from gl_helper import expDtofloat
from gl_class_files import *
from gl_centering import *

def run(gp):
    print('input: ', gpr.fil)
    x0,y0,z0,vz0,vb0,Mg0,PM0,comp0=np.genfromtxt(gpr.fil,skiprows=0,unpack=True,\
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
    pm1 = (comp0 == 1) # will be overwritten below if gp.metalpop
    pm2 = (comp0 == 2) # same same
    pm3 = (comp0 == 3)

    
    if gp.metalpop:
        # drawing of populations based on metallicity
        # get parameters from function in pymcmetal.py
        import pymcmetal as pmc
        p, mu1, sig1, mu2, sig2, M = pmc.bimodal_gauss(Mg0)   # [TODO]
        pm1, pm2 = pmc.assign_pop(Mg0, p, mu1, sig1, mu2, sig2)   # [1]
        # output: changed pm1, pm2
        # assumption: no component 3 stars are included

    # cutting pm_i to a maximum of ntracers particles:
    from random import shuffle
    ind = np.arange(len(x0))
    np.random.shuffle(ind)
    ind = ind[:gp.files.ntracer]
    x0=x0[ind]; y0=y0[ind]; z0=z0[ind]; comp0=comp0[ind]; vz0=vz0[ind]; vb0=vb0[ind]; Mg0=Mg0[ind]
    PM0 = PM0[ind]; pm1 = pm1[ind]; pm2 = pm2[ind]; pm3 = pm3[ind]; pm = pm1+pm2+pm3
    
    # get center of mass with means
    #com_x, com_y,com_z = com_mean(x0,y0,z0,PM0) # [TODO], and TODO: z component included if available
    
    # get COM with shrinking sphere method
    com_x, com_y, com_z = com_shrinkcircle(x0,y0,z0,PM0)
    print('COM [pc]: ', com_x, com_y, com_z)


    com_vz = np.sum(vz0*PM0)/np.sum(PM0) # [km/s]
    print('VOM [km/s]', com_vz)

    # from now on, continue to work with 3D data. store to different files
    
    x0 -= com_x; y0 -= com_y; z0 -= com_z # [pc]
    vz0 -= com_vz #[km/s]

    # but still get the same radii as from 2D method, to get comparison of integration routines right
    R0 = np.sqrt(x0**2+y0**2) # [pc]
    Rc = R0                   # [pc]
    Rc.sort()                 # [pc]
    for i in range(len(Rc)-1):
        if Rc[i]>Rc[i+1]:               # [pc]
            print('sorting error!')
            exit(1)
    Rhalf = Rc[floor(len(Rc)/2)]        # [pc]
    Rscale = Rhalf                       # or gpr.r_DM # [pc]
    print('Rscale = ', Rscale,  ' pc')
    print('max(R) = ', max(Rc) ,' pc')
    print('last element of R : ',Rc[-1],' pc')
    print('total number of stars: ',len(Rc))
    
    x0 = x0/Rscale; y0 = y0/Rscale; z0 = z0/Rscale              # [Rscale]
    
    i = -1
    for pmn in [pm,pm1,pm2,pm3]:
        pmr = (R0<(gpr.rprior*Rscale))  # TODO: read from gl_class_file
        pmn = pmn*pmr                  # [1]
        print("fraction of members = ",1.0*sum(pmn)/len(pmn))
        i = i+1
        x=x0[pmn]; y=y0[pmn]; z=z0[pmn]; vz=vz0[pmn]; vb=vb0[pmn];  # [1], [km/s]
        Mg=Mg0[pmn]; comp=comp0[pmn]; PMN=PM0[pmn]   # [ang], [1], [1]
        m = np.ones(len(pmn))
        
        R = np.sqrt(x*x+y*y)            # [Rscale]
        
        # print("x y z" on first line, to interprete data later on)
        crscale = open(gpr.get_params_file(i)+'_3D','w')
        print('# Rscale in [pc], surfdens_central (=dens0) in [munit/rscale**2], and in [munit/pc**2], and totmass [munit], and max(v_LOS) in [km/s]', file=crscale)
        print(Rscale, file=crscale)
        crscale.close()

        print('output: ',gpr.get_com_file(i)+'_3D')
        c = open(gpr.get_com_file(i)+'_3D','w')
        print('# x [Rscale],','y [Rscale],', 'z [Rscale]','vLOS [km/s],','Rscale = ',Rscale,' pc', file=c)
        for k in range(len(x)):
            print(x[k],y[k],z[k],vz[k], file=c)      # 3* [Rscale], [km/s]
        c.close()
        
        
        if not gpr.showplots: continue
        # plot x-z values
        ion(); subplot(111)
        res = (abs(x)<3)*(abs(z)<3)
        x = x[res]; z = z[res]           # [Rscale]
        en = len(x)
        if en == 0: continue
        scatter(x[:en], z[:en], c=pmn[:en], s=35, vmin=0.95, vmax=1.0, lw=0.0, alpha=0.2)
        # xscale('log'); yscale('log')
        if i == 0: colorbar()
        circ_HL=Circle((0,0), radius=Rscale/Rscale, fc='None', ec='g', lw=1)
        circ_DM=Circle((0,0), radius=gpr.r_DM/Rscale, fc='None', ec='r', lw=1)
        gca().add_patch(circ_HL); gca().add_patch(circ_DM)
        
        # visible region
        maxr = max(np.abs(x));  mayr = max(np.abs(z)) #[rscale]
        width2 = max([maxr,mayr]) #[rscale]
        xlim([-width2,width2]); ylim([-width2,width2])
        axes().set_aspect('equal')
    
        xlabel(r'$x [R_s]$'); ylabel(r'$z [R_s]$')
        # title(gpr.fil)
        savefig(gpr.get_com_png(i)+'_3D.png')
        if gpr.showplots:
            ioff();show();clf()
    
if __name__=='__main__':
    # gpr.showplots = True
    run(gp)
