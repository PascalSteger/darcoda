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
import gl_params as gp
import gr_params as gpr
from gl_helper import expDtofloat
from gl_class_files import *
from gl_centering import *

def run():
    print('input:')
    print(gpr.fil)
    x0,y0,z0,vz0,vb0,Mg0,PM0,comp0=np.genfromtxt(gpr.fil,skiprows=0,unpack=True,\
                                                 usecols=(0,1,2,11,12,13,19,20),\
                                                 dtype="d17",\
                                                 converters={0:expDtofloat,  # x0  in pc \
                                                             1:expDtofloat,  # y0  in pc \
                                                             2:expDtofloat,  # z0  in pc \
                                                             11:expDtofloat, # vz0 in km/s\
                                                             12:expDtofloat, # vb0(LOS due binary), km/s\
                                                             13:expDtofloat, # Mg0 in Angstrom\
                                                             19:expDtofloat, # PM0 [1]\
                                                             20:expDtofloat}) # comp0 1,2,3(background)
    
    
    # TODO: use component 6-1 instead of 12-1 for z velocity, to include observational errors

    # only use stars which are members of the dwarf: exclude pop3 by construction
    pm = (PM0 >= gpr.pmsplit) # exclude foreground contamination, outliers
    PM0 = PM0[pm]
    comp0 = comp0[pm]; x0=x0[pm]; y0=y0[pm]; z0=z0[pm]; vz0=vz0[pm]; vb0=vb0[pm]; Mg0=Mg0[pm]
    pm1 = (comp0 == 1) # will be overwritten below if gp.metalpop
    pm2 = (comp0 == 2) # same same
    
    
    if gp.metalpop:
        # drawing of populations based on metallicity
        # get parameters from function in pymcmetal.py
        import pymcmetal2 as pmc
        p,mu1,sig1,mu2,sig2, M = pmc.bimodal_gauss(Mg0)
        pm1, pm2 = pmc.assign_pop(Mg0,p,mu1,sig1,mu2,sig2)
        # output: changed pm1, pm2

    # cutting pm_i to a maximum of ntracers particles:
    from random import shuffle
    ind = np.arange(len(x0))
    np.random.shuffle(ind)
    ind = ind[:gp.files.ntracer]
    x0=x0[ind]; y0=y0[ind]; z0=z0[ind]; comp0=comp0[ind]; vz0=vz0[ind]; vb0=vb0[ind]; Mg0=Mg0[ind]
    PM0 = PM0[ind]; pm1 = pm1[ind]; pm2 = pm2[ind]; pm = pm1+pm2
    
    # get center of mass with means
    #com_x, com_y,com_z = com_mean(x0,y0,z0,PM0) # [TODO], and TODO: z component included if available
    
    # get COM with shrinking sphere method
    com_x, com_y, com_z, com_vz = com_shrinkcircle(x0,y0,z0,vz0,PM0)
    print('COM [pc]: ', com_x, com_y, com_z)
    print('VOM [km/s]', com_vz)

    # from now on, work with 2D data only; z0 was only used to get center in (x,y) better
    
    x0 -= com_x; y0 -= com_y # [pc]
    vz0 -= com_vz #[km/s]
    
    R0 = np.sqrt(x0**2+y0**2) # [pc]
    Rc = R0                   # [pc]
    Rc.sort()                 # [pc]
    for i in range(len(Rc)-1):
        if Rc[i]>Rc[i+1]:               # [pc]
            print('sorting error!')
            exit(1)
    Rhalf = Rc[floor(len(Rc)/2)]        # [pc]
    # Rcore = Rhalf                       # or gpr.r_DM # [pc]
    Rcore = gpr.r_DM                    # TODO: delete, only work with data
    print('Rcore = ',Rcore,' pc')
    print('max(R) = ',max(Rc),' pc')
    print('last element of R : ',Rc[-1],' pc')
    print('total number of stars: ',len(Rc))
    
    x0 = x0/Rcore; y0 = y0/Rcore           # [Rcore]
    
    i = -1
    for pmn in [pm,pm1,pm2]:
        pmr = (R0<(gpr.rprior*Rcore))  # TODO: read from gl_class_file
        pmn = pmn*pmr                  # [1]
        print("fraction of members = ",1.0*sum(pmn)/len(pmn))
        i = i+1
        x=x0[pmn]; y=y0[pmn]; vz=vz0[pmn]; vb=vb0[pmn];  # [1], [km/s]
        Mg=Mg0[pmn]; comp=comp0[pmn]; PMN=PM0[pmn]   # [ang], [1], [1]
        m = np.ones(len(pmn))
        
        R = np.sqrt(x*x+y*y)            # [Rcore]
        
        # print("x y z" on first line, to interprete data later on)
        crcore = open(gpr.get_params_file(i),'w')
        print('# Rcore in [pc], surfdens_central (=dens0) in [munit/rcore**2], and in [munit/pc**2], and totmass [munit], and max(v_LOS) in [km/s]', file=crcore)
        print(Rcore, file=crcore)
        crcore.close()

        print('output: ',gpr.get_com_file(i))
        c = open(gpr.get_com_file(i),'w')
        print('# x [Rcore],','y [Rcore],','vLOS [km/s],','Rcore = ',Rcore,' pc', file=c)
        for k in range(len(x)):
            print(x[k],y[k],vz[k], file=c)      # [Rcore], [Rcore], [km/s]
        c.close()
        
        
        if not gpr.showplots: continue
        
        ion(); subplot(111)
        res = (abs(x)<3)*(abs(y)<3)
        x = x[res]; y = y[res]           # [Rcore]
        en = len(x)
        if en == 0: continue
        scatter(x[:en], y[:en], c=pmn[:en], s=35, vmin=0.95, vmax=1.0, lw=0.0, alpha=0.2)
        # xscale('log'); yscale('log')
        if i == 0: colorbar()
        circ_HL=Circle((0,0), radius=Rcore/Rcore, fc='None', ec='g', lw=1)
        circ_DM=Circle((0,0), radius=gpr.r_DM/Rcore, fc='None', ec='r', lw=1)
        gca().add_patch(circ_HL); gca().add_patch(circ_DM)
        
        # visible region
        maxr = max(np.abs(x));  mayr = max(np.abs(y)) #[rcore]
        width2 = max([maxr,mayr]) #[rcore]
        xlim([-width2,width2]); ylim([-width2,width2])
        axes().set_aspect('equal')
    
        xlabel(r'$x [R_s]$'); ylabel(r'$y [R_s]$')
        # title(gpr.fil)
        savefig(gpr.get_com_png(i))
        if gpr.showplots:
            ioff();show();clf()
    
if __name__=='__main__':
    # gpr.showplots = True
    run()
