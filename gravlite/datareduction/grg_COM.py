#!/usr/bin/env python3

##
# @file
# calculate approximative center of mass, assuming constant stellar mass
# 3D version, see grw_COM for 2D

# (c) 2013 ETHZ Pascal S.P. Steger, psteger@phys.ethz.ch

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
    print('input:',gpr.fil)
    x0,y0,z0,vx,vy,vz=np.transpose(np.loadtxt(gpr.fil))
    
    # cutting pm_i to a maximum of ntracers particles:
    from random import shuffle
    ind = np.arange(len(x0))
    np.random.shuffle(ind)
    ind = ind[:gp.files.ntracer]
    x0=x0[ind]; y0=y0[ind]; z0=z0[ind];
    vx = vx[ind]; vy = vy[ind]; vz=vz[ind]
    
    # get center of mass with means
    #com_x, com_y,com_z = com_mean(x0,y0,z0,PM0) 
    # [TODO], and TODO: z component included if available
    
    # get COM with shrinking sphere method
    PM = np.ones(len(x0)) # assign all particles the full probability of membership
    com_x, com_y, com_z, com_vz = com_shrinkcircle_v(x0,y0,z0,vz,PM)
    print('COM [pc]: ', com_x, com_y, com_z)
    print('VOM [km/s]', com_vz)

    # from now on, work with 2D data only; 
    # z0 was only used to get center in (x,y) better
    
    x0 -= com_x; y0 -= com_y # [pc]
    vz -= com_vz #[km/s]
    
    R0 = np.sqrt(x0**2+y0**2) # [pc]
    Rc = R0                   # [pc]
    Rc.sort()                 # [pc]
    for i in range(len(Rc)-1):
        if Rc[i]>Rc[i+1]:               # [pc]
            print('sorting error!')
            exit(1)
    Rhalf = Rc[floor(len(Rc)/2)]        # [pc]
    Rcore = Rhalf                       # or gpr.r_DM # [pc]
    #Rcore = gpr.r_DM                   # TODO: delete, only work with data
    print('Rcore = ',Rcore,' pc')
    print('max(R) = ',max(Rc),' pc')
    print('last element of R : ',Rc[-1],' pc')
    print('total number of stars: ',len(Rc))
    
    x0 = x0/Rcore; y0 = y0/Rcore           # [Rcore]
    
    i = -1
    for comp in range(2):
        pmr = (R0<(gpr.rprior*Rcore))  # TODO: read from gl_class_file
        i = i+1
        m = np.ones(len(R0))
        x = x0; y = y0
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
        figure(1)
        subplot(111)
        res = (abs(x)<3)*(abs(y)<3)
        x = x[res]; y = y[res]                     # [Rcore]
        en = len(x)
        if en == 0: continue
        scatter(x[:en], y[:en],\
                s=35, vmin=0.95, vmax=1.0, lw=0.0, alpha=0.2)
        # xscale('log'); yscale('log')
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
        show(block=True)
    
if __name__=='__main__':
    # gpr.showplots = True
    run()
