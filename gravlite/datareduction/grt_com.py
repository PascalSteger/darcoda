#!/usr/bin/env python3

##
# @file
# calculate approximative center of mass, assuming constant stellar mass, for triax. system

# (c) 2013 Pascal S.P. Steger, psteger@phys.ethz.ch


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
    x0,y0,vlos = np.genfromtxt(gpr.fil, skiprows=0, unpack =  True,
                               usecols = (0,1,5))

    # use only 3000 random particles:
    ind = np.arange(len(x0))
    np.random.shuffle(ind)
    ind = ind[:3000]
    x0 = x0[ind];    y0 = y0[ind];    vlos = vlos[ind]
    
    x0 *= 1000.                         # [pc]
    y0 *= 1000.                         # [pc]
    
    # center of mass with means
    #com_x, com_y = com_mean(x0,y0) # [TODO]
    
    # shrinking sphere method
    com_x, com_y, com_vz = com_shrinkcircle_v_2D(x0,y0, vlos)
    print('COM [pc]: ', com_x, com_y)
    print('VOM [km/s]', com_vz)

    x0 -= com_x; y0 -= com_y # [pc]
    vlos -= com_vz #[km/s]
    
    rc = np.sqrt(x0**2+y0**2) # [pc]
    rc.sort() # [pc]
    for i in range(len(rc)-1):
        if rc[i]>rc[i+1]: #[pc]
            print('sorting error!')
            exit(1)
    rhalf = rc[floor(len(rc)/2)] # [pc]
    rscale = rhalf # or gpr.r_DM # [pc]
    print('rscale = ',rscale,' pc')
    print('max(r) = ',max(rc),' pc')
    print('last element of r : ',rc[-1],' pc')
    print('total number of stars: ',len(rc))

    r0 = np.sqrt(x0**2+y0**2)/rscale
    sel = (r0<gpr.rprior)
    x = x0[sel]/rscale; y = y0[sel]/rscale # [r_scale]
    vz = vlos[sel]
    m = np.ones(len(x))
    r = np.sqrt(x*x+y*y) #[r_scale]
    # print("x y z") on first line, to interprete data later on
    crscale = open(gpr.get_params_file(0),'w')
    print('# rscale in [pc], surfdens_central (=dens0) in [munit/rscale**2], and in [munit/pc**2], and totmass [munit], and max(v_LOS) in [km/s]', file=crscale)
    print(rscale, file=crscale)
    crscale.close()

    print('output: ',gpr.get_com_file(0))
    c = open(gpr.get_com_file(0),'w')
    print('# x [rscale],','y [rscale],','vLOS [km/s],','rscale = ',rscale,' pc', file=c)
    for k in range(len(x)):
        print(x[k],y[k],vz[k], file=c) #[rscale], [rscale], [km/s]
    c.close()
        
        
    if not gpr.showplots: return
        
    ion(); ax = subplot(111)
    # res = (abs(x)<3)*(abs(y)<3)
    # x = x[res]; y = y[res] #[rscale]
    en = len(x)

    H, xedges, yedges = np.histogram2d(x, y,  bins=(30,30),  range=[[-3.,3.], [-3.,3.]])
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    subplots_adjust(bottom=0.15, left=0.15)


    imshow(H, interpolation='bilinear', origin='lower',
           cmap=cm.gray, extent=(-3,3,-3,3))
    
    levels = np.logspace(2,4,10)
    cset = contour(H, levels, origin='lower', extent=extent) #cmap=cm.gray
    clabel(cset, inline=1, fontsize=8, fmt='%1.0i')
    for c in cset.collections:
        c.set_linestyle('solid')
        
    # scatter(x[:en], y[:en], s=35, vmin=0.95, vmax=1.0, lw=0.0, alpha=0.2)
    # xscale('log'); yscale('log')
    circ_HL=Circle((0,0), radius=rscale/rscale, fc='None', ec='g', lw=3)
    circ_DM=Circle((0,0), radius=gpr.r_DM/rscale, fc='None', ec='r', lw=3)
    gca().add_patch(circ_HL); gca().add_patch(circ_DM)
    
    # visible region
    axes().set_aspect('equal')
    
    xlabel(r'$x [R_s]$'); ylabel(r'$y [R_s]$')
    # title(gpr.fil)
    savefig(gpr.get_com_png(0))
    if gpr.showplots:
        ioff();show();clf()
    
if __name__=='__main__':
    # gpr.showplots = True
    import gl_params
    gp = gl_params.Params()
    run(gp)
