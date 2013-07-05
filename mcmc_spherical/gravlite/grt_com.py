#!/usr/bin/python
# (c) 2013 Pascal S.P. Steger
'''calculate approximative center of mass, assuming constant stellar mass, for triax. system'''

import numpy as np
import sys
import pdb

from pylab import *
import gl_params as gp
import gr_params as gpr
from gl_helper import expDtofloat
from gl_class_files import *


def com_mean(x,y):
    '''mean COM, weighted by probability of membership'''
    pm = np.ones(len(x))
    com_x = 1.*np.sum(x*pm)/np.sum(pm) # [pc]
    com_y = 1.*np.sum(y*pm)/np.sum(pm) # [pc]
    return com_x, com_y


def com_shrinkcircle(x,y):
    print 'shrinking sphere'
    eps = 1e-6
    pm = np.ones(len(x))
    com_x = 1.*np.sum(x*pm)/np.sum(pm);    com_y = 1.*np.sum(y*pm)/np.sum(pm)
    bucom_x = com_x; bucom_y = com_y
    x -= com_x; y -= com_y
    dr = np.sqrt(com_x**2+com_y**2)
    r0 = np.sqrt(x**2+y**2)

    nit = 0; minlen = len(x)/2.
    while nit < 200 and len(x) > minlen:
        nit += 1
        print 'iteration ',nit,' with ',len(x), ' particles has overall COM of: ',bucom_x,bucom_y,' with remaining offset ',dr

        # shrink sphere:
        # 1) calc radius
        r0 = np.sqrt(x**2+y**2)
        # 2) sort remaining particles
        order = np.argsort(r0)
        r0 = np.array(r0)[order]; x = np.array(x)[order]; y = np.array(y)[order]; pm = np.array(pm)[order]

        # 3) cut x,y,z,pm after 1-10%
        end = len(r0)*0.95
        r0 = r0[:end]; x = x[:end]; y = y[:end]; pm = pm[:end]

        # calculate new COM
        com_x = 1.*np.sum(x*pm)/np.sum(pm);    com_y = 1.*np.sum(y*pm)/np.sum(pm)
        dr = np.sqrt(com_x**2+com_y**2)

        # add to bucom
        bucom_x -= com_x; bucom_y -= com_y

        # recenter particles
        x -= com_x; y -= com_y

    return bucom_x, bucom_y


def run():
    print 'input:'
    print gpr.fil
    x0,y0,vlos = np.genfromtxt(gpr.fil, skiprows=0, unpack =  True,
                               usecols = (0,1,5))
    x0 *= 1000.                         # [pc]
    y0 *= 1000.                         # [pc]
    
    # center of mass with means
    #com_x, com_y = com_mean(x0,y0) # [TODO]
    
    # shrinking sphere method
    com_x, com_y = com_shrinkcircle(x0,y0)
    print 'COM [pc]: ', com_x, com_y


    com_vz = np.sum(vlos)/(1.*len(vlos)) # [km/s]
    print 'VOM [km/s]', com_vz


    x0 -= com_x; y0 -= com_y # [pc]
    vlos -= com_vz #[km/s]
    
    rc = np.sqrt(x0**2+y0**2) # [pc]
    rc.sort() # [pc]
    for i in range(len(rc)-1):
        if rc[i]>rc[i+1]: #[pc]
            print 'sorting error!'
            exit(1)
    rhalf = rc[floor(len(rc)/2)] # [pc]
    rcore = rhalf # or gpr.r_DM # [pc]
    print 'rcore = ',rcore,' pc'
    print 'max(r) = ',max(rc),' pc'
    print 'last element of r : ',rc[-1],' pc'
    print 'total number of stars: ',len(rc)

    r0 = np.sqrt(x0**2+y0**2)/rcore
    sel = (r0<gpr.rprior)
    x = x0[sel]/rcore; y = y0[sel]/rcore # [r_core]
    vz=vlos[sel]
    m = np.ones(len(x))
    r = np.sqrt(x*x+y*y) #[r_core]
    # print "x y z" on first line, to interprete data later on
    crcore = open(gpr.get_params_file(0),'w')
    print >> crcore, '# rcore in [pc], surfdens_central (=dens0) in [munit/rcore**2], and in [munit/pc**2], and totmass [munit], and max(v_LOS) in [km/s]'
    print >> crcore, rcore
    crcore.close()

    print 'output: ',gpr.get_com_file(0)
    c = open(gpr.get_com_file(0),'w')
    print >> c,'# x [rcore],','y [rcore],','vLOS [km/s],','rcore = ',rcore,' pc'
    for k in range(len(x)):
        print >> c,x[k],y[k],vz[k] #[rcore], [rcore], [km/s]
    c.close()
        
        
    if not gp.testplot_read: return
        
    ion(); ax = subplot(111)
    # res = (abs(x)<3)*(abs(y)<3)
    # x = x[res]; y = y[res] #[rcore]
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
    circ_HL=Circle((0,0), radius=rcore/rcore, fc='None', ec='g', lw=3)
    circ_DM=Circle((0,0), radius=gpr.r_DM/rcore, fc='None', ec='r', lw=3)
    gca().add_patch(circ_HL); gca().add_patch(circ_DM)
    
    # visible region
    axes().set_aspect('equal')
    
    xlabel(r'$x [R_s]$'); ylabel(r'$y [R_s]$')
    # title(gpr.fil)
    savefig(gpr.get_com_png(0))
    if gpr.showplots:
        ioff();show();clf()
    
if __name__=='__main__':
    # gp.testplot_read = True
    run()
