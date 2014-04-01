#!/usr/bin/env python3

##
# @file
# calculate surface mass density falloff of circular rings around center of mass
# for triaxial systems
# TODO: enable bins in elliptical form

# (c) 2013 Pascal S.P. Steger


import sys
import pdb
import numpy as np
from pylab import *
import math
import gr_params as gpr
import gl_file as gfile
from gl_helper import expDtofloat, bin_r_linear, bin_r_log, bin_r_const_tracers
from gl_class_files import *

def run(gp):
    print('input: ',gpr.get_com_file(0))
    # start from data centered on COM already:
    x,y,v = np.loadtxt(gpr.get_com_file(0),\
                       skiprows=1,usecols=(0,1,2),unpack=True) #[rscale], [rscale], [km/s]
    
    # calculate 2D radius on the skyplane
    r = np.sqrt(x**2+y**2) #[rscale]
    
    # set number and size of (linearly spaced) bins
    rmin = 0. #[rscale]
    rmax = max(r) if gp.maxR < 0 else 1.0*gp.maxR #[rscale]
        
    print('rmax [rscale] = ', rmax)
    sel = (r<rmax)
    x = x[sel]; y = y[sel]; v = v[sel] #[rscale]
    totmass = 1.*len(x) #[munit], munit = 1/star
    
    if gp.lograd:
        # space logarithmically in radius
        binmin, binmax, rbin = bin_r_log(rmax/gp.nipol, rmax, gp.nipol)
    elif gp.consttr:
        binmin, binmax, rbin = bin_r_const_tracers(r, len(r)/gp.nipol)
    else:
        binmin, binmax, rbin = bin_r_linear(rmin, rmax, gp.nipol)
            
    #volume of a circular bin from binmin to binmax
    vol = np.zeros(gp.nipol)
    for k in range(gp.nipol):
        vol[k] = np.pi*(binmax[k]**2-binmin[k]**2) # [rscale**2]
            
    # rs = gpr.Rerror*np.random.randn(len(r))+r
    rs = r  #[rscale] # if no initial offset is whished
    
    tr = open(gp.files.get_ntracer_file(0),'w')
    print(totmass, file=tr)
    tr.close()

    de = open(gp.files.nufiles[0],'w')
    em = open(gp.files.massfiles[0],'w')
    
    print('r','nu(r)/nu(0)','error', file=de)
    print('r','M(<r)','error', file=em)

    # 30 iterations for getting random picked radius values
    density = np.zeros((gp.nipol,gpr.n))
    a       = np.zeros((gp.nipol,gpr.n))
    for k in range(gpr.n):
        rsi = gpr.Rerror * np.random.randn(len(rs)) + rs # [rscale]
        for j in range(gp.nipol):
            ind1 = np.argwhere(np.logical_and(rsi>=binmin[j],rsi<binmax[j])).flatten() # [1]
            density[j][k] = (1.*len(ind1))/vol[j]*totmass # [munit/rscale**2]
            a[j][k] = 1.*len(ind1) #[1]
            
    dens0 = np.sum(density[0])/(1.*gpr.n) # [munit/rscale**2]
    print('dens0 = ',dens0,'[munit/rscale**2]')
    crscale = open(gp.files.get_scale_file(0),'r')
    rscale = np.loadtxt(crscale, comments='#', skiprows=1, unpack=False)
    crscale.close()

    cdens = open(gp.files.get_scale_file(0),'a')
    print(dens0, file=cdens)               # [munit/rscale**2]
    print(dens0/rscale**2, file=cdens)      # [munit/pc**2]
    print(totmass, file=cdens)             # [munit]
    cdens.close()
    
    ab0   = np.sum(a[0])/(1.*gpr.n)     # [1]
    denserr0 = dens0/np.sqrt(ab0)       # [munit/rscale**2]

    p_dens  = np.zeros(gp.nipol);  p_edens = np.zeros(gp.nipol)
    
    for b in range(gp.nipol):
        dens = np.sum(density[b])/(1.*gpr.n) # [munit/rscale**2]
        ab   = np.sum(a[b])/(1.*gpr.n)       # [1]
        denserr = dens/np.sqrt(ab)       # [munit/rscale**2]
        denserror = np.sqrt((denserr/dens0)**2+(dens*denserr0/(dens0**2))**2) #[1]
        if(math.isnan(denserror)):
            denserror = 0. # [1]
            p_dens[b] = p_dens[b-1]  # [1]
            p_edens[b]= p_edens[b-1] # [1]
        else:
            p_dens[b] = dens/dens0   # [1]
            p_edens[b]= denserror    # [1] #100/rbin would be artificial guess

    for b in range(gp.nipol):
        print(rbin[b], binmin[b], binmax[b], p_dens[b], p_edens[b], file=de) # [rscale], [dens0], [dens0]
        indr = (r<binmax[b])
        menclosed = 1.0*np.sum(indr)/totmass # /totmass for normalization to 1 at last bin #[totmass]
        merror = menclosed/np.sqrt(ab) # artificial menclosed/10 gives good approximation #[totmass]
        print(rbin[b], binmin[b], binmax[b], menclosed, merror, file=em) # [rscale], [totmass], [totmass]
    de.close()
    em.close()


    if not gpr.showplots: return
    ion(); subplot(111)
    print('rbin = ', rbin)
    print('p_dens = ', p_dens)
    print('p_edens = ', p_edens)

    plot(rbin, p_dens, 'b', linewidth=3)
    lbound = p_dens-p_edens; lbound[lbound<1e-6] = 1e-6
    ubound = p_dens+p_edens; 
    fill_between(rbin,lbound,ubound,alpha=0.5,color='r')
    # xscale('log'); 
    yscale('log')
    xlim([np.min(rbin),np.max(rbin)])
    ylim([np.min(lbound),np.max(ubound)])
    # ylim([1e-3,3.])#ylim([1e-6,2*np.max(p_dens)])
    # ylim([0,1])
    xlabel(r'$r [r_c]$')
    ylabel(r'$\nu(r)/\nu(0)$')
    # plt.legend(['\rho','\rho'],'lower left')
    # title(fil)
    # axes().set_aspect('equal')
    savefig(gpr.get_dens_png(0))
    if gpr.showplots:
        ioff(); show(); clf()

if __name__ == '__main__':
    # gpr.showplots = True
    import gl_params
    gp = gl_params.Params()
    run(gp)
