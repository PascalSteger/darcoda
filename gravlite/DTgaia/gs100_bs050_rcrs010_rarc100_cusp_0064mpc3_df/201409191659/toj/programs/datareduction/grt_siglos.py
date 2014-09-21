#!/usr/bin/env ipython3

##
# @file
# calculate velocity dispersion of 2D rings from a Walker dataset'''
# for triaxial systems
### TODO: extend method to work with elliptical bins // introduce tensors (cf Agnello)

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

from pylab import *
import numpy as np
import sys
import math
from BiWeight import meanbiweight
import gr_params as gpr
import gl_file as gfile
from gl_class_files import *
from gl_helper import bin_r_linear, bin_r_log, bin_r_const_tracers

def run(gp):
    # get radius, used for all binning
    print('input: ', gpr.get_com_file(0))
    if gfile.bufcount(gpr.get_com_file(0))<2:
        return
    x,y,vlos = np.loadtxt(gpr.get_com_file(0), skiprows=1, unpack=True) #2*[rscale], [km/s]
    totmass = 1.*len(x)  # [Munit], [Munit], where each star is weighted with the same mass
    r = np.sqrt(x*x+y*y) # [rscale]
    
    #set binning
    #gp.nipol = (max - min)*N^(1/3)/(2*(Q3-Q1)) #(method of wand)
    rmin = 0.                                       # [rscale]
    rmax = max(r) if gp.maxR < 0 else 1.0*gp.maxR # [rscale]
    
    if gpr.lograd:
        # space logarithmically in radius
        binmin, binmax, rbin = bin_r_log(rmax/gp.nipol, rmax, gp.nipol)
    elif gp.consttr:
        binmin, binmax, rbin = bin_r_const_tracers(r, len(r)/gp.nipol)
    else:
        binmin, binmax, rbin = bin_r_linear(rmin, rmax, gp.nipol)
        
    # offset from the start!
    rs = gpr.Rerr*np.random.randn(len(r))+r #[rscale]
    vlos = gpr.vrerr*np.random.randn(len(vlos))+vlos #[km/s]
    vfil = open(gp.files.sigfiles[0], 'w')
    print('r', 'sigma_r(r)', 'error', file=vfil)

    # 30 iterations for drawing a given radius in bin
    dispvelocity = np.zeros((gp.nipol,gpr.n))
    a = np.zeros((gp.nipol,gpr.n))
    p_dvlos = np.zeros(gp.nipol)
    p_edvlos = np.zeros(gp.nipol)
    
    for k in range(gpr.n):
        rsi = gpr.Rerr*np.random.randn(len(rs))+rs #[rscale]
        vlosi = gpr.vrerr*np.random.randn(len(vlos))+vlos #[km/s]
        for i in range(gp.nipol):
            ind1 = np.argwhere(np.logical_and(rsi>binmin[i],rsi<binmax[i])).flatten()
            a[i][k] = len(ind1) #[1]
            vlos1 = vlosi[ind1] #[km/s]
            if(len(ind1)<=1):
                dispvelocity[i][k] = dispvelocity[i-1][k]
                # attention! should be 0, uses last value
            else:
                dispvelocity[i][k] = meanbiweight(vlos1,ci_perc=68.4,\
                                                  ci_mean=True,ci_std=True)[1]
                # [km/s], see BiWeight.py
                
    for i in range(gp.nipol):
        dispvel = np.sum(dispvelocity[i])/gpr.n #[km/s]
        ab = np.sum(a[i])/(1.*gpr.n) #[1]
        if ab == 0:
            dispvelerr = p_edvlos[i-1] #[km/s]
            # attention! uses last error
        else:
            dispvelerr = dispvel/np.sqrt(ab) #[km/s]
        p_dvlos[i] = dispvel      #[km/s]
        p_edvlos[i]= dispvelerr #[km/s]

    maxsiglos = max(p_dvlos) #[km/s]
    print('maxsiglos = ',maxsiglos,'[km/s]')
    fpars = open(gp.files.get_scale_file(0),'a')
    print(maxsiglos, file=fpars)          #[km/s]
    fpars.close()
    import shutil
    shutil.copy2(gp.files.get_scale_file(0), gp.files.get_scale_file(1))

    for i in range(gp.nipol):
        #             [rscale]  [maxsiglos]                  [maxsiglos]
        print(rbin[i], binmin[i], binmax[i], np.abs(p_dvlos[i]/maxsiglos),np.abs(p_edvlos[i]/maxsiglos), file=vfil) #/np.sqrt(n))
    vfil.close()

    if gpr.showplots:
        gpr.show_plots_vlos(rbin, p_dvlos, p_edvlos)


if __name__ == '__main__':
    gpr.showplots = True
    import gl_params
    gp = gl_params.Params()
    run(gp)

