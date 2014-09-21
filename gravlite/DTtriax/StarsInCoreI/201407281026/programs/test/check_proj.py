#!/usr/bin/env ipython3

##
# @file
# check deprojection and projection of nu

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import sys, pdb
import numpy as np
from scipy.stats import kurtosis
from pylab import *

import import_path as ip
ip.insert_sys_path('/home/psteger/sci/gravlite/programs/datareduction')

import gl_params as gp
import gr_params as gpr
import gl_file as gfile
import gl_helper as gp
import gl_project as glp


def run(gp):
    Rscale = []; Dens0Rscale = []; Dens0pc = []; Totmass = []; Maxsiglos = []
    rscale = []; dens0Rscale = []; dens0pc = []; totmass = []; maxsiglos = []

    for comp in range(3):
        A = np.loadtxt(gp.files.get_scale_file(comp), unpack=False, skiprows=1)
        Rscale.append(A[0])
        Dens0Rscale.append(A[1])
        Dens0pc.append(A[2])
        Totmass.append(A[3])
        
        B = np.loadtxt(gp.files.get_scale_file(comp)+'_3D', unpack=False, skiprows=1)
        rscale.append(B[0])
        dens0Rscale.append(B[1])
        dens0pc.append(B[2])
        totmass.append(B[3])
        
        print('#######  working on component ',comp)
        print('input: ',gpr.get_com_file(comp)+'_3D')
        # start from data centered on COM already:
        if gfile.bufcount(gpr.get_com_file(comp)+'_3D')<2: continue
        
        
        
        Rbin,Binmin,Binmax,Dens,Denserr = np.loadtxt(gp.files.Sigfiles[comp],\
                                                     skiprows=1,usecols=(0,1,2,3,4),\
                                                     unpack=True) # 3*[Rscale], [km/s]
        Rbin*=Rscale[comp]; Binmin*=Rscale[comp]; Binmax*=Rscale[comp]; Dens*=Dens0pc[comp]; Denserr*=Dens0pc[comp]
        
        
        
        rbin,binmin,binmax,dens,denserr = np.loadtxt(gp.files.Sigfiles[comp]+'_3D',\
                                                     skiprows=1,usecols=(0,1,2,3,4),\
                                                     unpack=True) # 3*[Rscale], [km/s]
        rbin*=rscale[comp]; binmin*=rscale[comp]; binmax*=rscale[comp]; dens*=dens0pc[comp]; denserr*=dens0pc[comp]
        
        
        ion()
        f=figure(figsize=(6,3))
        ax1 = f.add_subplot(121)                         # Sig
        ax2 = f.add_subplot(122, sharex=ax1)             # nu
        
        ax1.plot(Rbin, Dens,'b',lw=1)
        lbound = Dens-Denserr; lbound[lbound<1e-6] = 1e-6
        ubound = Dens+Denserr; 
        ax1.fill_between(Rbin,lbound,ubound,alpha=0.5,color='r')
        ax1.set_yscale('log')
        ax1.set_xlim([0,np.max(Binmax)])
        ax1.set_ylim([np.min(lbound),np.max(ubound)])
        ax1.set_xlabel(r'$R [R_c]$')
        ax1.set_ylabel(r'$\nu_{2D}(R)/\nu_{2D}(0)$')
        
        try:
            ax1.plot(Rbin, rho_INT_Rho(Rbin, dens,denserr))
            ax1.plot(Rbin, rho_INT_Rho(Rbin, Rho_INT_rho(Rbin,Dens,Denserr),denserr))
        except Exception as detail:
            print('rho_INT_Rho giving NaN in plotting')
        draw()
        
        
        ax2.plot(rbin, dens,'b',lw=1)
        lbound = dens-denserr; lbound[lbound<1e-6] = 1e-6
        ubound = dens+denserr;
        ax2.fill_between(rbin,lbound,ubound,alpha=0.5,color='r')
        
        ax2.set_yscale('log')
        ax2.set_xlim([0,np.max(binmax)])
        ax2.set_ylim([np.min(lbound),np.max(ubound)])
        ax2.set_xlabel(r'$r [R_c]$')
        ax2.set_ylabel(r'$\nu(r)/\nu(0)$')
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        draw()
        
        # projSig = rho_INT_Rho(rbin, dens)
        # projSig = test(rbin, binmin, binmax, dens)
        # ax1.plot(rbin, projSig)
        
        ax2.plot(rbin,Rho_INT_rho(Rbin,Dens,Denserr),color='green')
        draw()
        
        ioff(); show()
## run(gp)
# check for projection and deprojection artefacts
# @param gp global parameters


if __name__=="__main__":
    import gl_params
    gp = gl_params.Params()

    run(gp)

    
