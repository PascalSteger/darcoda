#!/usr/bin/env python3

##
# @file
# move centered positions to spherical coordinates

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import gl_params as gp
import gr_params as gpr
import pdb
import random

for comp in range(gpr.ncomp):
    print('grh_Pos: input:')
    print(gpr.fileposcenter[comp])
    xall,yall,vlosall = np.loadtxt(gpr.fileposcenter[comp],
                                   comments='#', unpack=True)

    Rall = np.sqrt(xall**2+yall**2)
    rest = ( Rall <= gpr.Rmax )
    xall = xall[rest];  yall = yall[rest]; vlosall= vlosall[rest]

    n = len(xall)
    x=xall; y=yall; vlos=vlosall

    # old output to cartesian variables: needed anywhere?
    # print('grh_Pos: output:')
    # print(gpr.fileposcartesian)
    # fileposcartesian = open(gpr.fileposcartesian[comp], 'w')
    # print('# x [], y []',file=fileposcartesian)
    # for k in range(n):
    #     print(x[k], y[k], file=fileposcartesian)
    # fileposcartesian.close()

    # print(gpr.filevelcartesian)
    # filevelcartesian = open(gpr.filevelcartesian[comp], 'w')
    # print('# vlos [km/s]',file=filevelcartesian)
    # for k in range(n):
    #     print(vlos[k],file=filevelcartesian)
    # filevelcartesian.close()

    x = np.array(x);  y = np.array(y)
    R = np.sqrt(x**2+y**2)

    frac = np.sum((R <= gpr.Rmax))*1./len(R)
    print('fraction of particles nearer than rmax: ',frac*100,'%')

    # phi (azimuthal angle [-pi,pi])
    Phi   = np.zeros(n)
    for k in range(n):
        if x[k]>0:
            Phi[k] = np.arctan(y[k]/x[k])
        elif x[k]==0:
            Phi[k] = np.sign(y[k])*np.pi/2
        elif y[k]>=0:
            Phi[k] = np.arctan(y[k]/x[k])+np.pi
        elif y[k]<0:
            Phi[k] = np.arctan(y[k]/x[k])-np.pi

    print(gpr.fileposspherical[comp])
    fileposspherical = open(gpr.fileposspherical[comp], 'w')
    print('# R [pc]', 'Phi [rad]', 'vlos [km/s]', file=fileposspherical)
    for k in range(n):
        print(R[k], Phi[k], vlos[k], file=fileposspherical)
    fileposspherical.close()

    Rcore = R[len(R)/2]
    crcore = open(gpr.get_params_file(comp),'w')
    print('# Rcore in [pc],', 'surfdens_central (=dens0) in [munit/rcore**2]',\
          ', and in [munit/pc**2], and totmass [munit],',\
          ' and max(v_LOS) in [km/s]', file=crcore)
    print(Rcore, file=crcore)
    crcore.close()
