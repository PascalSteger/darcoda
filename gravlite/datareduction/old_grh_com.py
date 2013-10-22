#!/usr/bin/env python3

##
# @file
# old version of grh_com

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import random
import gl_params as gp
import gr_params as gpr

print('input: ', gpr.simpos)
xall, yall, zall = np.loadtxt(gpr.simpos,skiprows=1,unpack=True)
vxall,vyall,vzall= np.loadtxt(gpr.simvel,skiprows=1,unpack=True)
nall = len(xall)

# shuffle and restrict to ntracer random points
if gpr.ntracers1>0:
  ndm = min(gpr.ntracers1,nall-1)
  trace = random.sample(range(nall), nall)
else:
  ndm = min(gpr.ntracers2,nall-1)
  trace = random.sample(range(nall), nall)
if gpr.ntracers1+gpr.ntracers2 == 0:
  ndm = nall
  trace = np.arange(nall)


x  = [ xall[i]    for i in trace ]
y  = [ yall[i]    for i in trace ]
z  = [ zall[i]    for i in trace ]
vz = [ vzall[i]   for i in trace ]

centreofmassx = np.sum(x)/(1.*ndm)
centreofmassy = np.sum(y)/(1.*ndm)
centreofmassz = np.sum(z)/(1.*ndm)
comvz = np.sum(vz)/(1.*ndm)

print(centreofmassx,centreofmassy,centreofmassz)

xnew = (x-centreofmassx)*gp.ascale
ynew = (y-centreofmassy)*gp.ascale
znew = (z-centreofmassz)*gp.ascale
vznew = vz-comvz


R0 = np.sqrt(xnew**2+ynew**2)   # [pc]
R0.sort()                       # [pc]
Rhalf = R0[len(R0)/2]           # [pc]
Rcore = Rhalf                   # or gpr.r_DM # [pc]

xnew /= Rcore; ynew /= Rcore

# only for 0 (all) and 1 (first and only population)
for comp in range(gpr.ncomp):
    crcore = open(gpr.get_params_file(i),'w')
    print('# Rcore in [pc], surfdens_central (=dens0) in [munit/rcore**2], and in [munit/pc**2], and totmass [munit], and max(v_LOS) in [km/s]', file=crcore)
    print(Rcore, file=crcore)
    crcore.close()


    print('output: ', gpr.fileposcenter[comp])
    filepos = open(gpr.fileposcenter[comp],'w')
    print('x [Rcore]','y [Rcore]','z [Rcore]','vz [km/s]', file=filepos)
    for k in range(ndm):
        print(xnew[k], ynew[k], znew[k], vznew[k], file=filepos)
    filepos.close()
