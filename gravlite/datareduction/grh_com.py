#!/usr/bin/env python3

##
# @file
# read simulation files (x,y,z, vx,vy,vz), center, cut to N particles, write 2D x,y,vlos

# (c) 2013 ETHZ, Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import random
import pdb
import gl_params as gp
import gr_params as gpr

# TODO: run() function

print('grh_com: input: ', gpr.simpos)
xall, yall, zall = np.loadtxt(gpr.simpos,skiprows=1,unpack=True) # 3*[ascale]
vxall,vyall,vzall= np.loadtxt(gpr.simvel,skiprows=1,unpack=True) # 3*[ascale]
nall = len(xall)                                                 # [1]

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

PM = [1. for i in trace]           # [1], constant, no probability of membership information from dataset
x  = [ xall[i]    for i in trace ] # [ascale]
y  = [ yall[i]    for i in trace ] # [ascale]
z  = [ zall[i]    for i in trace ] # [ascale]
vz = [ vzall[i]   for i in trace ] # [km/s]
PM = np.array(PM); x=np.array(x); y=np.array(y); z=np.array(z); vz=np.array(vz)

from gl_centering import *
com_x, com_y, com_z, com_vz = com_shrinkcircle_v(x,y,z,vz,PM) # 3*[ascale], [velocity]
print('COM [ascale]: ', com_x, com_y, com_z, com_vz)

# com_vz = np.sum(vz)/(1.*ndm) # TODO: include this in the shrinking algorithm

xnew = (x-com_x)*gp.ascale      # [pc]
ynew = (y-com_y)*gp.ascale      # [pc]
znew = (z-com_z)*gp.ascale      # [pc]
vznew = (vz-com_vz)*np.sqrt(gp.G1*gp.Mscale/gp.ascale) # [km/s], from conversion from system with L=G=M=1


R0 = np.sqrt(xnew**2+ynew**2)   # [pc]
R0.sort()                       # [pc]
Rhalf = R0[len(R0)/2]           # [pc]
Rcore = Rhalf                   # or gpr.r_DM # [pc]

xnew /= Rcore; ynew /= Rcore    # [rcore]

# only for 0 (all) and 1 (first and only population)
for comp in range(gpr.ncomp):
    crcore = open(gpr.get_params_file(comp),'w')
    print('# Rcore in [pc],',' surfdens_central (=dens0) in [munit/rcore**2],',\
          ' and in [munit/pc**2],',' and totmass [munit],',\
          ' and max(v_LOS) in [km/s]',file=crcore)
    print(Rcore,file=crcore)
    crcore.close()


    print('grh_com: output: ', gpr.fileposcenter[comp])
    filepos = open(gpr.fileposcenter[comp],'w')
    print('# x [Rcore]','y [Rcore]','vLOS [km/s]', file=filepos)
    for k in range(ndm):
        print(xnew[k], ynew[k], vznew[k], file=filepos)
    filepos.close()
    print('')
