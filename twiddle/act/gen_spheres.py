#!/usr/bin/env python2

## \file
# read particles of an AMR simulation snapshot, in a sphere around x y
# z within radius r for both dm and stars only

# (c) 2015 ETHZ, Pascal Steger, pascal@steger.aero

import sys
import os
import pdb
os.nice(1)
import lib.initialize as my
import lib.mysql as mys

show = True
run = True
loop = True

gsd = "get_sphere_dm"
gss = "get_sphere_stars"

i = len(sys.argv)
if(i!=3):
    print("usage: gen_spheres.py snap typ")
    exit(1)

facr = 1 # * rvir from AHF for maximal radial distance to include particles

snap = sys.argv[1]
typ = int(sys.argv[2])
sim = mys.d(snap)

# get xc,yc,zc,mvir,rvir
xcl, ycl, zcl, mvirl,rvirl = mys.getxyzmr(snap, typ)
xsl, ysl, zsl, rsl = mys.getxyzrstars(snap, typ)
pdb.set_trace()
#print('select ',xcl[0])
#print(len(xcl))

for i in range(3):
    xc=str(xcl[i]);       xs=str(xsl[i]);
    yc=str(ycl[i]);       ys=str(ysl[i]);
    zc=str(zcl[i]);       zs=str(zsl[i]);
    r=str(rvirl[i]*facr); rs=str(rsl[i]*facr*0.9);
    sj=str(i+1)

    # Dark Matter only
    cmd  = gsd+" -inp "+sim
    cmd +=" -xc "+xc+" -yc "+yc+" -zc "+zc+" -rc "+r
    cmd +=">"+sim+"dm/dm_"+sj+".dat && "
    cmd +="octreef "+sim+"dm/dm_"+sj+".dat > "+sim+"dm/rho_"+sj+".dat"
    if(show): print(cmd)
    if(run): my.thread(cmd)

    if(mys.is_dmonly):
        continue

    # stars only
    cmd = gss+" -inp "+sim
    # center on star center, but with virial radius
    if(typ==1):
        cmd +=" -xc "+xc+" -yc "+yc+" -zc "+zc+" -rc "+r
    if(typ==2): # second step: use shrinked sphere value for centering
        cmd +=" -xc "+xs+" -yc "+ys+" -zc "+zs+" -rc "+r
    cmd +=">"+sim+"stars/stars_"+sj+".dat"
    if(show): print(cmd)
    if(run): my.thread(cmd)

    if(not loop):
        break # for debug purposes only: break after first halo
