#!/usr/bin/env python2

## \file
# given a halo position, do:
# - get_sphere
# - get_prof
# - plot_prof
# for dm only

# (c) 2014 ETHZ, Pascal Steger, pascal@steger.aero

import sys, os
import lib.mysql as mys

i=len(sys.argv)
if i!=6:
    print("usage: gen_prof_dm.py sim xc yc zc r")
    print("example: gen_prof_dm.py 270 0.5 0.5 0.5 0.05")
    sys.exit()

rout = 1 # * rvir from AHF for maximal radial distance to include particles from
show = True
run  = True

ddm = mys.simdir()+"/ana/dm/"
snap = sys.argv[1].zfill(5)
outp = mys.d(snap)

xc = sys.argv[2]; yc = sys.argv[3]; zc = sys.argv[4]
r = str(float(sys.argv[5])*rout)

os.nice(0)

# Dark Matter
cmd = "get_sphere_dm -inp "+outp
cmd +=" -xc "+xc+" -yc "+yc+" -zc "+zc+" -rc "+r+"> "+ddm+"dm_"+snap+".dat"
cmd += " && octreef "+ddm+"dm_"+snap+".dat > "+ddm+"rho_"+snap+".dat"
if(show): print(cmd)
if(run): os.system(cmd)

cmd = "get_prof_sph.py "+ddm+"rho_"+snap+".dat "+xc+" "+yc+" "+zc+" "+r+" > "+ddm+"prof_"+snap+".dat"
if(show): print(cmd)
if(run): os.system(cmd)

#cmd = "plot_prof_sph.py "+ddm+"prof_"+snap+".dat "+ddm+"rho_"+snap+".dat "+ddm+"prof_"+snap+".png "+xc+" "+yc+" "+zc
#if(show): print(cmd)
#if(run): os.system(cmd)
