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
if i!=7:
    print("usage: gen_prof_gas.py sim xc yc zc r")
    print("example: gen_prof_gas.py 270 0.5 0.5 0.5 0.05")
    sys.exit()

rout = 1 # * rvir from AHF for maximal radial distance to include particles from
show = True
run  = True

dgas = mys.simdir()+"/ana/gas/"
snap = sys.argv[1].zfill(5)
outp = mys.d(snap)

xc = sys.argv[2]; yc = sys.argv[3]; zc = sys.argv[4]
r = str(float(sys.argv[5])*rout)
lma = sys.argv[6]

os.nice(0)

# gas
cmd = "get_sphere_gas -inp "+outp
cmd +=" -xc "+xc+" -yc "+yc+" -zc "+zc+" -rc "+r+" -out "+dgas+"rho_"+snap+".dat -typ 0 -lma "+lma
if(show): print(cmd)
if(run): os.system(cmd)

cmd = "get_prof_sph.py "+dgas+"rho_"+snap+".dat "+xc+" "+yc+" "+zc+" "+r+" > "+dgas+"prof_"+snap+".dat"
if(show): print(cmd)
if(run): os.system(cmd)

#cmd = "plot_prof_sph.py "+dgas+"prof_"+snap+".dat "+dgas+"rho_"+snap+".dat "+dgas+"prof_"+snap+".png "+xc+" "+yc+" "+zc
#if(show): print(cmd)
#if(run): os.system(cmd)
