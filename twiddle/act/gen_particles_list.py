#!/usr/bin/env python2

## \file
# read particles IDs of an AMR simulation snapshot, in a sphere around
# x y z within radius r for both dm and stars only

# (c) 2014 ETHZ, Pascal Steger, pascal@steger.aero

import sys, os
import lib.initialize as my
import lib.mysql as mys

os.nice(10)

show= True; run = True; loop=True

# choose: get_particles    for all particles
#         get_particles_dm for dm  particles only
gps = "get_particles"

i = len(sys.argv)
if(i!=2):
    print("usage: gen_particles_list.py snap")
    exit(1)
    
facr = 1 # * rvir from AHF for maximal radial distance to include particles
snap = int(sys.argv[1])
d    = mys.d(snap)
pdm = d+"particles_dm"
xl,yl,zl,ml,rl = mys.getxyzmr(snap,2)

nhalo= mys.get_nhalo(snap)
os.system("rm "+d+"particles_dm")
os.system("echo '"+str(nhalo)+"' >> "+pdm)

for i in range(len(xl)):
    xc=str(xl[i]); yc=str(yl[i]); zc=str(zl[i]); r=str(facr*rl[i])
    tmpdat = d+"tmp_"+str(i+1)+".dat"

    # all particles (DM/stars mixed)
    cmd1 = gps+" -inp "+d
    cmd1 += " -xc " + xc + " -yc " + yc + " -zc " + zc + " -r " + r
    cmd1 += ">" + tmpdat
    cmd2 = "true"
    my.threadif(cmd1,cmd2,run,run,show,run)
    if(not loop):
        break
    #TODO: stars only
