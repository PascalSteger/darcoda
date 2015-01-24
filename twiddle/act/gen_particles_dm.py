#!/usr/bin/env python2

## \file
# read particles IDs of an AMR simulation snapshot, in a sphere around
# x y z within radius r for both dm and stars only

# (c) 2014 ETHZ, Pascal Steger, pascal@steger.aero

import sys, os
import lib.mysql as mys

os.nice(10)

show= True; run = True; loop=True

i = len(sys.argv)
if(i!=2):
    print("usage: gen_particles_dm.py snap")
    exit(1)
    
snap = int(sys.argv[1])
d=mys.d(snap)
pdm=d+"particles_dm"
nhalo= mys.get_nhalo(snap)
os.system("rm "+d+"particles_dm")
os.system("echo '"+str(nhalo)+"' >> "+pdm)
for c in range(nhalo):
    tmpdat = d+"tmp_"+str(c+1)+".dat"
    cmd = "echo $(wc -l "+tmpdat+"|cut -d' ' -f1) >> "+pdm+";"
    cmd = cmd + "cat "+tmpdat+">>"+pdm+";"
    cmd = cmd + "rm "+tmpdat
    if(show): print(cmd)
    if(run): os.system(cmd)
    if(not loop): break
    # TODO: stars only
