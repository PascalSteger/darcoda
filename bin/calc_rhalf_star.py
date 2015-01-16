#!/usr/bin/python
# find radius including half of all star particles

import sys
import os
import numpy as np
from numpy import *
import math
import mys
import initialize as my

i=len(sys.argv)
if(i!=3):
    print "find_rhalf_star.py snap hid"
    exit(0)

snap=int(sys.argv[1])
hid=int(sys.argv[2])

#print "missing file"
if(not mys.exists_snap(snap)):
    print "snapshot "+str(snap)+" missing"
    exit(0)

xc,yc,zc,rvir=mys.getxyzrstars_hid(snap,hid)

file_name=mys.d(snap)+"stars/stars_"+str(hid)+".dat"
if(os.path.getsize(file_name)==0):
    print "empty file"
    exit()

halo = my.open_file(file_name,"r")
x=[];y=[];z=[]; N=0
for line in halo:
    N=N+1
    val=line.split()
    x.append(float(val[1]))
    y.append(float(val[2]))
    z.append(float(val[3]))

x = array(x); y = array(y); z = array(z)
r  = np.sqrt((x-xc)**2+(y-yc)**2+(z-zc)**2)
order = r.argsort()
r = r[order]
# ASS: all star particles have the same mass
mys.set_rhalf_star(snap,hid,r[N/2])
