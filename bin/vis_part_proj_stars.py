#!/usr/bin/python
''' Plot raw positions of star particles, with metallicity encoded in color, projected, and a contour plot of dm density in the background'''

import os
import sys
import math
from math import *
import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
rc('text',usetex=True)
import pylab
from pylab import *
import numpy as np
import matplotlib.pyplot as plt

i = len(sys.argv)
if i!=7:
    print "usage: vis_part_proj_stars.py x y z r instar indm"
    print "ex: vis_part_proj_stars.py 0.5 0.5 0.5 0.01 stars/part/stars_1.dat dm/part/dm_1.dat"
    exit(1)

xc = float(sys.argv[1]); yc = float(sys.argv[2]); zc = float(sys.argv[3])
r  = float(sys.argv[4])

fstar = open(sys.argv[5],'r')
m=[];x=[];y=[];z=[]; age=[]; metal=[]
for line in fstar:
    val = line.split()
    if(float(val[2-1])>xc):
        continue
    m.append(log(float(val[1-1])))
    x.append((float(val[2-1])-xc)*1000)
    y.append((float(val[3-1])-yc)*1000)
    z.append((float(val[4-1])-zc)*1000)
    metal.append((float(val[5-1])))
    age.append(float(val[6-1]))
fstar.close()

# sort such that small x values show up on top
x = array(x);y=array(y);z=array(z);metal=array(metal);age=array(age)
order = x.argsort()
x = x[order];y=y[order];z=z[order];metal=metal[order];age=age[order]

blowup = 1
mmi = min(metal); mma = max(metal)
for i in xrange(len(metal)):
    metal[i] = blowup*(metal[i]-mmi)/(mma-mmi)
print min(metal), max(metal)

fs   = 20
fig  = plt.figure()
ax = fig.add_subplot(111,aspect='equal')

fdm   = open(sys.argv[6],'r')
mdm=[];xdm=[];ydm=[];zdm=[];
for line in fdm:
    val = line.split()
#    if(float(val[2-1])>xc):
#        continue
    mdm.append(log(float(val[1-1])))
    xdm.append((float(val[2-1])-xc)*1000)
    ydm.append((float(val[3-1])-yc)*1000)
    zdm.append((float(val[4-1])-zc)*1000)
fdm.close()

ax.hexbin(ydm,zdm,bins='log',gridsize=200,cmap=cm.Reds)
ax.scatter(y,z,s=m,c=metal,vmin=0.0,vmax=blowup,alpha=1.0,cmap=cm.binary)

#scatter(ycm,zcm,"bo")
xlim([-2,2])
ylim([-2,2])
grid(True)

xlabel(r'$y\quad[{\rm kpc}/h]$',fontsize=fs); xticks(fontsize=fs)
ylabel(r'$z\quad[{\rm kpc}/h]$',fontsize=fs); yticks(fontsize=fs)

savefig('dm_sphere.png')
