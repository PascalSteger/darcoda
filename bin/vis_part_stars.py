#!/usr/bin/python
''' Plot raw 3D positions of star particles'''
import os
import sys
import math
import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
rc('text',usetex=True)
import pylab
from pylab import *

i = len(sys.argv)
if i!=6:
    print "usage: vis_part_stars.py x y z r file_star_pos file_dm_pos"
    print "ex: vis_part_stars.py 0.5 0.5 0.5 0.001 stars/part/stars_1.dat dm/part/dm_1.dat"
    exit(1)

xc = float(sys.argv[1])
yc = float(sys.argv[2])
zc = float(sys.argv[3])
r  = float(sys.argv[4])

fstar = open(sys.argv[5],'r')
fdm   = open(sys.argv[6],'r')
x=[];y=[];z=[]; age=[]; metal=[]
for line in fstar:
    val = line.split()
    x.append(float(val[2-1]))
    y.append(float(val[3-1]))
    z.append(float(val[4-1]))
    metal.append(float(val[5-1]))
    age.append(float(val[6-1]))
fstar.close()

#print metal,age
xlabel('y',fontsize=fs)
ylabel('z',fontsize=fs)
title('all stars, blue: low metal, red: high metal')

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.gca(projection='3d')

mmi = min(metal); mma = max(metal)
for i in xrange(len(metal)):
    metal[i] = (metal[i]-mmi)/(mma-mmi)
ax.scatter(x, y, z, zdir='z', label='zdir=z', color=metal, marker='o')
#ax.legend()
#ax.set_xlim3d(xc-r/2,xc+r/2)
#ax.set_ylim3d(yc-r/2,yc+r/2)
#ax.set_zlim3d(zc-r/2,zc+r/2)

savefig('star_sphere.png')
