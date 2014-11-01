#!/usr/bin/python
''' Plot raw positions of DM particles'''

import os
import sys
import math
import numpy
import matplotlib
matplotlib.use('Agg')
import pylab
from pylab import *

i = len(sys.argv)
if i!=6:
    print "usage: vis_part_dm.py xc yc zc r filename"
    print "ex: vis_part_dm.py 0.5 0.5 0.5 0.01 dm.dat"
    exit(1)


xc = float(sys.argv[1])
yc = float(sys.argv[2])
zc = float(sys.argv[3])
r  = float(sys.argv[4])

f = open(sys.argv[5],'r')
x=[];y=[];z=[]
for line in f:
    values = line.split()
    x.append(float(values[2-1])-xc)
    y.append(float(values[3-1])-yc)
    z.append(float(values[4-1])-zc)
f.close()

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = Axes3D(fig)
    #fig.gca(projection='3d')
ax.scatter(x, y, z, zdir='z', label='zdir=z', color='blue', marker='o',alpha=0.7)

#ax.legend()
#plot(x,y,'b.')
#grid(True)
#axes().set_aspect('equal')
#TODO: fit to xc+-r/2

savefig('dm_sphere.png')
