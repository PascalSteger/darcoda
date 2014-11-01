#!/usr/bin/python
'''not working!'''
''' Plot raw positions of DM particles'''

import os
import sys
import math
import matplotlib
matplotlib.use('Agg')
import pylab
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import *
from scipy.interpolate import griddata
import numpy.ma as ma
from numpy.random import uniform, seed

i = len(sys.argv)
if i!=6:
    print "usage: vis_part_dm.py xc yc zc r filename"
    print "ex: vis_part_dm.py 0.5 0.5 0.5 0.01 dm.dat"
    exit(1)

xc = float(sys.argv[1]); yc = float(sys.argv[2]); zc = float(sys.argv[3])
r  = float(sys.argv[4])

f = open(sys.argv[5],'r')
x=[];y=[];z=[];m=[]
for line in f:
    values = line.split()
    mcheck = float(values[1-1])
    xcheck = float(values[2-1])-xc
    ycheck = float(values[3-1])-yc
    zcheck = float(values[4-1])-zc
    if(abs(xcheck)<r and abs(ycheck) < r and abs(zcheck)<r):
        m.append(mcheck);x.append(xcheck); y.append(ycheck); z.append(zcheck)
f.close()


# make up some randomly distributed data
npts = 200
# define grid.
xi = np.linspace(-r,r,100)
yi = np.linspace(-r,r,100)
# grid the data.
zi = griddata((x, y), m, (xi[None,:], yi[:,None]), method='cubic')
# contour the gridded data, plotting dots at the randomly spaced data points.
CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
plt.colorbar() # draw colorbar
# plot data points.
plt.scatter(x,y,marker='o',c='b',s=5)
plt.xlim(-r,r)
plt.ylim(-r,r)
plt.title('griddata test (%d points)' % npts)
plt.show()

savefig('dm_sphere.png')
