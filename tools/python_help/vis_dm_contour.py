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
if i!=7:
    print "usage: vis_dm_contour.py xc yc zc r infile outfile"
    print "ex: vis_dm_contour.py 0.5 0.5 0.5 0.01 dm.dat dm.png"
    exit(1)

infile = sys.argv[5]; outfile=sys.argv[6]
xc = float(sys.argv[1]); yc = float(sys.argv[2]); zc = float(sys.argv[3])
r  = float(sys.argv[4])

f = open(infile,'r')
x=[];y=[];z=[];m=[];rho=[];h=[]
for line in f:
    values = line.split()
    xcheck = float(values[1-1])-xc
    ycheck = float(values[2-1])-yc
    zcheck = float(values[3-1])-zc
    mcheck = float(values[4-1])
    rhocheck=float(values[5-1])
    hcheck = float(values[6-1])
    if(abs(xcheck)<r and abs(ycheck) < r and abs(zcheck)<r/10.0):
        m.append(mcheck);x.append(xcheck); y.append(ycheck); z.append(zcheck)
        rho.append(rhocheck);h.append(hcheck);
f.close()
print min(rho), max(rho), mean(rho), median(rho)

# make up some randomly distributed data

# define grid.
ngrid = 30
xi = linspace(-r,r,ngrid)
yi = linspace(-r,r,ngrid)

# grid the data.
zi = griddata((x, y), rho, (xi[None,:], yi[:,None]))#, method='cubic')

#print zi

# contour the gridded data, plotting dots at the randomly spaced data points.
nstep = 6
CS = plt.contour(xi,yi,zi,nstep,linewidths=0.5,colors='k')
CS = plt.contourf(xi,yi,zi,nstep,cmap=plt.cm.jet)
plt.colorbar() # draw colorbar
# plot data points.
plt.scatter(x,y,marker='o',c='b',s=1)
plt.xlim(-r,r)
plt.ylim(-r,r)
plt.title('rho (%d points)' % len(x))
plt.show()

savefig(outfile)
