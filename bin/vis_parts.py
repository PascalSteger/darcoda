#!/usr/bin/python
''' Plot raw positions of DM and star particles, with centers from AHF
and after shrinking sphere'''

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
import numpy as np
import matplotlib.pyplot as plt

i = len(sys.argv)
if i!=10:
    print "usage: vis_parts.py dm.dat stars.dat r xsp ysp zsp xahf yahf zahf"
    exit(1)


xsp = float(sys.argv[4]);     xahf = float(sys.argv[7])
ysp = float(sys.argv[5]);     yahf = float(sys.argv[8])
zsp = float(sys.argv[6]);     zahf = float(sys.argv[9])
r  = float(sys.argv[3])

fdm = open(sys.argv[1],'r')
xdm=[];ydm=[];zdm=[]
for line in fdm:
    values = line.split()
    xcheck = float(values[2-1])-xsp
    ycheck = float(values[3-1])-ysp
    zcheck = float(values[4-1])-zsp
    if(abs(xcheck)<r and abs(ycheck) < r and abs(zcheck)<r):
        xdm.append(xcheck)
        ydm.append(ycheck)
        zdm.append(zcheck)
fdm.close()

fs = open(sys.argv[2],'r')
xs=[];ys=[];zs=[]
for line in fs:
    values = line.split()
    xcheck = float(values[2-1])-xsp
    ycheck = float(values[3-1])-ysp
    zcheck = float(values[4-1])-zsp
    if(abs(xcheck)<r and abs(ycheck) < r and abs(zcheck)<r):
        xs.append(xcheck)
        ys.append(ycheck)
        zs.append(zcheck)
fs.close()

fig = pylab.figure()
ax  = fig.add_subplot(111,aspect='equal')
xticks([]); yticks([])

ax.scatter(xdm, ydm,alpha=0.3,marker='o',color='b')
ax.scatter(xs, ys,alpha=0.1,marker='o',color='g')
print xahf-xsp,yahf-ysp

a = fig.gca()

cir = Circle( (xsp-xsp,ysp-ysp), radius=0.00005,color='black')
a.add_patch(cir)

cir = Circle( (xahf-xsp,yahf-ysp), radius=0.00005,color='red')
a.add_patch(cir)

xticks([-0.002,-0.001,0.0,0.001,0.002],fontsize=20)
yticks([-0.002,-0.001,0.0,0.001,0.002],fontsize=20)
grid()
ax.set_xlabel(r'$r\quad[{\rm Mpc}/h]$',fontsize=20)
ax.set_ylabel(r'$r\quad[{\rm Mpc}/h]$',fontsize=20)
plt.savefig('dm_sphere.png')
