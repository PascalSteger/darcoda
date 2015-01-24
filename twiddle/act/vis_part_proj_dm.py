#!/usr/bin/env python2

## \file
# Plot raw positions of DM particles

# (c) 2014 ETHZ, Pascal Steger, pascal@steger.aero

import sys, matplotlib
matplotlib.use('Agg')
from matplotlib import rc
rc('text',usetex=True)
#import pylab
from pylab import figure,xticks,yticks,grid,savefig
import matplotlib.pyplot as plt

i = len(sys.argv)
if i!=7:
    print("usage: vis_part_dm.py xc yc zc r infile outfile")
    print("ex: vis_part_dm.py 0.5 0.5 0.5 0.01 dm.dat mov/dm.png")
    exit(1)

xc = float(sys.argv[1]); yc = float(sys.argv[2]); zc = float(sys.argv[3]); r  = float(sys.argv[4])
infile=sys.argv[5]
outfile=sys.argv[6]

f = open(infile,'r')
x=[];y=[];z=[]
for line in f:
    values = line.split()
    xcheck = float(values[1])-xc
    ycheck = float(values[2])-yc
    zcheck = float(values[3])-zc
    if(abs(xcheck)<r and abs(ycheck) < r and abs(zcheck)<r):
        x.append(xcheck)
        y.append(ycheck)
        z.append(zcheck)
f.close()

fig = figure()
#xticks([-0.002,-0.001,0.0,0.001,0.002])
#yticks([-0.002,-0.001,0.0,0.001,0.002])
grid()
#ax  = fig.add_subplot(111)
#imshow(x,y)
print(len(x))
plt.scatter(x, y,alpha=0.3,marker='d')
plt.xlabel(r'r\quad[{\rm Mpc}/h]')
plt.ylabel(r'r\quad[{\rm Mpc}/h]')
savefig(outfile)
