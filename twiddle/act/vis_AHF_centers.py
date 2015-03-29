#!/usr/bin/env python2

## \file
# This program plots an image of the dm distribution together with AHF
# generated halo positions.

# (c) 2014 ETHZ, Pascal Steger, pascal@steger.aero

import os, sys
import numpy as np
import array as ar
from math import log
import struct, matplotlib
matplotlib.use('Agg')
from matplotlib import rc
rc('text',usetex=True)
import pylab
import lib.mysql as mys

part2map="part2map_8byte"

# file-dependent variables.
VMIN=-.1
VMAX=3.0
RM=4
IMGSIZE='f'
factor=1e-30
datamin=1e-30

if RM == 4:
	SRM = "i"
elif RM == 8:
	SRM = "l"
else:
	print("Your choice of RM is unreasonable.")

i = len(sys.argv)

if i != 2:
	print("Incorrect Usage. Received ", i-1, " arguments.")
	print("$: vis_centers.py nsnap")
	sys.exit()

nsnap=int(sys.argv[1])
fname=mys.d(nsnap)+"/ov.dat"
cmd = part2map+" -inp "+mys.d(nsnap)+" -out "+fname+" -dir z -nx 1024 -ny 1024"
print(cmd)
os.system(cmd)

f = open(fname,'rb')
s = f.read(2*RM+struct.calcsize("2i"))
par,nx,ny,par = struct.unpack(SRM+"i"+"i"+SRM,s)

print("Found image size. ", nx,ny)
s = f.read(RM)
binvalues = ar.array(IMGSIZE)
binvalues.read(f,nx*ny)
data = np.array(binvalues,'float')
data = np.reshape(data,(ny,nx))

data2=np.zeros((ny,nx),'float')
for i in xrange(nx):
	for j in xrange(ny):
		if(abs(data[ny-j-1,i])<1E-90):
			data2[j,i]=0
		else:
			data2[j,i]=log(abs(data[ny-j-1,i]))

print(np.min(data2),np.max(data2))

pylab.figure()
pylab.xticks([])
pylab.yticks([])
pylab.imshow(data2,aspect=1,cmap='hot')
xc,yc,zc,m,r=mys.getxyzmr(nsnap,1)

for i in range(len(xc)):
	xs=xc[i]*nx;ys=ny-yc[i]*ny;s=r[i]*nx
	print(xs, ys, s)
	cir = pylab.Circle((xs, ys), 3.0*s, facecolor='none', edgecolor=(0.0,0,0.1), linewidth=1, alpha=0.5)
	pylab.gca().add_patch(cir)

pylab.savefig(fname+".png")

print(fname + ".png")
