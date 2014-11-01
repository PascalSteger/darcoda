#!/usr/bin/python
''' This program plots an image of the dm distribution together 
with sod generated substructure positions. Usage: fim.py xc yc zc rvir'''

import os
import sys
import numpy
import array as ar
from array import *
import struct
import matplotlib
matplotlib.use('Agg')
import pylab
from pylab import *

# file-dependent variables.
#VMIN=-12
#VMAX=-3
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
	print "Your choice of RM is unreasonable."

i = len(sys.argv)

if i != 6:
	print "Incorrect Usage. Received ", i-1, " arguments."
	print "$: fim.py xc yc zc rvir mmin"
	print "where xc,yc,zc give the position of centering, in Mpc/h"
	print " and rvir is half of the box size for zoom, in kpc/h"
	print " and mmin is minimal mass for AHF halos, in Msun/h"
        sys.exit()

xc = float(sys.argv[2-1])
yc = float(sys.argv[3-1])
zc = float(sys.argv[4-1])
rvir = float(sys.argv[5-1])
mmin = float(sys.argv[6-1])

cmd = "part2map -inp /home/psteger/sci/run/output_00270/ -out map.dat -dir z"
xmi = xc-rvir
xma = xc+rvir
ymi = yc-rvir
yma = yc+rvir
zmi = zc-rvir
zma = zc+rvir
cmd = cmd + " -xmi " + str(xmi) + " -xma " + str(xma)
cmd = cmd + " -ymi " + str(ymi) + " -yma " + str(yma)
cmd = cmd + " -zmi " + str(zmi) + " -zma " + str(zma)
print cmd
os.system(cmd)

f = open("map.dat",'rb')
s = f.read(2*RM+struct.calcsize("2i"))
par,nx,ny,par = struct.unpack(SRM+"i"+"i"+SRM,s)

print "Found image size. ", nx,ny
s = f.read(RM)
binvalues = ar.array(IMGSIZE)
binvalues.read(f,nx*ny)
data = numpy.array(binvalues,'float')
data = numpy.reshape(data,(ny,nx))
#if int(sys.argv[2]) == 1: 
#	for i in xrange(nx):
#        	for j in xrange(ny):
#			data[j,i] = numpy.log10(max(data[j,i],datamin))

data2=numpy.zeros((ny,nx),'float')
for i in xrange(nx):
	for j in xrange(ny):
		data2[j,i]=data[ny-j-1,i]

print numpy.min(data2),numpy.max(data2)

pylab.figure()
pylab.xticks([])
pylab.yticks([])
#pylab.imshow(data2,aspect=1,vmin=VMIN,vmax=VMAX)
pylab.imshow(data2,aspect=1)
# loop over all sod positions, draw circles with radius proportional mass
f=open("xcm/sod_nr.dat","r")
zoom = 1/(2*rvir); # 1Mpc big box/boxlength
for line in f:
        values=line.split()
	xs=(((float(values[4-1]))-xc)*zoom+0.5)*nx
	ys=(0.5-((float(values[5-1]))-yc)*zoom)*ny
	zs=float(values[6-1])
	s=float(values[7])*zoom*nx
	if(xs-s>0 and ys-s>0 and xs+s<nx and ys+s<ny and zs>zmi and zs<zma):
		cir = pylab.Circle((xs, ys), s, facecolor='none', edgecolor=(0.9,0,0.9), linewidth=1, alpha=0.5)
		pylab.gca().add_patch(cir)

h=0.719;
# loop over all ahf positions, draw crosses at positions
print "#plot only ahf positions for subhalos with mvir>"+str(mmin)+" Msun/h"
#ahf=open("xcm/ahf_halos","r")
ahf=open("halo","r")
for line in ahf:
	if(line[0]=="#"):
		continue
	values=line.split()
	x = float(values[2-1]); #3 for ahf_halos
	y = float(values[3-1]); #3 for ahf_halos
	z = float(values[4-1]); #3 for ahf_halos

	xa=((x-xc)/zoom+0.5)*zoom*nx
	ya=((y-yc)/zoom+0.5)*zoom*ny
	za=z

	#ma=float(values[9-1])
	#ra=float(values[10-1])
	if(xa>0 and ya>0 and xa<nx and ya<ny and za>zmi and za<zma):# and ma>mmin):
		pylab.plot([xa],[ya],'rx',lw=2)
	
pylab.savefig("map.dat.png")
#cmd = "cp map.dat.png map_"+str(xc)+"_"+str(yc)+"_"+str(zc)+"_dirz_"+str(rvir)+"_"+str(mmin)+".png"
#os.system(cmd)
