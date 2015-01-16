#!/usr/bin/python
'''plot all 3 projections of dm/stars/gas density, with centers'''

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
VMIN=-.1; VMAX=3.0
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
	sys.exit()

i = len(sys.argv)
if i != 5:
	print "$: vis_proj_dgs.py xc yc zc rvir"
        sys.exit()
	
xc = float(sys.argv[1]); yc = float(sys.argv[2]); zc = float(sys.argv[3])
rvir = float(sys.argv[4])

Ntick = 5; # must be odd (so we have center already indicated)
Nspace = Ntick - 1
BL = 2*rvir

xmi = xc-rvir; xma = xc+rvir; ymi = yc-rvir; yma = yc+rvir; zmi = zc-rvir; zma = zc+rvir
border = " -xmi " + str(xmi) + " -xma " + str(xma) + " -ymi " + str(ymi) + " -yma " + str(yma) + " -zmi " + str(zmi) + " -zma " + str(zma)

if True:
	os.system("part2map -inp /home/psteger/sci/sim/sim_ABoley/output_00270/ -out map_x.dat -dir x " + border)
	os.system("part2map -inp /home/psteger/sci/sim/sim_ABoley/output_00270/ -out map_y.dat -dir y " + border)
	os.system("part2map -inp /home/psteger/sci/sim/sim_ABoley/output_00270/ -out map_z.dat -dir z " + border)

	os.system("part2map -inp /home/psteger/sci/sim/sim_ABoley/output_00270/ -out map_x_star.dat -str true -dir x " + border)
	os.system("part2map -inp /home/psteger/sci/sim/sim_ABoley/output_00270/ -out map_y_star.dat -str true -dir y " + border)
	os.system("part2map -inp /home/psteger/sci/sim/sim_ABoley/output_00270/ -out map_z_star.dat -str true -dir z " + border)
	lma = " -lma 14 "
	os.system("amr2map -inp /home/psteger/sci/sim/sim_ABoley/output_00270/ -out map_x_amr.dat -dir x -typ 1 " + lma + border)
	os.system("amr2map -inp /home/psteger/sci/sim/sim_ABoley/output_00270/ -out map_y_amr.dat -dir y -typ 1 " + lma + border)
	os.system("amr2map -inp /home/psteger/sci/sim/sim_ABoley/output_00270/ -out map_z_amr.dat -dir z -typ 1 " + lma + border)

zoom = 1/(2*rvir); # 1Mpc big box/boxlength

# get ahf positions
# ahf=open("halosp","r")
ahf=open("halo","r")
xa = []; ya = []; za = []; Na = 0;
for line in ahf:
	if(line[0]=="#"):
		continue
	Na = Na + 1
	val=line.split()
	xa.append(float(val[2]))
	ya.append(float(val[3]))
	za.append(float(val[4]))
	
ahf.close()

fig=pylab.figure(figsize=(16,16),dpi=100)
fig.subplots_adjust(wspace=0.2,hspace=0.2)
fs=20

###############################################################################
#DM

###############################################################################
# along x
print "DM along x"
ax=pylab.subplot(3,3,1)
pylab.xticks([])
pylab.yticks([])

xt = []; yt = []
for i in arange(Ntick):
	xval = 1000*(-(i-Nspace/2)*BL/Nspace)
	yval = 1000*(-(i-Nspace/2)*BL/Nspace)
	xt.append("%.0f"%xval)
	yt.append("%.0f"%yval)

f = open("map_x.dat",'rb')
s = f.read(2*RM+struct.calcsize("2i"))
par,nx,ny,par = struct.unpack(SRM+"i"+"i"+SRM,s)

pylab.yticks(ny/Nspace*arange(Ntick),yt,fontsize=fs)
pylab.ylabel(r"$r\quad[{\rm kpc}/h]$",fontsize=fs)

s = f.read(RM)
binvalues = ar.array(IMGSIZE); binvalues.read(f,nx*ny)
data = numpy.array(binvalues,'float'); data = numpy.reshape(data,(ny,nx))
data2= numpy.zeros((ny,nx),'float')
f.close()

for i in xrange(nx):
	for j in xrange(ny):
		data2[j,i]=data[ny-j-1,i]
ax.imshow(data2,aspect=1,cmap=cm.hsv)

# loop over all ahf positions, draw crosses at positions
for i in range(Na):
	yp = (ya[i]-yc)*zoom+0.5 #0 for halosp, 2
	zp = 1-((za[i]-zc)*zoom+0.5) #1 for halosp, 3
	xp = xa[i] #2 for halosp, 4
	
	if(yp>0 and yp<1 and\
	   zp>0 and zp<1 and\
	   xp>xmi and xp<xma):# and ma>mmin):
		pylab.plot(yp*nx,zp*ny,'bo')


###############################################################################


###############################################################################
# along y
print "DM along y"
ax = pylab.subplot(3,3,2)
pylab.xticks([])
pylab.yticks([])
f = open("map_y.dat",'rb')
s = f.read(2*RM+struct.calcsize("2i"))
par,nx,ny,par = struct.unpack(SRM+"i"+"i"+SRM,s)

s = f.read(RM)
binvalues = ar.array(IMGSIZE); binvalues.read(f,nx*ny)
data = numpy.array(binvalues,'float'); data = numpy.reshape(data,(ny,nx))
data2= numpy.zeros((ny,nx),'float')
f.close()

for i in xrange(nx):
	for j in xrange(ny):
		data2[j,i]=data[ny-j-1,i]
pylab.imshow(data2,aspect=1,cmap=cm.hsv)

# loop over all ahf positions, draw crosses at positions
for i in range(Na):
	zp = (za[i]-zc)*zoom+0.5 #0 for halosp, 2
	zp = 1-zp;
	xp = 1-((xa[i]-xc)*zoom+0.5) #1 for halosp, 3
	xp = 1-xp;
	yp = ya[i] #2 for halosp, 4
	
	if(zp>0 and zp<1 and\
	   xp>0 and xp<1 and\
	   yp>ymi and yp<yma):# and ma>mmin):
		pylab.plot(xp*nx,zp*ny,'bo',lw=2)



###############################################################################
# along z
print "DM along z"
ax = pylab.subplot(3,3,3)
pylab.xticks([])
pylab.yticks([])
f = open("map_z.dat",'rb')
s = f.read(2*RM+struct.calcsize("2i"))
par,nx,ny,par = struct.unpack(SRM+"i"+"i"+SRM,s)

s = f.read(RM)
binvalues = ar.array(IMGSIZE); binvalues.read(f,nx*ny)
data = numpy.array(binvalues,'float'); data = numpy.reshape(data,(ny,nx))
data2= numpy.zeros((ny,nx),'float')
f.close()

for i in xrange(nx):
	for j in xrange(ny):
		data2[j,i]=data[ny-j-1,i]
pylab.imshow(data2,aspect=1,cmap=cm.hsv)

# loop over all ahf positions, draw crosses at positions
for i in range(Na):
	xp = (xa[i]-xc)*zoom+0.5 #0 for halosp, 2
	yp = 1-((ya[i]-yc)*zoom+0.5) #1 for halosp, 3
	zp = za[i] #2 for halosp, 4
	
	if(xp>0 and xp<1 and\
	   yp>0 and yp<1 and\
	   zp>zmi and zp<zma):# and ma>mmin):
		pylab.plot(xp*nx,yp*ny,'bo',lw=2)
		

################################################################################
#		stars

###############################################################################
# along x
print "stars along x"
pylab.subplot(3,3,4)
pylab.xticks([])
pylab.yticks([])
xt = []; yt = []
for i in arange(Ntick):
	yval = 1000*(zc+(i-Nspace/2)*BL/Nspace-zc)
	yt.append("%.0f"%yval)
pylab.yticks(ny/Nspace*arange(Ntick),yt,fontsize=fs)
pylab.ylabel(r"$r\quad[{\rm kpc}/h]$",fontsize=fs)

f = open("map_x_star.dat",'rb')
s = f.read(2*RM+struct.calcsize("2i"))
par,nx,ny,par = struct.unpack(SRM+"i"+"i"+SRM,s)

s = f.read(RM)
binvalues = ar.array(IMGSIZE); binvalues.read(f,nx*ny)
data = numpy.array(binvalues,'float'); data = numpy.reshape(data,(ny,nx))
data2= numpy.zeros((ny,nx),'float')
f.close()

for i in xrange(nx):
	for j in xrange(ny):
		data2[j,i]=data[ny-j-1,i]
pylab.imshow(data2,aspect=1,cmap=cm.hsv)

# loop over all ahf positions, draw crosses at positions
for i in range(Na):
	yp = (ya[i]-yc)*zoom+0.5 #0 for halosp, 2
	zp = 1-((za[i]-zc)*zoom+0.5) #1 for halosp, 3
	xp = xa[i] #2 for halosp, 4
	
	if(yp>0 and yp<1 and\
	   zp>0 and zp<1 and\
	   xp>xmi and xp<xma):# and ma>mmin):
		pylab.plot(yp*nx,zp*ny,'bo',lw=2)
###############################################################################
		
###############################################################################
# along y
print "stars along y"
pylab.subplot(3,3,5)
pylab.xticks([])
pylab.yticks([])

f = open("map_y_star.dat",'rb')
s = f.read(2*RM+struct.calcsize("2i"))
par,nx,ny,par = struct.unpack(SRM+"i"+"i"+SRM,s)

s = f.read(RM)
binvalues = ar.array(IMGSIZE); binvalues.read(f,nx*ny)
data = numpy.array(binvalues,'float'); data = numpy.reshape(data,(ny,nx))
data2= numpy.zeros((ny,nx),'float')
f.close()

for i in xrange(nx):
	for j in xrange(ny):
		data2[j,i]=data[ny-j-1,i]
pylab.imshow(data2,aspect=1,cmap=cm.hsv)

# loop over all ahf positions, draw crosses at positions
for i in range(Na):
	zp = (za[i]-zc)*zoom+0.5 #0 for halosp, 2
	zp = 1-zp;
	xp = 1-((xa[i]-xc)*zoom+0.5) #1 for halosp, 3
	xp = 1-xp;
	yp = ya[i] #2 for halosp, 4
	
	if(zp>0 and zp<1 and\
	   xp>0 and xp<1 and\
	   yp>ymi and yp<yma):# and ma>mmin):
		pylab.plot(xp*nx,zp*ny,'bo',lw=2)

###############################################################################
# along z
print "stars along z"
pylab.subplot(3,3,6)
pylab.xticks([])
pylab.yticks([])

f = open("map_z_star.dat",'rb')
s = f.read(2*RM+struct.calcsize("2i"))
par,nx,ny,par = struct.unpack(SRM+"i"+"i"+SRM,s)

s = f.read(RM)
binvalues = ar.array(IMGSIZE); binvalues.read(f,nx*ny)
data = numpy.array(binvalues,'float'); data = numpy.reshape(data,(ny,nx))
data2= numpy.zeros((ny,nx),'float')
f.close()

for i in xrange(nx):
	for j in xrange(ny):
		data2[j,i]=data[ny-j-1,i]
pylab.imshow(data2,aspect=1,cmap=cm.hsv)

# loop over all ahf positions, draw crosses at positions
for i in range(Na):
	xp = (xa[i]-xc)*zoom+0.5 #0 for halosp, 2
	yp = 1-((ya[i]-yc)*zoom+0.5) #1 for halosp, 3
	zp = za[i] #2 for halosp, 4
	
	if(xp>0 and xp<1 and\
	   yp>0 and yp<1 and\
	   zp>zmi and zp<zma):# and ma>mmin):
		pylab.plot(xp*nx,yp*ny,'bo',lw=2)

###############################################################################
## gas
print "gas along x"
# along x
pylab.subplot(3,3,7)
pylab.xticks([])
pylab.yticks([])
f = open("map_x_amr.dat",'rb')
s = f.read(2*RM+struct.calcsize("2i"))
par,nx,ny,par = struct.unpack(SRM+"i"+"i"+SRM,s)

for i in arange(Ntick):
	xval = 1000*(yc+(i-Nspace/2)*BL/Nspace-yc)
	yval = 1000*(zc+(i-Nspace/2)*BL/Nspace-zc)
	xt.append("%.0f"%xval)
	yt.append("%.0f"%yval)
pylab.xticks(nx/Nspace*arange(Ntick),xt,fontsize=fs)
pylab.yticks(ny/Nspace*arange(Ntick),yt,fontsize=fs)
pylab.xlabel(r"$r\quad[{\rm kpc}/h]$",fontsize=fs)
pylab.ylabel(r"$r\quad[{\rm kpc}/h]$",fontsize=fs)

s = f.read(RM)
binvalues = ar.array(IMGSIZE); binvalues.read(f,nx*ny)
data = numpy.array(binvalues,'float'); data = numpy.reshape(data,(ny,nx))
data2= numpy.zeros((ny,nx),'float')
f.close()

for i in xrange(nx):
	for j in xrange(ny):
		data2[j,i]=data[ny-j-1,i]
pylab.imshow(data2,aspect=1,cmap=cm.hsv)

# loop over all ahf positions, draw crosses at positions
for i in range(Na):
	yp = (ya[i]-yc)*zoom+0.5 #0 for halosp, 2
	zp = 1-((za[i]-zc)*zoom+0.5) #1 for halosp, 3
	xp = xa[i] #2 for halosp, 4
	
	if(yp>0 and yp<1 and\
	   zp>0 and zp<1 and\
	   xp>xmi and xp<xma):# and ma>mmin):
		pylab.plot(yp*nx,zp*ny,'bo',lw=2)

###############################################################################
# along y
print "gas along y"
pylab.subplot(3,3,8)
pylab.xticks([])
pylab.yticks([])
f = open("map_y_amr.dat",'rb')
s = f.read(2*RM+struct.calcsize("2i"))
par,nx,ny,par = struct.unpack(SRM+"i"+"i"+SRM,s)

for i in arange(Ntick):
	xval = 1000*((1-xc)+(i-Nspace/2)*BL/Nspace-xc)
	xt.append("%.0f"%xval)
pylab.xticks(nx/Nspace*arange(Ntick),xt,fontsize=fs)
pylab.xlabel(r"$r\quad[{\rm kpc}/h]$",fontsize=fs)

s = f.read(RM)
binvalues = ar.array(IMGSIZE); binvalues.read(f,nx*ny)
data = numpy.array(binvalues,'float'); data = numpy.reshape(data,(ny,nx))
data2= numpy.zeros((ny,nx),'float')
f.close()

for i in xrange(nx):
	for j in xrange(ny):
		data2[j,i]=data[ny-j-1,i]
pylab.imshow(data2,aspect=1,cmap=cm.hsv)

# loop over all ahf positions, draw crosses at positions
for i in range(Na):
	zp = (za[i]-zc)*zoom+0.5 #0 for halosp, 2
	zp = 1-zp;
	xp = 1-((xa[i]-xc)*zoom+0.5) #1 for halosp, 3
	xp = 1-xp;
	yp = ya[i] #2 for halosp, 4
	
	if(zp>0 and zp<1 and\
	   xp>0 and xp<1 and\
	   yp>ymi and yp<yma):# and ma>mmin):
		pylab.plot(xp*nx,zp*ny,'bo',lw=2)

###############################################################################
# along z
print "gas along z"
pylab.subplot(3,3,9)
pylab.xticks([])
pylab.yticks([])
f = open("map_z_amr.dat",'rb')
s = f.read(2*RM+struct.calcsize("2i"))
par,nx,ny,par = struct.unpack(SRM+"i"+"i"+SRM,s)

xt = []; yt = []
for i in arange(Ntick):
	xval = 1000*(xc+(i-Nspace/2)*BL/Nspace-xc)
	xt.append("%.0f"%xval)
pylab.xticks(nx/Nspace*arange(Ntick),xt,fontsize=fs)
pylab.xlabel(r"$r\quad[{\rm kpc}/h]$",fontsize=fs)

s = f.read(RM)
binvalues = ar.array(IMGSIZE); binvalues.read(f,nx*ny)
data = numpy.array(binvalues,'float'); data = numpy.reshape(data,(ny,nx))
data2= numpy.zeros((ny,nx),'float')
f.close()

for i in xrange(nx):
	for j in xrange(ny):
		data2[j,i]=data[ny-j-1,i]
pylab.imshow(data2,aspect=1,cmap=cm.hsv)

# loop over all ahf positions, draw crosses at positions
for i in range(Na):
	xp = (xa[i]-xc)*zoom+0.5 #0 for halosp, 2
	yp = 1-((ya[i]-yc)*zoom+0.5) #1 for halosp, 3
	zp = za[i] #2 for halosp, 4
	
	if(xp>0 and xp<1 and\
	   yp>0 and yp<1 and\
	   zp>zmi and zp<zma):# and ma>mmin):
		pylab.plot(xp*nx,yp*ny,'bo',lw=2)
###############################################################################

print "save image"
pylab.savefig("ov.png")
#cmd = "cp map.dat.png map_"+str(xc)+"_"+str(yc)+"_"+str(zc)+"_dirz_"+str(rvir)+"_"+str(mmin)+".png"
#os.system(cmd)
