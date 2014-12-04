#!/usr/bin/python
'''Plots surface temperature based on gas and pressure maps.'''
import numpy
import array
import struct
import pylab
import sys
import math
import toolKit2 as TK

# file-dependent variables.
RM=4
SMALL=1e-16
AUTOSCALE=1
VMAX=5.
VMIN=1.
gasConstant=8.254e7
OUTDIR="../run/"

RM=4
IMGSIZE='f'

if RM == 4:
        SRM = "i"
elif RM == 8:
        SRM = "l"
else:
        print "Your choice of RM is unreasonable."


i = len(sys.argv)

if i != 6:
        print "Incorrect Usage. Received ", i-1, " arguments."
        print "$: temp.py pressure gas number outfile log[0|1]"
        sys.exit()

NOUT=int(sys.argv[3])


header=TK.getRamsesInfo(NOUT,OUTDIR+"./output_")
ncpu,ndim,levelmin,levelmax,ngridmax,nstep_coarse,boxlen,time,aexp,H0,omega_m,omega_l,omega_k,omega_b,unit_l,unit_d,unit_t=header
print header

f=open(sys.argv[1],"rb")
h=open(sys.argv[2],"rb")

s = f.read(2*RM+struct.calcsize("2i"))
par,nx,ny,par = struct.unpack(SRM+"i"+"i"+SRM,s)

print "Found image size. ", nx,ny

s = h.read(2*RM+struct.calcsize("2i"))
par,nx2,ny2,par = struct.unpack(SRM+"i"+"i"+SRM,s)

print "Found image size. ", nx2,ny2

assert nx2==nx,ny2==ny

s = f.read(RM)
binvalues = array.array(IMGSIZE)
binvalues.read(f,nx*ny)
data=numpy.array(binvalues,'float')
data=numpy.reshape(data,(ny,nx))

s = h.read(RM)
binvalues1 = array.array(IMGSIZE)
binvalues1.read(h,nx*ny)
data1=numpy.array(binvalues1,'float')
data1=numpy.reshape(data1,(ny,nx))

f.close()
h.close()

for i in xrange(nx):
        for j in xrange(ny):
                data[j,i]*=1./data1[j,i]*(unit_l/unit_t)**2/gasConstant

for i in xrange(nx):
        for j in xrange(ny):
                data1[j,i]=data[ny-j-1,i]


if int(sys.argv[5]) == 1: data1 = numpy.log10(data1)
print numpy.min(data1),numpy.max(data1)

pylab.figure()
pylab.xticks([])
pylab.yticks([])
if AUTOSCALE==0:IM=pylab.imshow(data1,aspect=1,vmin=VMIN,vmax=VMAX)
if AUTOSCALE==1:IM=pylab.imshow(data1,aspect=1)
#pylab.colorbar(IM)
pylab.savefig(sys.argv[4])
