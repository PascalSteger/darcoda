#!/usr/bin/python
''' This program plots an image of a fortran image file.  The data are assumed to be float32, but that can
be changed by altering the IMGSIZE variable.  Record marker length is also taken to be 4, but that to
can be changed with the RM variable. Limits are set by VMIN and VMAX'''

import numpy
import array
import struct
import pylab
import sys

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

if i != 3:
	print "Incorrect Usage. Received ", i-1, " arguments."
	print "$: fortranimage.py infile log10[0|1]"
        sys.exit()

f = open(sys.argv[1],'rb')

s = f.read(2*RM+struct.calcsize("2i"))
par,nx,ny,par = struct.unpack(SRM+"i"+"i"+SRM,s)

print "Found image size. ", nx,ny

s = f.read(RM)
binvalues = array.array(IMGSIZE)
binvalues.read(f,nx*ny)
data=numpy.array(binvalues,'float')
data = numpy.reshape(data,(ny,nx))
if int(sys.argv[2]) == 1: 
	for i in xrange(nx):
        	for j in xrange(ny):
			data[j,i] = numpy.log10(max(data[j,i],datamin))

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
pylab.savefig(sys.argv[1]+".png")


