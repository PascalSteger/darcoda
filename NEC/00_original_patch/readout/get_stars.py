#!/usr/bin/python
'''Get stars'''
import sys
import os
import numpy

fstart=270
fstop =271
OUTDIR="../run/"


xc = .5
yc = .5
zc = .5
box1=1.

box1*=.5
xmin=repr(xc-box1)
xmax=repr(xc+box1)
ymin=repr(yc-box1)
ymax=repr(yc+box1)
zmin=repr(zc-box1)
zmax=repr(zc+box1)
box1*=2.

for i in xrange(fstart,fstop):
	number=repr(i).zfill(5)
	infile=OUTDIR+"output_"+number
	outfile="star"+number+".dat"
	command='./p2gnuplot -inp '+infile+' -zmin '+zmin+' -zmax '+zmax+' -ymin '+ymin+' -ymax '+ymax+' -xmin '+xmin+' -xmax '+xmax+' -str true >'+outfile 
	print command
	os.system(command)
