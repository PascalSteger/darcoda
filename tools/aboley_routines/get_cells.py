#!/usr/bin/python
'''Runs amr2cell for producing cell distributions.'''
import sys
import os
import numpy
import pylab
import toolKit2 as TK

fstart=1
fstop =266
OUTDIR="../run/"
# initial guess for location.
# search length
xc=.5
yc=.5
zc=.5
box1=.25

pc=3.08e18
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
	amr2cell_outfile="cellROI"+number+".dat"
	print "Using amr2cell for output "+number+"."
	command='./amr2cell -inp '+OUTDIR+'output_'+number+' -zmi '+zmin+' -zma '+zmax+' -ymi '+ymin+' -yma '+ymax+' -xmi '+xmin+' -xma '+xmax+' -out '+amr2cell_outfile
	os.system(command)
