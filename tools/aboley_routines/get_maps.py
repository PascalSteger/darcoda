#!/usr/bin/python
import sys
import os
import numpy
import toolKit2 as TK

istart=265
istop =266
box1=0.25
directory="MAPS/"
typ="5"
firstname="p_diry"
lma="-lma 12"
dir="y"
OUTDIR="../run/"
xc=0.964; yc=0.627; zc=0.8276
xc=0.962139;yc= 0.627026; zc=0.82346
xc=0.961838; yc=0.627493;zc= 0.820452
xc=0.5;yc=0.5;zc=0.5
#xc=0.509667385135;yc=0.507039940851;zc=0.493635728636 #ROI1

for i in xrange(istart,istop):
        number=repr(i).zfill(5)
#	f=open("sod_"+number+".dat","r")
#	for line in f:
#		val=line.split()
#		if 0.5<float(val[3])<0.6:
#			xc=float(val[3])
#			yc=float(val[4])
#			zc=float(val[5])
#	f.close()
	box1*=.5
	xmin=repr(xc-box1)
	xmax=repr(xc+box1)
	ymin=repr(yc-box1)
	ymax=repr(yc+box1)
	zmin=repr(zc-box1)
	zmax=repr(zc+box1)
	box1*=2.
	
        filename=OUTDIR+"output_"+number
	command='./amr2map -inp '+filename+' -zmin '+zmin+' -zmax '+zmax+' -ymin '+ymin+' -ymax '+ymax+' -xmin '+xmin+' -xmax '+xmax+' -dir '+dir+' -typ '+typ+' -out '+directory+firstname+number+".dat "+lma
	print command
	os.system(command)
