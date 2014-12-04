#!/usr/bin/python
'''Get dark matter particles. Plot using plot_dm.py. Reads in from sod output.
Finds appropriate halo based on user-supplied criteria.'''
import sys
import os
import numpy
import toolKit2 as TK

istart=270
istop =271
box1=0.1
firstname="DM/dmROI"
#box1=.01
#firstname="DM/dmROI"

xc=0.96212;yc= 0.627266;zc= 0.822716
xc=0.961838; yc=0.627493; zc=0.820452

xc=.5;yc=.5;zc=.5

#xc=0.512181801766; yc=0.497567898833; zc=0.492664528769
#xc=.465907748587; yc=0.549363715653; zc=0.488315598565
#xc=0.503872566719; yc=0.510072762492; zc=0.500204699594
#xc=0.522065456497; yc=0.517368513465; zc=0.508104226162
#xc=0.500877356013; yc=0.497868332451;zc= 0.498453450227
#xc=0.51159960012; yc= 0.499975848709; zc=0.479329563805

for i in xrange(istart,istop):
        number=repr(i).zfill(5)
	box1*=.5
	xmin=repr(xc-box1)
	xmax=repr(xc+box1)
	ymin=repr(yc-box1)
	ymax=repr(yc+box1)
	zmin=repr(zc-box1)
	zmax=repr(zc+box1)
	box1*=2.
	
        outfile=firstname+number+".dat"
	infile="../run/output_"+number
        command='./p2gnuplot -per true -inp '+infile+' -zmin '+zmin+' -zmax '+zmax+' -ymin '+ymin+' -ymax '+ymax+' -xmin '+xmin+' -xmax '+xmax+' >'+outfile 
	print command
	os.system(command)
