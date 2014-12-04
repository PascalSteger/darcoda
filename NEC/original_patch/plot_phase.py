#!/usr/bin/python
'''Calculate the radial distribution from any data cube.'''

import numpy
import array
import struct
import pylab
import math
import sys
import toolKit2 as TK

# file-dependent variables.
istart=33
istop =34
OUTDIR="../run/"

RM=4
mp=1.67e-24
kb=1.3801e-16
USECENTERING=0

for i in xrange(istart,istop):
	number=repr(i).zfill(5)
	filename="cell"+number+".dat"
	f = open(filename,'r')
	header=TK.getRamsesInfo(i,OUTDIR+"./output_")
	ncpu,ndim,levelmin,levelmax,ngridmax,nstep_coarse,boxlen,time,aexp,H0,omega_m,omega_l,omega_k,omega_b,unit_l,unit_d,unit_t=header

	norm=unit_d/mp

	dd=[]
	tk=[]
	h2=[]

	l=0
	ddMax=0.
	offSet=0.
	for line in f:
		values=line.split()
		dd.append(float(values[6]))
		tk.append(float(values[10]))
		tk[l]=tk[l]/dd[l]*(unit_l/unit_t)**2*mp/kb
		h2.append(float(values[19]))
		dd[l]*=norm
		l+=1
        f.close()
	filename="phase"+number+".png"         

	print numpy.min(dd),numpy.max(dd)
	print numpy.min(tk),numpy.max(tk)

	pylab.figure()
	pylab.ylabel("T/mu")
	pylab.xlabel("#/cc")
	CP=pylab.loglog(dd,tk,'.',markersize=0.5)
        pylab.axis([1e-4,1e6,1,1e6])
	pylab.savefig(filename)
        

