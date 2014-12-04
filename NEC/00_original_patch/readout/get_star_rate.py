#!/usr/bin/python
'''Plot stars.  Also runs amr2cell to get density center.  To do so,
set RUNCELLS=1. Nothing is actually done with this information at the moment.'''
import sys
import os
import numpy
import pylab
import math
import toolKit2 as TK

fstart=265
fstop=fstart+1
PRINTOUT=0
DIR="./"
extension=""
prefix="clusterROI1"
ZMAX=0
ZMIN=-7.
DZ=.1
BINN=(ZMAX-ZMIN)/DZ
NBINN=int(BINN)
# initial guess for location.
# search length
# offset for plotting
extent=0.001

pc=3.08e18


for i in xrange(fstart,fstop):
	number=repr(i).zfill(5)
	header=TK.getRamsesInfo(i,"../run/output_")
	ncpu,ndim,levelmin,levelmax,ngridmax,nstep_coarse,boxlen,time,aexp,H0,omega_m,omega_l,omega_k,omega_b,unit_l,unit_d,unit_t=header
	infile=DIR+prefix+number+".dat"
	f=open(infile,"r")
	outfile=DIR+prefix+number+extension+".png"
	print "#open file ",f
	for j in xrange(8):
		line=f.readline()
		if line.find("time")>0:
			key,val=line.split()
			simtime=float(val)
		if line.find("simu")>0:
			key,val=line.split("=")
			simage=float(val)*1e3 # 1e3 is to convert to Myr
	print "#found simage of ", simage
	mass=[]
	age=[]
	metal=[]
	isort=[]
	ii=0
        for line in f:
                if line[0]=="#" or line[1]=="#":continue
                val = line.split()
		if int(val[6])==0:continue
                mass.append(float(val[0]))
                metal.append(math.log10(max(abs(float(val[4])),1e-99)))
		print float(val[5])
		age.append(simage-float(val[5])*1e3)
		isort.append(ii)
		ii+=1

        n=len(mass)
	tmpage=age
	age_index=zip(tmpage,isort)
	age_index.sort()
	tmpage,isort=zip(*age_index)
	mtot=[]
	mtot.append(mass[isort[0]])
	print "# min max age",min(age),max(age)
	for ii in xrange(n):
		if(ii>0):
			mtot.append(mtot[ii-1]+mass[isort[ii]])
		


	f.close()
	
	if PRINTOUT:
		for ii in xrange(n):
			print tmpage[ii],mtot[ii]
	else:
		pylab.figure()
		#pylab.bar(metal_bin,metal_bin_mass,width=DZ,bottom=0)
		pylab.plot(tmpage,mtot)
		pylab.xlabel("Myr")
		pylab.ylabel("Solar Masses")
		pylab.savefig("star_rateROI100265.png")


