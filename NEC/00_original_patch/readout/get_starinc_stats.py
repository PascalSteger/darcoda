#!/usr/bin/python
'''Plot stars.  Also runs amr2cell to get density center.  To do so,
set RUNCELLS=1. Nothing is actually done with this information at the moment.'''
import sys
import os
import numpy
import scipy
import scipy.optimize
import pylab
import math
import toolKit2 as TK

fstart=265
fstop =266
RUNCELLS=0
DIR="./"
extension="SCI"
prefix="clusterROI1"
ZMAX=0
ZMIN=-7.
DZ=.1
BINN=(ZMAX-ZMIN)/DZ
NBINN=int(BINN)
# initial guess for location.
# search length
# offset for plotting
extent=0.003

pc=3.08e18

for i in xrange(fstart,fstop):
	number=repr(i).zfill(5)
	header=TK.getRamsesInfo(i,"../run/output_")
	ncpu,ndim,levelmin,levelmax,ngridmax,nstep_coarse,boxlen,time,aexp,H0,omega_m,omega_l,omega_k,omega_b,unit_l,unit_d,unit_t=header
	infile=DIR+prefix+number+".dat"
	f=open(infile,"r")
	print "#open file ",f
	mass=[]
	x=[]
	y=[]
	z=[]
	metal=[]
	xc_star=0.
	yc_star=0.
	zc_star=0.
	metal_bin_num=[]
	metal_bin_mass=[]
	metal_bin=[]
	for i in xrange(NBINN):
		metal_bin_num.append(0)
		metal_bin_mass.append(0.)
		metal_bin.append(ZMIN+i*DZ)
	ii=0
        for line in f:
                if line[0]=="#" or line[1]=="#":continue
                val = line.split()
		if int(val[6])==0:continue
                #if ii>0:
                #        if (val[1])==(x[ii-1]) and (val[2])==(y[ii-1]) and (val[3])==(z[ii-1]): continue
                mass.append(float(val[0]))
                x.append(((val[1])))
                y.append(((val[2])))
                z.append(((val[3])))
                metal.append(math.log10(max(abs(float(val[4])),1e-99)))
                if int(val[6])<0:metal[ii]=max(metal[ii]-math.log10(2e-2),ZMIN)
                else:metal[ii]=max(metal[ii]-math.log10(2e-2),ZMIN+1.)
		bin_indx=int(min((metal[ii]-ZMIN)/DZ,BINN-1))
		metal_bin_num[bin_indx]+=1
		metal_bin_mass[bin_indx]+=mass[ii]/1e3
		#print bin_indx,bin_indx*DZ+ZMIN,metal[ii]
		ii+=1

        n=len(mass)
        for ii in xrange(n):
                x[ii]=float(x[ii])
                y[ii]=float(y[ii])
                z[ii]=float(z[ii])
                xc_star+=x[ii]
                yc_star+=y[ii]
                zc_star+=z[ii]

	f.close()

	a=1.313
	b=-1.5
	c=.22
	func=[]
	for i in xrange(NBINN):
		func.append(a/((metal_bin[i]-b)**2+c**2))

	#fitfunc = lambda p, x: p[0]/((p[1]-x)**2+p[2]**2)
	fitfunc = lambda p, x: p[0]*numpy.exp(-(p[1]-x)**2/(2.*p[2]**2))
	errfunc = lambda p, x, y: fitfunc(p, x)-y
	p0=[1000.,-1.5,.3]
	p1,success = scipy.optimize.leastsq(errfunc, p0[:], args=(metal_bin,metal_bin_mass))
	print p1,success

	x=scipy.linspace(min(metal_bin),max(metal_bin),1000)
	pylab.figure(figsize=(5,5))
	b=pylab.bar(metal_bin,metal_bin_mass,width=DZ,bottom=0)
	l=pylab.plot(x,fitfunc(p1,x),color="red")
	pylab.xlabel("[Fe/H]")
	pylab.ylabel("1000 Msun")
	pylab.savefig("clusterROI100265.eps")


