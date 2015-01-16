#!/usr/bin/python2.5
'''Plot stars.  Also runs amr2cell to get density center.  To do so,
set RUNCELLS=1. Nothing is actually done with this information at the moment.'''
import sys
import os
import numpy
import pylab
import math
import toolKit2 as TK

fstart=431
fstop =432
RUNCELLS=0
DIR="./STAR/"
extension="SCI"
NR=20
# initial guess for location.
# search length
# offset for plotting
extent=0.002
refine=extent/1.
mindf=1.4
# first guess
xc=0.547045824143; yc=0.512265650452; zc=0.458181241497
xc=0.541716549315; yc=0.509773442409; zc=0.461798114621
xc=0.534073102248; yc=0.506075144996; zc=0.467231232454

pc=3.08e18

for i in xrange(fstart,fstop):
	number=repr(i).zfill(5)
	header=TK.getRamsesInfo(i,"../output_")
	ncpu,ndim,levelmin,levelmax,ngridmax,nstep_coarse,boxlen,time,aexp,H0,omega_m,omega_l,omega_k,omega_b,unit_l,unit_d,unit_t=header
	infile=DIR+"star"+number+".dat"
	f=open(infile,"r")
	outfile=DIR+"star"+number+extension+".png"
	print "#open file ",f
	mass=[]
	x=[]
	y=[]
	z=[]
	metal=[]
	xc_star=0.
	yc_star=0.
	zc_star=0.
	ii=0
	for line in f:
		if line[0]=="#" or line[1]=="#":continue
		val = line.split()

		if int(val[6])==0:continue

                mass.append(float(val[0]))
                x.append(((val[1])))
                y.append(((val[2])))
                z.append(((val[3])))
                metal.append(float(val[4]))
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
        n=len(x)
	xc_star/=float(n)
	yc_star/=float(n)
	zc_star/=float(n)

	xci=0.
	yci=0.
	zci=0.
	mass_cluster=[]
	cluster_x=[]
	cluster_y=[]
	cluster_z=[]
	cluster_metal=[]
	flag=[]
	mtot=0.
	for i in xrange(n):
		r=((x[i]-xc)**2+(y[i]-yc)**2+(z[i]-zc)**2)**.5
		if r<=extent:
			mass_cluster.append(mass[i])
			mtot+=mass[i]
			cluster_metal.append(metal[i])
			cluster_x.append(x[i])
			cluster_y.append(y[i])
			cluster_z.append(z[i])
			flag.append(0)
			xci+=x[i]
			yci+=y[i]
			zci+=z[i]

		#xci=x[i];yci=y[i];zci=z[i]
		#print mass[i],x[i],y[i],z[i],id[i]
  	m=len(mass_cluster)
	xci/=float(m)
	yci/=float(m)
	zci/=float(m)
	DR=mindf/float(NR)
	rbin=[]
	mbin=[]
	dbin=[]
	for i in xrange(NR+1):
		rbin.append(refine*10.**(float(i)*(DR)-mindf)*unit_l/pc)
		mbin.append(0.)
	rbin[0]=0.
	for i in xrange(m):
		r=((cluster_x[i]-xci)**2+(cluster_y[i]-yci)**2+(cluster_z[i]-zci)**2)**.5*unit_l/pc
		for j in xrange(NR):
			if flag[i]==0 and (rbin[j]<r<rbin[j+1]):
				mbin[j]+=mass_cluster[i]
				flag[i]=1
	mtally=0.
	for j in xrange(NR):
		mtally+=mbin[j]
		dbin.append(mbin[j]/(rbin[j+1]**3-rbin[j]**3)*3./4./math.pi)
		print rbin[j+1],mtally,dbin[j],mtally/mtot


							

