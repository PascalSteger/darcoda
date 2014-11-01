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
fstop =fstart+1
RUNCELLS=0
DIR="./"
extension="SCI"
# initial guess for location.
# search length
# offset for plotting
extent=0.002
refine=extent/4.
mindf=1.75
# first guess
xc=0.547045824143; yc=0.512265650452; zc=0.458181241497
xc=0.429377013232; yc= 0.396090123252; zc=0.420656749876
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
	age=[]
	x=[]
	y=[]
	z=[]
	metal=[]
	xc_star=0.
	yc_star=0.
	zc_star=0.

	for j in xrange(8):
		line=f.readline()
		if line.find("time")>0:
			key,val=line.split()
			simtime=float(val)
		if line.find("simu")>0:
			key,val=line.split("=")
			simage=float(val)*1e3 # 1e3 is to convert to Myr
	print "#found simtime of ", simtime

	ii=0
	for line in f:
		if line[0]=="#" or line[1]=="#":continue
		val = line.split()

		if int(val[6])==0:continue

                mass.append(float(val[0]))
                x.append(((val[1])))
                y.append(((val[2])))
                z.append(((val[3])))
		metal.append(math.log10(max(abs(float(val[4]))/2e-2,1e-6)))
		age.append(-float(val[5])*1e3+simage)
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
	cluster_age=[]
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
			cluster_age.append(age[i])
			flag.append(0)
			xci+=x[i]
			yci+=y[i]
			zci+=z[i]

  	m=len(mass_cluster)
	xci/=float(m)
	yci/=float(m)
	zci/=float(m)
	for i in xrange(m):
		print cluster_age[i],cluster_metal[i]

