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
fstop =266
RUNCELLS=0
DIR="./"
extension="SCI"
# initial guess for location.
# search length
# offset for plotting
extent=0.005
refine=extent/4.
# first guess
#ROI1
xc=0.509734634382; yc=0.507350727893; zc=0.493557386116
#ROI2
#xc=0.521898092018; yc=0.50365273389; zc=0.479760326627
#ROI3
#xc=0.5022598061; yc=0.498520087996; zc=0.496897785565

pc=3.08e18

for i in xrange(fstart,fstop):
	number=repr(i).zfill(5)
	header=TK.getRamsesInfo(i,"../run/output_")
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
	age=[]
	vx=[]
	vy=[]
	vz=[]
	id=[]
	ii=0
	buffer=""
	for line in f:
		if line[0]=="#" or line[1]=="#":
			if line.find("time")>0:buffer+=line
			if line.find("simu")>0:buffer+=line
			continue
		val = line.split()

		if int(val[6])==0:continue

                mass.append(float(val[0]))
                x.append((float(val[1])))
                y.append((float(val[2])))
                z.append((float(val[3])))
                metal.append(float(val[4]))
		age.append(float(val[5]))
		id.append(int(val[6]))
		vx.append(float(val[7]))
		vy.append(float(val[8]))
		vz.append(float(val[9]))
		ii+=1
	n=len(mass)

	f.close()
        n=len(x)

	xci=0.
	yci=0.
	zci=0.
	mass_cluster=[]
	cluster_x=[]
	cluster_y=[]
	cluster_z=[]
	cluster_metal=[]
	c_age=[]
	c_id=[]
	c_vx=[]
	c_vy=[]
	c_vz=[]
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
			c_age.append(age[i])
			c_id.append(id[i])
			c_vx.append(vx[i])
			c_vy.append(vx[i])
			c_vz.append(vx[i])
			flag.append(0)
  	m=len(mass_cluster)
	print buffer.rstrip('\n')
	for i in xrange(m):
		print mass_cluster[i],cluster_x[i],cluster_y[i],cluster_z[i],cluster_metal[i],c_age[i],c_id[i],c_vx[i],c_vy[i],c_vz[i]
	print "#total cluster mass, ",mtot


							

