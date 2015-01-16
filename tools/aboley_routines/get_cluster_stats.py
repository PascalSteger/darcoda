#!/usr/bin/python
'''Plot stars.  Also runs amr2cell to get density center.  To do so,
set RUNCELLS=1. Nothing is actually done with this information at the moment.'''
import sys
import os
import numpy
import pylab
import toolKit2 as TK

fstart=265
fstop =266
RUNCELLS=0
DIR="./"
OUTPUT_DIR="../run/"
# initial guess for location.
# search length
# offset for plotting
extent=0.004

pc=3.08e18

for i in xrange(fstart,fstop):
	number=repr(i).zfill(5)
	header=TK.getRamsesInfo(i,OUTPUT_DIR+"./output_")
	ncpu,ndim,levelmin,levelmax,ngridmax,nstep_coarse,boxlen,time,aexp,H0,omega_m,omega_l,omega_k,omega_b,unit_l,unit_d,unit_t=header
	infile=DIR+"star"+number+".dat"
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
	ii=0
	for line in f:
		if line[0]=="#" or line[1]=="#":continue
		val = line.split()
		if int(val[6])==0:continue
		#if float(val[4])<1e-8:continue
#                if ii>0:
#                        if (val[1])==(x[ii-1]) and (val[2])==(y[ii-1]) and (val[3])==(z[ii-1]): continue
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

	xci=x[0]
	yci=y[0]
	zci=z[0]
	mass_cluster=[]
	cluster_x=[]
	cluster_y=[]
	cluster_z=[]
	cluster_metal=[]
	id=[]
	mctot=0.
	total_mass=0.
	ztot=0.
	id_indx=0
	tmpx=0.
	tmpy=0.
	tmpz=0.
	count=0.
	for i in xrange(n):
		r=((x[i]-xci)**2+(y[i]-yci)**2+(z[i]-zci)**2)**.5
		if r<=extent:
			mctot+=mass[i]
			ztot+=metal[i]
			count+=1.
			tmpx+=x[i]
			tmpy+=y[i]
			tmpz+=z[i]
		else:
			mass_cluster.append(mctot)
			mctot=mass[i]
			id_indx+=1
			cluster_metal.append(ztot/count)
			cluster_x.append(tmpx/count)
			cluster_y.append(tmpy/count)
			cluster_z.append(tmpz/count)
			ztot=metal[i]
			tmpx=x[i]
			tmpy=y[i]
			tmpz=z[i]
			count=1.
		id.append(id_indx)	
		total_mass+=mass[i]
		xci=x[i];yci=y[i];zci=z[i]
		#print mass[i],x[i],y[i],z[i],id[i]
	mass_cluster.append(mctot)
	cluster_metal.append(ztot/count)
        cluster_x.append(tmpx/count)
        cluster_y.append(tmpy/count)
        cluster_z.append(tmpz/count)

	n=len(mass_cluster)
	print n
	mtot=0.
	flag=[]
	for i in xrange(n):
		flag.append(0)
	for i in xrange(n):
		for j in xrange(n):
			if i!=j:
				r=((cluster_x[i]-cluster_x[j])**2+(cluster_y[i]-cluster_y[j])**2+(cluster_z[i]-cluster_z[j])**2)**.5
				if r<extent:
					if not flag[j]:
						mass_cluster[i]+=mass_cluster[j]
						cluster_metal[i]=(cluster_metal[i]+cluster_metal[j])*.5
						cluster_x[i]=(cluster_x[i]+cluster_x[j])*.5
						cluster_y[i]=(cluster_y[i]+cluster_y[j])*.5
						cluster_z[i]=(cluster_z[i]+cluster_z[j])*.5
						flag[j]=1
		if not flag[i]: print i, mass_cluster[i],cluster_x[i],cluster_y[i],cluster_z[i],cluster_metal[i]
	print "# Total mass ",total_mass
							

