#!/usr/bin/python
'''Plot profiles'''
import math
import sys
import os
import toolKit2 as TK

mH=1.67e-24
kpc=3.086e21
NR=100
mindf=3.1
OUTDIR="../run/"

i=len(sys.argv)
if i < 5 or (i > 5 and i < 8):
	print "Usage: get_profile.py filename input max_dist(kpc)  value [x0 y0 z0]"
	sys.exit()

f=open(sys.argv[1],"r")
index = 5+int(sys.argv[4])

array=[]
x=[]
y=[]
z=[]
dx=[]
rad_bin=[]
mass_bin=[]
mass=[]

header=TK.getRamsesInfo(int(sys.argv[2]),OUTDIR+"./output_")
ncpu,ndim,levelmin,levelmax,ngridmax,nstep_coarse,boxlen,time,aexp,H0,omega_m,omega_l,omega_k,omega_b,unit_l,unit_d,unit_t=header

array_conv=1.
if index >= 6:
	array_conv=unit_d/mH

amax=0.
i=0
for line in f:
	values=line.split()
	x.append(float(values[0])*unit_l/kpc)
	y.append(float(values[1])*unit_l/kpc)
	z.append(float(values[2])*unit_l/kpc)
	dx.append(float(values[3])*unit_l/kpc)
	array.append(float(values[index])*array_conv)
	mass.append(array[i]*mH*(dx[i]*kpc)**3)
	if array[i]>amax:
		amax=array[i];zm=z[i];ym=y[i];xm=x[i]
	i+=1

i=len(sys.argv)
if i > 5:
	x0=float(sys.argv[5]);y0=float(sys.argv[6]);z0=float(sys.argv[7])
	x0*=unit_l/kpc
	y0*=unit_l/kpc
	z0*=unit_l/kpc
else:
	x0=xm;y0=ym;z0=zm

print "#Using x0,y0,z0 ",x0,y0,z0,x0/unit_l*kpc,y0/unit_l*kpc,z0/unit_l*kpc

n=len(x)

max_dist=float(sys.argv[3])
DR=mindf/float(NR)
flag=[]
mass_save=[]
r_save=[]
rho_bin=[]
for i in xrange(n):
	dist=(math.sqrt( (x[i]-x0)**2+(y[i]-y0)**2+(z[i]-z0)**2 ))
	if dist <= max_dist:
		flag.append(1)
		r_save.append(dist)
		mass_save.append(mass[i]/2e39)
n=len(flag)
for j in xrange(NR+1):
	rad_bin.append(10.**(float(j)*DR-mindf)*max_dist)
	rad_bin[0]=0.
for j in xrange(NR):
	mass_bin.append(0.)
        #if j>0:mass_bin[j]+=mass_bin[j-1]
	for i in xrange(n):
		if flag[i]==1:
			if rad_bin[j]<r_save[i]<rad_bin[j+1]:
				mass_bin[j]+=mass_save[i]
				flag[i]=0
	rho_bin.append(mass_bin[j]*1e6/(4.*math.pi/3.*((1e3*rad_bin[j+1])**3-(1e3*rad_bin[j])**3))) # Msun/pc^3
mtot=0.
for j in xrange(NR):
	mtot+=mass_bin[j]
	print 1e3*rad_bin[j+1],mtot,rho_bin[j]

f.close()










