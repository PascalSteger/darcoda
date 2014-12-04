#!/usr/bin/python
'''Plot profiles'''
import math
import sys
import os
import toolKit2 as TK
import pylab as P

mH=1.67e-24
kpc=3.08e21
pc=kpc*1e-3
NR=100
#mindf=3.1
mindf=2.5
OUTDIR="../run/"

i=len(sys.argv)
if i < 4 or (i > 4 and i < 7):
	print "Usage: get_profile.py filename input max_dist(kpc)  [x0 y0 z0]"
	sys.exit()

f=open(sys.argv[1],"r")

array=[]
x=[]
y=[]
z=[]
dx=[]
rad_bin=[]
mass_bin=[]

header=TK.getRamsesInfo(int(sys.argv[2]),OUTDIR+"/output_")
ncpu,ndim,levelmin,levelmax,ngridmax,nstep_coarse,boxlen,time,aexp,H0,omega_m,omega_l,omega_k,omega_b,unit_l,unit_d,unit_t=header

rho0=3./8./6.67e-8/math.pi*(H0*1e5/kpc/1e3)**2*(omega_m-omega_b)/2e33*(3.08e18)**3
print "# rho0 ", rho0

array_conv=1e-6 # Mega solar mass units

amax=0.
i=0
for line in f:
	if line[0]=="#" or line[1]=="#":continue
	values=line.split()
	x.append(float(values[1])*unit_l/kpc)
	y.append(float(values[2])*unit_l/kpc)
	z.append(float(values[3])*unit_l/kpc)
	array.append(float(values[0])*array_conv)
	if array[i]>amax:
		amax=array[i];zm=z[i];ym=y[i];xm=x[i]
	i+=1

i=len(sys.argv)
if i > 5:
	x0=float(sys.argv[4]);y0=float(sys.argv[5]);z0=float(sys.argv[6])
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
mass_bin=[]
rad_bin=[]
for i in xrange(n):
	dist=(math.sqrt( (x[i]-x0)**2+(y[i]-y0)**2+(z[i]-z0)**2 ))
	if dist <= max_dist:
		flag.append(1)
		r_save.append(dist)
		mass_save.append(array[i])
n=len(flag)
print "# number of particles",n
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
	rho_bin.append((mass_bin[j])*1e6/(4.*math.pi/3.*((1e3*rad_bin[j+1])**3-(1e3*rad_bin[j])**3)))  # Msun/pc^3
mtot=0.
for j in xrange(NR):
	mtot+=mass_bin[j]
	print 1e3*rad_bin[j+1],rad_bin[j+1]*kpc/unit_l,mtot,rho_bin[j],mtot/(4.*math.pi/3.*(1e3*rad_bin[j+1])**3)*1e6/rho0*aexp**3

P.figure()
P.semilogy(rad_bin[0:NR],rho_bin[0:NR])
P.show()

f.close()


