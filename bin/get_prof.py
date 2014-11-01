#!/usr/bin/python
'''get profile, write out to prof.dat'''
import math
import sys
import os
import toolKit2 as TK
import matplotlib

mH=1.67e-24
kpc=3.0857e21 #cm
pc=kpc*1e-3   #cm
NR=30 # number of bins
#mindf=3.1 or 2.5
OUTDIR="/scratch/psteger/sim_aboley/"

i=len(sys.argv)
if i!=7:
	print "Usage: get_prof.py filename input_Nr xc yc zc r (in Mpc/h)"
	print "example: get_prof.py dm.dat 50 0.47 0.49 0.54 0.001"
	sys.exit()

f=open(sys.argv[1],"r")

array=[]
x=[]; y=[]; z=[];
dx=[]
rad_bin=[]; mass_bin=[]

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
		amax=array[i]; xm=x[i]; ym=y[i]; zm=z[i]
	i+=1

x0=float(sys.argv[3]); y0=float(sys.argv[4]); z0=float(sys.argv[5]); r_max=float(sys.argv[6])
x0*=unit_l/kpc
y0*=unit_l/kpc
z0*=unit_l/kpc
r_max*=unit_l/kpc

print "#Using x0,y0,z0, r_max ",x0,y0,z0,r_max,x0/unit_l*kpc,y0/unit_l*kpc,z0/unit_l*kpc,r_max/unit_l*kpc

n=len(x)

mindf=3.1;
DR=mindf/float(NR) # step size
flag=[]
mass_save=[]
r_save=[]
rho_bin=[]; mass_bin=[]; rad_bin=[]; n_bin=[]
for i in xrange(n):
	dist=(math.sqrt( (x[i]-x0)**2+(y[i]-y0)**2+(z[i]-z0)**2 ))
	if dist <= r_max:
		flag.append(1)
		r_save.append(dist)
		mass_save.append(array[i])
n=len(flag)
print "# number of particles",n
for j in xrange(NR+1):
	rad_bin.append(10.**(float(j)*DR-mindf)*r_max)
	rad_bin[0]=0.
for j in xrange(NR):
	mass_bin.append(0.)
	n_bin.append(0)
        #if j>0:mass_bin[j]+=mass_bin[j-1]
	for i in xrange(n):
		if flag[i]==1:
			if rad_bin[j]<r_save[i]<rad_bin[j+1]:
				mass_bin[j]+=mass_save[i]
				n_bin[j]=n_bin[j]+1
				flag[i]=0
	rho_bin.append((mass_bin[j])*1e6/(4.*math.pi/3.*((1e3*rad_bin[j+1])**3-(1e3*rad_bin[j])**3)))  # Msun/pc^3
mtot=0.
print "# rad_bin, rad_bin, mtot, rho_bin, rho_inside, n_bin"
for j in xrange(NR):
	mtot+=mass_bin[j]
	print 1e3*rad_bin[j+1],rad_bin[j+1]*kpc/unit_l,mtot,rho_bin[j],mtot/(4.*math.pi/3.*(1e3*rad_bin[j+1])**3)*1e6/rho0*aexp**3,n_bin[j]


f.close()


