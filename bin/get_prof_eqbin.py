#!/usr/bin/python
'''get profile with equal number of particles per bin, write out to prof.dat'''
import math
import sys
import os
import toolKit2 as TK
import matplotlib
import numpy

NR=24 # number of bins
OUTDIR="/scratch/psteger/sim/"
mH=1.67e-24
kpc=3.0857e21 #cm
pc=kpc*1e-3   #cm

i=len(sys.argv)
if i!=7:
	print "Usage: get_prof.py filename input_Nr xc yc zc r (in Mpc/h)"
	print "example: get_prof.py dm.dat 50 0.47 0.49 0.54 0.001"
	sys.exit()

f=open(sys.argv[1],"r")

array=[]
x=[]; y=[]; z=[]
dx=[]

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
x0*=unit_l/kpc; y0*=unit_l/kpc; z0*=unit_l/kpc; r_max*=unit_l/kpc

print "#Using x0,y0,z0, r_max ",x0,y0,z0,r_max

n=len(x)

flag=[]
r_save=[];m_save=[]
rho_bin=[]; mass_bin=[]; rad_bin=[]; n_bin=[]
for i in xrange(n):
	dist=(math.sqrt( (x[i]-x0)**2+(y[i]-y0)**2+(z[i]-z0)**2 ))
	if dist <= r_max:
		flag.append(1)
		r_save.append(dist)
		m_save.append(array[i])

n=len(flag)
print "# number of particles",len(x)
print "# number of particles inside rvir",n
r_bin=[]; m_bin=[]

# sort, increasing radius
order = numpy.argsort(r_save)

# get bin radii
# determine mean radius on loglog scale, volume, density, errors
#for j in range(n):
#	print r_save[order[j]]

for j in xrange(NR):
	r_bin.append(r_save[order[int(1.0*n/NR*j)]])
	m_bin.append(0.0)
	for i in xrange(n/NR):
		m_bin[len(m_bin)-1] += m_save[order[i+n/NR*j]]
		
	rho_bin.append((j+1)*m_bin[j]/(4.*math.pi/3.*(r_bin[j]**3)))

for j in range(NR):
	print math.log10(r_bin[j]), math.log10(m_bin[j]), math.log10(rho_bin[j])


exit(0)

mtot=0.
print "# rad_bin, rad_bin, mtot, rho_bin, rho_inside, n_bin"
for j in xrange(NR):
	mtot+=mass_bin[j]
	print 1e3*rad_bin[j+1],rad_bin[j+1]*kpc/unit_l,mtot,rho_bin[j],mtot/(4.*math.pi/3.*(1e3*rad_bin[j+1])**3)*1e6/rho0*aexp**3,n_bin[j]


f.close()


