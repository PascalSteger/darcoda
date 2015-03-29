#!/usr/bin/env python2

## \file
# get profile, write out to prof.dat

# (c) 2014 ETHZ, Pascal Steger, pascal@steger.aero

import sys
from numpy import array,sqrt

NBIN = 20

i=len(sys.argv)
if i!=6:
	print("Usage: get_prof_sph.py filename xc yc zc r (in Mpc/h)")
	print("example: get_prof_sph.py sph.dat 0.47 0.49 0.54 0.001")
	sys.exit()

f=open(sys.argv[1],"r")
xc=float(sys.argv[2]);yc=float(sys.argv[3]);zc=float(sys.argv[4]);
rc=float(sys.argv[5])

x=[]; y=[]; z=[]; r=[]; m=[]; rho=[]

for line in f:
	if line[0]=="#" or line[1]=="#":continue
	values=line.split()
	print(values)
	xch=float(values[1-1])-xc
	ych=float(values[2-1])-yc
	zch=float(values[3-1])-zc
	mch=float(values[4-1])
	rhoch=float(values[5-1])
	rch=sqrt(xch*xch+ych*ych+zch*zch)
	print(rch)
	if(rch<rc):
		x.append(xch);y.append(ych);z.append(zch);
		m.append(mch);rho.append(rhoch)
		r.append(rch)
f.close()
print("#",len(r),"particles")

# sort according radius
x = array(x);y = array(y);z = array(z);
m = array(m); r = array(r); rho=array(rho);

order = r.argsort()
r = r[order]
x = x[order]; y = y[order]; z = z[order];
rho = rho[order]
#print('#rmin',r[0])
#print('#rho',rho)

#logminr = log10(10./1e6)   # 10 pc
#logmaxr = log10(rc) # 2kpc
#rbin = np.logspace(logminr,logmaxr,num=NBIN,base=10.0)
##rbin = np.linspace(10**logminr,10**logmaxr,num=NBIN)
#rhobin=[]
#for i in range(len(rbin)-1):
#	tmp=[]
#	for j in range(len(r)):
#		if(rbin[i]<r[j] and r[j]<rbin[i+1]):
#			tmp.append(rho[j])
#	print(len(tmp))
#	rhobin.append(sum(tmp)/len(tmp))
#print(rhobin)

# use same number of particles in each bin
nstep=int(int(len(r)/NBIN+0.5)-0.5)
rmin=[]; rmax=[]; rhobin=[]; err=[]
for i in range(NBIN):
	tmp=[]
	for j in range(nstep):
		tmp.append(rho[i*nstep+j])
	rhot = sum(tmp)/len(tmp)
	rhobin.append(rhot)
	err.append(rhot/sqrt(len(tmp)))
	rmin.append(r[i*nstep])
	rmax.append(r[(i+1)*nstep])

#print('#rmin rmax rho err')
for i in range(len(rmin)):
	print(rmin[i], rmax[i], rhobin[i], err[i])

