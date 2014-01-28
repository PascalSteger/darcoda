#!/usr/bin/env python3

## @file
# calculate density falloff of circular rings around center of mass

# (c) 2013 ETHZ Pascal Steger psteger@phys.ethz.ch

#start from data centered on COM already:
import sys
if(len(sys.argv)<2):
    print("use: densityfalloff.py [car,scl,sex,for]")
    exit(1)
dwarf=sys.argv[1]
dir = gp.files.machine+"/DTobs"
print(dir+dwarf+"/centerpos.txt")
import numpy as np
x,y,v = np.loadtxt(dir+dwarf+'/centerpos.txt',skiprows=1,usecols=(0,1,2),unpack=True)

totmass=len(x)

# calculate radius
r = np.sqrt(x*x+y*y)

# set number and size of (linearly spaced) bins
bins = 12
rmin = 0
rmax = max(r)
binlength=(rmax-rmin)/bins
binmin=np.zeros(bins)
binmax=np.zeros(bins)
rbin=np.zeros(bins)
for i in range(bins):
  binmin[i]=rmin+i*binlength
  binmax[i]=rmin+(i+1)*binlength
  rbin[i]=binmin[i]+0.5*binlength

#distance error: assuming constant error, up to now: from sim
rerror=0.1

#volume of a circular bin with dr=binlength
vol=np.zeros(bins)
for i in range(bins):
    vol[i] = np.pi*(binmax[i]**2-binmin[i]**2)

rs = rerror*np.random.randn(len(r))+r
d = open(dir+dwarf+'/densityfalloff.txt','w')
print('r','nu(r)/nu(0)','error', file=d)
d.close()

d = open(dir+dwarf+'/enclosedmass.txt','w')
print('r','M(<r)','error', file=d)
d.close()

#1000 iterations for getting random picked radius values
n = 1000
density=np.zeros((bins,n))
a=np.zeros((bins,n))
for k in range(n):
  rsi = rerror*np.random.randn(len(rs))+rs
  for i in range(bins):
    ind1 = np.argwhere(np.logical_and(rsi>binmin[i],rsi<binmax[i])).flatten()
    density[i][k] = 1.*totmass*len(ind1)/vol[i]
    a[i][k] = 1.*len(ind1)
    
dens0 = np.sum(density[0])/n
ab0 = np.sum(a[0])/n
denserr0 = dens0/np.sqrt(ab0)

p_dens  = np.zeros(bins)
p_edens = np.zeros(bins)
import math
for i in range(bins):
    dens = np.sum(density[i])/(1.*n)
    ab = np.sum(a[i])/(1.*n)
    if i==0:
        denserror = dens/dens0/10.
        p_dens[i] = dens/dens0
        p_edens[i]= denserror
    else:
        denserr = dens/np.sqrt(ab)
        denserror = np.sqrt((denserr/dens0)**2+(dens/(dens0**2)*denserr0)**2)
        if(math.isnan(denserror)):
            denserror=0
            p_dens[i]  = p_dens[i-1]
            p_edens[i] = p_edens[i-1]
        else:
            ## [PS]: TODO: change bin sizes to include same number of stars in each bin, not assigning wrong density as below
            p_dens[i]  = dens/dens0
            p_edens[i] = denserror

    d=open(dir+dwarf+'/densityfalloff.txt','a')
    print(rbin[i],p_dens[i],p_edens[i], file=d)
    d.close()
    
    indr = (r<binmax[i])
    menclosed = 1.0*np.sum(indr)/totmass
    merror = 10.0/totmass
    d=open(dir+dwarf+'/enclosedmass.txt','a')
    print(rbin[i],menclosed,merror, file=d)    # TODO: check: take rbinmax for MCMC?
    d.close()


from pylab import *
ion(); subplot(111)
print('rbin = ',rbin)
print('p_dens = ',p_dens)
print('p_edens = ',p_edens)

#linear
plot(rbin,p_dens,'b',linewidth=3)
lbound = p_dens-p_edens; lbound[lbound<1e-6] = 1e-6
ubound = p_dens+p_edens; 
fill_between(rbin,lbound,ubound,alpha=0.5,color='r')
#xscale('log')
yscale('log')

# vertical line
#axvline(x=4.0)
#ylim([0,1])
xlabel(r'$r [pc]$');    ylabel(r'$\nu(r)/\nu(0)$')
#plt.legend(['\rho','\rho'],'lower left')
title(dwarf)
#axes().set_aspect('equal')
savefig(dir+dwarf+"/densityfalloff.png")
ioff();show()

