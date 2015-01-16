#!/usr/bin/python
'''plot profile from sph-like densities'''

import sys
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import pylab
#from pylab import *
from matplotlib import rc
rc('text',usetex=True)# use latex in figures

# check
i=len(sys.argv)
if i!=7:
    print "Usage: plot_prof.py in_prof in_rho outfile xc yc zc"
    print "Example: plot_prof.py prof.dat rho.dat prof.png xc yc zc"
    exit(1)
    
infile  = sys.argv[1];  xc = float(sys.argv[4])
rhofile = sys.argv[2];  yc = float(sys.argv[5])
outfile = sys.argv[3];  zc = float(sys.argv[6])

# read in span in radii, over which rho is considered, rho, and Poisson error
f=open(infile,"r")
rmin=[];rmax=[];rho=[];err=[]
NR = 0
for line in f:
    if(line[0]=="#"):
        continue
    val=line.split()
    rmin.append(float(val[0])*1000)
    rmax.append(float(val[1])*1000)
    rho.append(float(val[2]))
    err.append(float(val[3]))
    NR += 1

rmean=[]; rspan=[]
for j in range(len(rmin)):
    rspan.append((rmax[j]-rmin[j])/2)
    rmean.append((rmax[j]-rmin[j])/2+rmin[j])


# read in particle positions with sph densities
rh=open(rhofile,"r")
pr=[]; prho=[]
for line in rh:
    val=line.split()
    x = float(val[0])-xc
    y = float(val[1])-yc
    z = float(val[2])-zc
    pr.append(math.sqrt(x*x+y*y+z*z)*1000)
    prho.append(float(val[4]))
#print pr[100:110]
#print prho[100:110]

# plot log-log
fig=plt.figure()
ax = fig.add_subplot(111)
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(pr,prho,'r.',linewidth=1)
ax.errorbar(rmean,rho,xerr=rspan,yerr=err,linewidth=3)


# vertical line @ resolution limit 4pc given by A. Boley
# (CIC @ max level 14 for z~12, 1Mpc/h box, 256^3 equivalent part.
#plt.axvline(x=4.0)

# visible region
#plt.xlim([10**-1,2*1])
#plt.ylim([10**0,10**2])
plt.xlabel(r'$r [{\rm kpc}/h]$')
plt.ylabel(r'$\rho [m_H/{\rm cm}^3]$')
#plt.legend(['\rho','\rho'],'lower left')
#plt.title('z=11.7')
plt.savefig(outfile)
