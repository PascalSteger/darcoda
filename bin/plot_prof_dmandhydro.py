#!/usr/bin/python
'''plot profile from sph-like densities, for both dmonly and hydro sims'''

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
if i!=5:
    print "Usage: plot_prof.py infile1 infile2 outfile z"
    print "Example: plot_prof.py prof_dm.dat prof_hydro.dat prof.png 10.1"
    exit(1)
    
infile1 = sys.argv[1]
infile2 = sys.argv[2]
outfile = sys.argv[3]
z=sys.argv[4]

# read in span in radii, over which rho is considered, rho, and Poisson error
f1=open(infile1,"r")
rmin1=[];rmax1=[];rho1=[];err1=[]
NR1 = 0
for line in f1:
    if(line[0]=="#"):
        continue
    val=line.split()
    rmin1.append(float(val[0])*1000)
    rmax1.append(float(val[1])*1000)
    rho1.append(float(val[2])*0.826) # /((o_m-o_b)/o_m)
    err1.append(float(val[3]))
    NR1 += 1
f1.close()

f2=open(infile2,"r")
rmin2=[];rmax2=[];rho2=[];err2=[]
NR2 = 0
for line in f2:
    if(line[0]=="#"):
        continue
    val=line.split()
    rmin2.append(float(val[0])*1000)
    rmax2.append(float(val[1])*1000)
    rho2.append(float(val[2]))
    err2.append(float(val[3]))
    NR2 += 1
f2.close()

rmean1=[]; rspan1=[]
for j in range(len(rmin1)):
    rspan1.append((rmax1[j]-rmin1[j])/2)
    rmean1.append((rmax1[j]-rmin1[j])/2+rmin1[j])
rmean2=[]; rspan2=[]
for j in range(len(rmin2)):
    rspan2.append((rmax2[j]-rmin2[j])/2)
    rmean2.append((rmax2[j]-rmin2[j])/2+rmin2[j])
    
# plot log-log
fig=plt.figure()
ax = fig.add_subplot(111)
ax.set_xscale('log')
ax.set_yscale('log')

ax.errorbar(rmean1,rho1,xerr=rspan1,yerr=err1,label='DM only',linewidth=2)
ax.errorbar(rmean2,rho2,xerr=rspan2,yerr=err2,label='full hydro',linewidth=2,color='red')
plt.legend()

# vertical line @ resolution limit 4pc given by A. Boley
# (CIC @ max level 14 for z~12, 1Mpc/h box, 256^3 equivalent part.
ax.axvline(x=0.004)

# visible region
plt.xlim([10**-1,2*10**0])
plt.ylim([10**0,10**2])
plt.grid(False)
plt.xlabel(r'$r\quad[{\rm kpc}/h]$',size=20)
plt.ylabel(r'$\rho\quad[m_H/{\rm cm}^3]$',size=20)
plt.title(r'z = '+z)
plt.savefig(outfile)
