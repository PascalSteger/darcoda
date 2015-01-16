#!/usr/bin/python
'''plot profiles (DM/stars together)'''

import sys
import math
import matplotlib
matplotlib.use('Agg')
import pylab as pl
from pylab import *

# check
i=len(sys.argv)
if i!=3:
    print "Usage: plot_profs.py prof_dm prof_stars"
    print "Example: plot_profs.py prof_dm.dat prof_stars.dat"
    exit()

# DM: read in rad_bin, rho_bin
f=open(sys.argv[1],"r")
rad_bin_dm=[]; rho_bin_dm=[]; n_dm=[]
NR_dm = 0
for line in f:
    if(line[0]=="#"):
        continue
    val=line.split()
    rad_bin_dm.append(float(val[0]))
    rho_bin_dm.append(float(val[3]))
    n_dm.append(float(val[5]))
    NR_dm += 1
f.close()

# stars: read in rad_bin, rho_bin
f=open(sys.argv[2],"r")
rad_bin_stars=[]; rho_bin_stars=[]; n_stars=[]
NR_stars = 0
for line in f:
    if(line[0]=="#"):
        continue
    val=line.split()
    rad_bin_stars.append(float(val[0]))
    rho_bin_stars.append(float(val[3]))
    n_stars.append(float(val[5]))
    NR_stars += 1
f.close()

import math
ne_dm=[]
for i in xrange(NR_dm):
    if(n_dm[i-1]==0):
        ne_dm.append(0);
    else:
        ne_dm.append(rho_bin_dm[i-1]/math.sqrt(n_dm[i-1]))

ne_stars=[]
for i in xrange(NR_stars):
    if(n_stars[i-1]==0):
        ne_stars.append(0);
    else:
        ne_stars.append(rho_bin_stars[i-1]/math.sqrt(n_stars[i-1]))

# plot log-log
fig = pl.figure()
ax = fig.add_subplot(111)
ax.loglog(rad_bin_dm,rho_bin_dm,linewidth="2",c="b")
ax.loglog(rad_bin_stars,rho_bin_stars,linewidth="2",c="r")

# vertical line @ resolution limit 4pc given by A. Boley
# (CIC @ max level 14 for z~12, 1Mpc/h box, 256^3 equivalent part.
ax.axvline(x=4.0,c="b")

# vertical line @ radius where trelax = tsim
pc = 3.08e16                    #[m]
G = 6.67e-11                    #[m3/kg s2]

tsim = 0.466e9*365.25*24*3600   #[s]
lnl = 5                         #[1]

#M = m * 1.99e30                 #[kg]
#N = 29580                       #[1]

#r = (tsim*8*math.sqrt(G*M)*lnl/N)**(2./3.)/pc
#ax.axvline(x=r,c="g")

# error bars, assume Poissonian errors
ax.errorbar(rad_bin_dm[0:NR_dm],rho_bin_dm[0:NR_dm],yerr=ne_dm[0:NR_dm],c="b")
ax.errorbar(rad_bin_stars[0:NR_stars],rho_bin_stars[0:NR_stars],yerr=ne_stars[0:NR_stars],c="r")

# restrict display area
ax.set_xbound(3,300)
ax.set_ybound(0.01,100)

# legends
ax.legend(('dm','stars'),#'stars'),\
           'upper right',shadow=True)
ax.set_xlabel('$r$ [pc/h]')
ax.set_ylabel('$rho [M_\odot/pc^3]$') #TODO: units?

pl.savefig("map.dat.png")
