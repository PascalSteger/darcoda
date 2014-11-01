#!/usr/bin/python
'''plot profiles (all DM together for 5 most massive halos)'''

import sys
import math
import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
rc('text',usetex=True)
import pylab as pl
from pylab import *


# check syntax
i=len(sys.argv)
if i!=1:
    print "Usage: plot_profs_dm_all.py"
    print "Example: plot_profs_dm_all.py"
    exit()


# prepare figure
fig = pl.figure()
ax = fig.add_subplot(111)

# legends
fs = 20
ax.set_xlabel(r'$r\quad[{\rm pc}/h]$',fontsize=fs)
ax.set_ylabel(r'$\rho\quad[M_\odot/{\rm pc}^3]$',fontsize=fs) #TODO: units?
xticks(fontsize=fs)
yticks(fontsize=fs)

jj=0
# DM: read in rad_bin, rho_bin
f=open("prof_dm_"+str(jj+1)+".dat","r")
rad_bin_dm_0=[]; rho_bin_dm_0=[]; n_dm_0=[]
NR_dm_0 = 0
for line in f:
    if(line[0]=="#"):
        continue
    val=line.split()
    rad_bin_dm_0.append(float(val[0]))
    rho_bin_dm_0.append(float(val[3]))
    n_dm_0.append(float(val[5]))
    NR_dm_0 += 1
f.close()

ne_dm_0=[]
for i in xrange(NR_dm_0):
    if(n_dm_0[i-1]==0):
        ne_dm_0.append(0);
    else:
        ne_dm_0.append(rho_bin_dm_0[i-1]/math.sqrt(n_dm_0[i-1]))
        
# plot log-log
ax.loglog(rad_bin_dm_0,rho_bin_dm_0,linewidth="2",c='g')



jj=1
# DM: read in rad_bin, rho_bin
f=open("prof_dm_"+str(jj+1)+".dat","r")
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

ne_dm=[]
for i in xrange(NR_dm):
    if(n_dm[i-1]==0):
        ne_dm.append(0);
    else:
        ne_dm.append(rho_bin_dm[i-1]/math.sqrt(n_dm[i-1]))
        
# plot log-log
ax.loglog(rad_bin_dm,rho_bin_dm,linewidth="2",c='r')


# error bars, assume Poissonian errors
ax.errorbar(rad_bin_dm[0:NR_dm],rho_bin_dm[0:NR_dm],yerr=ne_dm[0:NR_dm],c='r')

# error bars, assume Poissonian errors
ax.errorbar(rad_bin_dm_0[0:NR_dm_0],rho_bin_dm_0[0:NR_dm_0],yerr=ne_dm_0[0:NR_dm_0],c='g')

# vertical line @ resolution limit 4pc given by A. Boley
# (CIC @ max level 14 for z~12, 1Mpc/h box, 256^3 equivalent part.
ax.axvline(x=8.0,c='b')
ax.axvline(x=28.0,c='g')
ax.axvline(x=35.0,c='r')

# restrict display area
ax.set_xbound(3,500)
ax.set_ybound(0.01,100)
legend(['dm 1','dm 2'],'upper right')
show()
pl.savefig("map.dat.png")
