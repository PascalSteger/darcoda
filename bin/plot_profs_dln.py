#!/usr/bin/python
'''plot profiles (DM/stars together)'''

import sys
import math
from math import *
import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
rc('text',usetex=True)
import pylab
from pylab import *

# check
i=len(sys.argv)
if i!=4:
    print "Usage: plot_profs_dln.py prof_dm prof_stars prof_dmonly"
    print "Example: plot_profs_dln.py prof_dm.dat prof_stars.dat prof_dmonly.dat"
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

# DM only: read in rad_bin, rho_bin
f=open(sys.argv[3],"r")
rad_bin_dmonly=[]; rho_bin_dmonly=[]; n_dmonly=[]
NR_dmonly = 0
for line in f:
    if(line[0]=="#"):
        continue
    val=line.split()
    rad_bin_dmonly.append(float(val[0]))
    rho_bin_dmonly.append(float(val[3]))
    n_dmonly.append(float(val[5]))
    NR_dmonly += 1
f.close()


# errors
ne_dm=[]
for i in xrange(NR_dm):
    if(n_dm[i-1]==0):
        ne_dm.append(0);
    else:
        ne_dm.append(rho_bin_dm[i-1]/sqrt(n_dm[i-1]))

ne_dmonly=[]
for i in xrange(NR_dmonly):
    if(n_dmonly[i-1]==0):
        ne_dmonly.append(0);
    else:
        ne_dmonly.append(rho_bin_dmonly[i-1]/sqrt(n_dmonly[i-1]))

ne_stars=[]
for i in xrange(NR_stars):
    if(n_stars[i-1]==0):
        ne_stars.append(0);
    else:
        ne_stars.append(rho_bin_stars[i-1]/sqrt(n_stars[i-1]))

# dln\rho/dln r
dlrho = []; dlr = []
for i in xrange(NR_dm-1):
    dlrho.append(log(rho_bin_dm[i])-log(rho_bin_dm[i+1]))
    dlr.append(log(rad_bin_dm[i])-log(rad_bin_dm[i+1]))
dlorho = []; dlro = []
for i in xrange(NR_dmonly-1):
    dlorho.append(log(rho_bin_dmonly[i])-log(rho_bin_dmonly[i+1]))
    dlro.append(log(rad_bin_dmonly[i])-log(rad_bin_dmonly[i+1]))

dd_dm = array(dlrho)/array(dlr)
dd_dmonly = array(dlorho)/array(dlro)
    
# plot log-log
fig = figure()
ax = fig.add_subplot(111)
ax.semilogx(rad_bin_dm[0:NR_dm-1],dd_dm,linewidth="2",c="b")
ax.semilogx(rad_bin_dmonly[0:NR_dmonly-1],dd_dmonly,linewidth="2",c="g")
#ax.loglog(rad_bin_dm,rho_bin_dm,linewidth="2",c="b")

# vertical line @ resolution limit 4pc given by A. Boley
# (CIC @ max level 14 for z~12, 1Mpc/h box, 256^3 equivalent part.
ax.axvline(x=4.0,c="g")
ax.axvline(x=8.0,c="b")
ax.axvline(x=28.0,c="r")

# error bars, assume Poissonian errors
#ax.errorbar(rad_bin_dm[0:NR_dm],rho_bin_dm[0:NR_dm],yerr=ne_dm[0:NR_dm],c="b")

# restrict display area
ax.set_xbound(3,300)
ax.set_ybound(-3,1)

# legends
fs = 20
ax.legend(('dm (hydro sim)','dm (dm only sim)'),#'stars'),\
           'upper right',shadow=True)
ax.set_xlabel(r'$r [{\rm pc}/h]$',fontsize=fs)
ax.set_ylabel(r'$d \ln \rho/d\ln r$',fontsize=fs) #TODO: units?
xticks(fontsize=fs)
yticks(fontsize=fs)

savefig("map.dat.png")
