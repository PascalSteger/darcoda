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
if i!=1:
    print "Usage: plot_profs_dln.py"
    print "Example: plot_profs_dln.py"
    exit()

fig = figure()
ax = fig.add_subplot(111)

ii=1
# DM: read in rad_bin, rho_bin
f=open("prof_dm_"+str(ii)+".dat","r")
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

# errors
ne_dm=[]
for i in xrange(NR_dm):
    if(n_dm[i-1]==0):
        ne_dm.append(0);
    else:
        ne_dm.append(rho_bin_dm[i-1]/sqrt(n_dm[i-1]))
            
# dln\rho/dln r
dlrho = []; dlr = []
for i in xrange(NR_dm-1):
    dlrho.append(log(rho_bin_dm[i])-log(rho_bin_dm[i+1]))
    dlr.append(log(rad_bin_dm[i])-log(rad_bin_dm[i+1]))
    
dd_dm = array(dlrho)/array(dlr)
    
# plot semi log in x direction
ax.semilogx(rad_bin_dm[:NR_dm-1],dd_dm[:NR_dm-1],linewidth="2",c='g')


ii=2
# DM: read in rad_bin, rho_bin
f=open("prof_dm_"+str(ii)+".dat","r")
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

# errors
ne_dm=[]
for i in xrange(NR_dm):
    if(n_dm[i-1]==0):
        ne_dm.append(0);
    else:
        ne_dm.append(rho_bin_dm[i-1]/sqrt(n_dm[i-1]))
            
# dln\rho/dln r
dlrho = []; dlr = []
for i in xrange(NR_dm-1):
    dlrho.append(log(rho_bin_dm[i])-log(rho_bin_dm[i+1]))
    dlr.append(log(rad_bin_dm[i])-log(rad_bin_dm[i+1]))
    
dd_dm = array(dlrho)/array(dlr)
    
# plot semi log in x direction
ax.semilogx(rad_bin_dm[:NR_dm-1],dd_dm[:NR_dm-1],linewidth="2",c='r')


# vertical line @ resolution limit 4pc given by A. Boley
# (CIC @ max level 14 for z~12, 1Mpc/h box, 256^3 equivalent part.
ax.axvline(x=8.0,c='b')
ax.axvline(x=28.0,c='g')
ax.axvline(x=35.0,c='r')

# restrict display area
ax.set_xbound(3,500)
ax.set_ybound(-3,2)

# legends
fs=20
legend(['dm 1','dm 2'],'upper right')
ax.set_xlabel(r'$r\quad[{\rm pc}/h]$',fontsize=fs)
ax.set_ylabel(r'$d \ln \rho/d\ln r$',fontsize=fs) #TODO: units?
xticks(fontsize=fs)
yticks(fontsize=fs)

savefig("map.dat.png")
