#!/usr/bin/python
#plot mstar(z) from mysql DB
#output d N/dz

import mys
import sys
import math
import matplotlib
matplotlib.use('Agg')
import pylab
from pylab import *
from matplotlib import rc
rc('text',usetex=True)# use latex in figures

# check
i=len(sys.argv)
if i!=3:
    print "Usage: plot_sfr_halo.py snap_start snap_stop"
    print "Example: plot_sfr_halo.py 200 270"
    exit(1)

snapstart = int(sys.argv[1])
snapstop  = int(sys.argv[2])
dnz=[]; z=[]
i=0;

for snap in range(snapstart,snapstop):
    if(not mys.exists_snap(snap)):
        continue
    i=i+1
    mstar=mys.get_M_star(snap)[0]
    zsnap=mys.get_z(snap)
    print mstar,zsnap
    if(i==1):
        zold=zsnap
        nold=mstar
        continue
    znew=zsnap
    nnew=mstar
    z.append(znew)
    dnz.append((nnew-nold)/(zold-znew)) #looks clearer without log, for whole simulation
    zold=znew
    nold=nold

#sort by z so lines between points in plot are useful
z=array(z); dnz=array(dnz)
order=z.argsort()
z=z[order]; dnz=dnz[order]

# plot
pylab.figure()

pylab.plot(z,dnz,marker="o",linewidth="2",color="blue")

pylab.grid(True)
pylab.xlabel(r'$z$')
pylab.ylabel(r'$dn/dz$')

pylab.savefig("halo1_sfr.png")
