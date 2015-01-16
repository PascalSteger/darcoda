#!/usr/bin/python
#read mstar with lines: z, nstar_tot, mstar_tot, nstar_actual, mstar_actual
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
if i!=4:
    print "Usage: calc_sfr.py file column1x column2y"
    print "Example: calc_sfr.py mstar 1 4"
    exit(1)

inf = open(sys.argv[1],"r")
col1 = int(sys.argv[2])
col2 = int(sys.argv[3])

dnz=[]; z=[]
i=0;
for line in inf:
    i=i+1
    val=line.split()
    if(i==1):
        zold=float(val[col1])
        nold=float(val[col2])
        continue
    znew=float(val[col1])
    nnew=float(val[col2])
    z.append(znew)
    dnz.append((nnew-nold)/(zold-znew)) #looks clearer without log, for whole simulation
    zold=znew
    nold=nold
inf.close()

# TODO: sort by z for nice output
z=array(z); dnz=array(dnz)
order=z.argsort()
z=z[order]; dnz=dnz[order]

# plot
pylab.figure()

pylab.plot(z,dnz,marker="o",linewidth="2",color="blue")

pylab.grid(True)
pylab.xlabel(r'$z$')
pylab.ylabel(r'$dn/dz$')

pylab.savefig(sys.argv[1]+"_sfr.png")
