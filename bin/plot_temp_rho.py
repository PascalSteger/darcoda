#!/usr/bin/python
# plot phase diagram (cut) for a snapshot
# needs temperature and density for all cells of a given level
# as produced by onsim/get_temp


import initialize as my
import matplotlib
import numpy as np
import random
import sys
matplotlib.use('Agg')
import pylab
from matplotlib import rc
rc('text',usetex=True)

i=len(sys.argv)
if i!=3:
    print "Usage: plot_temp_rho.py infile outfile"
    print "Example: plot_temp_rho.py sample.dat sample.png"
    exit(1)


infile=sys.argv[1]
outfile=sys.argv[2]

f=open(infile,"r")
tmp_bin=[]; rho_bin=[]
NR = 0
for line in f:
    if(random.random()>0.1):
        continue
    print NR
    val=line.split()
    tmp_bin.append(np.log10(float(val[0])))
    rho_bin.append(np.log10(float(val[1])))
    NR += 1
    
pylab.figure()
#pylab.scatter(rho_bin[0:NR],tmp_bin[0:NR],'.',markersize=0.5)
pylab.scatter(rho_bin[0:NR],tmp_bin[0:NR],marker='o',s=3,alpha=0.5,edgecolors='none')
pylab.grid(True)
pylab.xlabel(r'$\log \rho$')
pylab.ylabel(r'$\log T/\mu\,[K]$')
#pylab.axis([-1,6,1,4])
pylab.savefig(outfile)
