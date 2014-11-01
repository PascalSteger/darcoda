#!/usr/bin/python
'''plot profile'''

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
    print "Usage: plot_prof.py infile outfile"
    print "Example: plot_prof.py prof.dat prof.png"

infile = sys.argv[1]
outfile = sys.argv[2]

# read in rad_bin, rho_bin
f=open(infile,"r")
rad_bin=[]; rho_bin=[]
NR = 0
for line in f:
    if(line[0]=="#"):
        continue
    val=line.split()
    rad_bin.append(float(val[1])*10**6) # conversion to pc/h
    rho_bin.append(float(val[3]))
    NR += 1
    
# plot log-log
pylab.figure()
# TODO: loglog plot command
#pylab.semilogy(rad_bin[0:NR],rho_bin[0:NR])
pylab.loglog(rad_bin[0:NR],rho_bin[0:NR],linewidth="2")

# vertical line @ resolution limit 4pc given by A. Boley
# (CIC @ max level 14 for z~12, 1Mpc/h box, 256^3 equivalent part.
#pylab.axvline(x=0.004)

# error bars, assume Poissonian errors
# TODO

# visible region
pylab.xlim([10**1,10**3])
pylab.ylim([10**-2,10**2])
pylab.grid(True)
pylab.xlabel(r'$r [{\rm pc}/h]$')
pylab.ylabel(r'$\rho [h^2 M_\odot/{\rm pc}^3]$')

pylab.savefig(outfile)
