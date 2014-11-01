#!/usr/bin/python
'''plot Mstar vs rshalf_star'''
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
    print "Usage: plot_Mstar_rhalf.py snap outfile"
    print "Example: plot_Mstar_rhalf.py 270 mr2_00270.png"
    exit(1)

snap=int(sys.argv[1])
outfile=sys.argv[2]

# read in Mstar, rhalf
mstar= array(mys.get_M_star(snap))
mdm  = array(mys.get_M_dm(snap))
issub= array(mys.get_is_sub(snap))
rhalf= array(mys.get_rhalf_star(snap))*1000
mtot=mdm+mstar
cc=mstar/mtot

print cc

rh1=[];ms1=[]; rh2=[]; ms2=[]

pylab.figure()

for i in range(len(cc)):
    
    if(cc[i]>0.5):
        co="red"
    else:
        co="blue"

    if(issub[i]>0):
        marker="^"
    else:
        marker="o"
    # plot log-log
    pylab.loglog(rhalf[i],mstar[i],marker,linewidth="2",color=co)

# visible region
#pylab.xlim([10**1,10**3])
#pylab.ylim([10**-2,10**2])
pylab.grid(True)
pylab.xlabel(r'$r_{1/2} [{\rm kpc}/h]$')
pylab.ylabel(r'$M_* [h^2 M_\odot/{\rm pc}^3]$')

pylab.savefig(outfile)
