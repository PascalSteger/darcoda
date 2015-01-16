#!/usr/bin/python
'''plot one column of a file against another'''
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
    print "Usage: plot_columns.py file column1x column2y"
    print "Example: plot_columns.py mstar 1 4"
    exit(1)

fi=sys.argv[1]
col1=int(sys.argv[2])
col2=int(sys.argv[3])

# read in file

x=[]; y=[];
fil=open(fi,"r")

for line in fil:
  var=line.split()
  x.append(var[col1])
  y.append(var[col2])

fil.close()

pylab.figure()

pylab.plot(x,y,marker="o",linewidth="2",color="blue")

# visible region
#pylab.xlim([10**1,10**3])
#pylab.ylim([10**-2,10**2])
pylab.grid(True)
pylab.xlabel(r'$col 1$')
pylab.ylabel(r'$col 2$')

pylab.savefig(fi+".png")
