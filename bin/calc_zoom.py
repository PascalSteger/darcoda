#!/usr/bin/python
# calc smooth transition between x,y,z,r1 to x,y,z,r2 in N steps

import os
import sys
import numpy
import math
import pylab
from pylab import *
import matplotlib
matplotlib.use('Agg') # use png output
from matplotlib import rc
rc('text',usetex=True)# use latex in figures
import Image
import fortranfile    # for fortran data import
import initialize as my

# check syntax
i=len(sys.argv)
if i!=7:
    print "Usage: calc_zoom.py x,y,z,r1,r2,N"
    exit(1)
    
x = float(sys.argv[1])
y = float(sys.argv[2])
z = float(sys.argv[3])
r1= float(sys.argv[4])
r2= float(sys.argv[5])
N = int(sys.argv[6])


dr = r1-r2
for i in range(N+1):
    print x,y,z,r2+dr/2+dr/2*(math.cos(i*math.pi/N))

    
