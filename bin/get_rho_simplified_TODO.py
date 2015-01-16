#!/usr/bin/python
# input:  star particles
# output: Plummer profile, concentration, stellar cluster

import os
import sys
#import array
from numpy import *
import scipy
import math
import matplotlib
matplotlib.use('Agg') # use png output
from matplotlib import rc
rc('text',usetex=True)# use latex in figures

import Image
import fortranfile    # for fortran data import
import initialize as my

# check syntax
i=len(sys.argv)
if i!=5:
    print "Usage: fit_plummer.py prof_stars.dat xc yc zc"
    exit(1)

filename=sys.argv[1]
f = my.open_file(filename,"r")
c = loadtxt(f)
f.close()

xc = float(sys.argv[2])
yc = float(sys.argv[3])
zc = float(sys.argv[4])

x=array(c[:,0])-xc
y=array(c[:,1])-yc
z=array(c[:,2])-zc
m=array(c[:,3])
rho=array(c[:,4])

r=sqrt(x*x+y*y+z*z)
order=r.argsort()
r=r[order]; rho=rho[order]
print r,rho
