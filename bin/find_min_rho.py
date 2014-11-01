#!/usr/bin/python
# determine center of mass by searching the smallest sph density

import sys
import os
import numpy as np
from numpy import *
import math

i=len(sys.argv)
if(i!=2):
    print "find_min_rho.py infile"
    exit(1)


in_file=sys.argv[1]
try:
    f = open(in_file,"r")
except(IOError),e:
    exit(1)

x=[]; y=[]; z=[]; rho=[]

# missing snapshot
if(os.path.getsize(in_file)==51):
    print "0 0 0 0"
    exit(1)

for line in f:
    val = line.split()
    x.append(float(val[0])); y.append(float(val[1])); z.append(float(val[2]))
    rho.append(float(val[4]))

f.close()


x = array(x); y = array(y); z = array(z)
rho = array(rho)

xmin = min(x); xmax = max(x); r = 0.9*(xmax-xmin)/2

minrho = -1e99
minx = 0.5; miny = 0.5; minz = 0.5
for i in range(len(x)):
    if(rho[i]>minrho):
        minx = x[i]; miny = y[i]; minz = z[i]
        minrho = rho[i]

print str(minx)+" "+str(miny)+" "+str(minz)+" "+str(r)
