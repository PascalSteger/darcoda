#!/usr/bin/python
# read redshift from a info file

import os
import sys
import numpy
import math
import Image
import fortranfile    # for fortran data import
import initialize as my

# check syntax
i=len(sys.argv)
if i!=2:
    print "Usage: read_z.py snap"
    exit(1)

snap = int(sys.argv[1])
num5=str(snap).zfill(5)
filename="output_"+num5+"/info_"+num5+".txt"
if(not os.path.exists(filename)): exit(1)
file = open(filename)
i=0; z=100.0
for line in file:
    i=i+1
    if(i==10):
        val = line.split()
        a = float(val[2])
        z = 1/a-1.0

print z
