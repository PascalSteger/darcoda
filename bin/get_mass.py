#!/usr/bin/python
# gets cumulative mass of all particle in dm/star spheres

import os
import threading

import sys
import numpy
import math

from matplotlib import rc
rc('text',usetex=True)# use latex in figures

import Image
import fortranfile    # for fortran data import

# check syntax
i=len(sys.argv)
if i!=2:
    print "Usage: get_mass.py dm.dat"
    exit(1)

filename = sys.argv[1]
if(os.path.getsize(filename) == 0):
    print "0.0"
    exit(0)

nproc = 48
os.nice(8)

semaphore = threading.Semaphore(nproc)
def run_command(cmd):
    with semaphore:
        os.system(cmd)

# TODO
#threading.Thread(target=run_command, args=("sample command -help", )).start()

f = open(filename,"r")
mtot = 0.0
for line in f:
    val = line.split()
    mtot += float(val[0])

print mtot
