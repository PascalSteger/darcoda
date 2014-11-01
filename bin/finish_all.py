#!/usr/bin/python
# finish all box / slice pictures inside mov

import os
import threading

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
if i!=3:
    print "Usage: finish_all.py snapstart snapstop"
    exit(1)

os.nice(1)

nstart = int(sys.argv[1])
nstop = int(sys.argv[2])

show=True
run=True
radius=500 #[kpc/h]
simdir = "/scratch/psteger/sim/"
for ncounter in range(nstop-nstart+1):
    nc = ncounter + nstart
    num5 = str(nc).zfill(5)
    d  = simdir + "/output_"+num5+"/"

    if(not my.file_exists(d+"info_"+num5+".txt")):
        continue
    os.system("head "+d+"/info_"+num5+".txt | tail -n1 | cut -d'=' -f2>"+d+"aexp")
    af = open(d+"aexp","r")
    a  = af.readline()
    print float(a)

    cmd = "finish_gas_box.py "+simdir+"mov/gas/gas_boxall_"+num5+".png "+simdir+"mov/gas/fgas_boxall_"+num5+".png "+str(nc)+" "+str(float(a))+" "+str(radius)
    if(show): print cmd
    if(run):  my.thread(cmd)
