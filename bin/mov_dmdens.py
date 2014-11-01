#!/usr/bin/python
# follow pos_sss, make dmdens movie

import os
import threading

import sys
import numpy
import math

import matplotlib
matplotlib.use('Agg') # use png output
from matplotlib import rc
rc('text',usetex=True)# use latex in figures

import pylab
from pylab import *


import Image
import fortranfile    # for fortran data import

# check syntax
i=len(sys.argv)
#if i!=TODO:
#    print "Usage: TODO.py TODO"

nproc = 48
os.nice(8)

semaphore = threading.Semaphore(nproc)
def run_command(cmd):
    with semaphore:
        os.system(cmd)

f = open('mt/cen','r')
i = 0
for line in f:
    i=i+1
    print i
#    if(i>10):
#        sys.exit(0)
    si = str(i).zfill(5)
    val = line.split()
    x=float(val[0]); y=float(val[1]); z=float(val[2]); r=float(val[3])
#    r=1
    xmi = str(x-r/2); xma=str(x+r/2)
    ymi = str(y-r/2); yma=str(y+r/2)
    zmi = str(z-r/2); zma=str(z+r/2)
#    threading.Thread(target=run_command, args=("part2map -inp output_"+si+" -out vid/dmdens/map_d_"+si+".dat -dir z -xmi "+xmi+" -ymi "+ymi+" -zmi "+zmi+" -xma "+xma+" -yma "+yma+" -zma "+zma+" -nx 512 -ny 512", )).start()

    cmd = "get_sphere_dm_8 -inp output_"+si+" -xc "+str(x)+" -yc "+str(y)+" -zc "+str(z)+" -r "+str(r)+" > vid/dmdens/part_"+si+".dat"
    print cmd
    #os.system(cmd)
    cmd = "get_prof.py vid/dmdens/part_"+si+".dat "+str(i)+" "+str(x)+" "+str(y)+" "+str(z)+" "+str(r)+" > vid/dmdens/prof_"+si+".dat"
    print cmd
    #os.system(cmd)

for i in range(270):
    si = str(i+1).zfill(5)
    #threading.Thread(target=run_command, args=("map2img.py -o vid/dmdens/map_d_"+si+".png vid/dmdens/map_d_"+si+".dat", )).start()
    os.system("plot_prof.py vid/dmdens/prof_"+si+".dat && mv map.dat.png vid/dmdens/prof_"+si+".png")
