#!/usr/bin/python
# plot all gas in simulation box

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

# check syntax
i=len(sys.argv)
if i!=1:
    print "Usage: no arguments"

nproc = 48
os.nice(1)

import initialize as my
calc=True
show=True
run=True
lma=str(9)
for i in range(80):
    stri=str(i+1).zfill(5)
    folder = "output_"+stri
    # plot gas in whole box
    bndry=" -xmi 0.0 -xma 1.0 -ymi 0.0 -yma 1.0 -zmi 0.0 -zma 1.0 "
    cmd1 = "amr2map -typ 1 -lma "+lma+" -inp "+folder
    cmd1 = cmd1 + " -out mov/gas/box/box_"+stri
    cmd1 = cmd1 + ".dat -dir z "+bndry
    cmc = " && "
    cmd2 = "map2img.py -l --colormap=hot "
    cmd2 = cmd2 + "mov/gas/box/box_"+stri+".dat "
    cmd2 = cmd2 + "-o mov/gas/box/box_"+stri+".png"
    if(calc):
        cmd = cmd1 + cmc + cmd2
    else:
        cmd = cmd2
    if(show):print cmd
    if(run): my.thread(cmd)

    # plot gas in slice
    bndry=" -xmi 0.0 -xma 1.0 -ymi 0.0 -yma 1.0 -zmi 0.4 -zma 0.6 "
    cmd1 = "amr2map -typ 1 -lma "+lma+" -inp "+folder
    cmd1 = cmd1 + " -out mov/gas/slice/slice_"+stri
    cmd1 = cmd1 + ".dat -dir z "+bndry
    cmc = " && "
    cmd2 = "map2img.py -l --colormap=hot "
    cmd2 = cmd2 + "mov/gas/slice/slice_"+stri+".dat "
    cmd2 = cmd2 + "-o mov/gas/slice/slice_"+stri+".png"
    if(calc):
        cmd = cmd1 + cmc + cmd2
    else:
        cmd = cmd2
    if(show):print cmd
    if(run): my.thread(cmd)
