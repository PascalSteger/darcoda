#!/usr/bin/python

import sys
import os
import math
import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
rc('text',usetex=True)
import pylab as pl
from pylab import *


for i in range(270):
    folder = "/scratch/psteger/sim_aboley/output_" + str(i+1).zfill(5) + "/"
    if not os.path.exists(folder):
            os.makedirs(folder)
    apar   = folder + "a.par"
    f=open(apar,"w")
    f.write(folder + "r2g/r2g. 61 1\n")
    f.write("ao\n")
    f.write("16\n")
    f.write("4\n")
    f.write("4\n")
    f.write("0\n")
    f.write("0\n")
    f.write("0\n")
    f.write("0\n")
    f.close()

    os.system("cd "+folder+" && AHFstep a.par || true")
