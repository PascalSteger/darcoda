#!/usr/bin/python
# TODO

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
if i!=TODO:
    print "Usage: TODO.py TODO"

nproc = 48
os.nice(8)

semaphore = threading.Semaphore(nproc)
def run_command(cmd):
    with semaphore:
        os.system(cmd)

# TODO
threading.Thread(target=run_command, args=("sample command -help", )).start()
