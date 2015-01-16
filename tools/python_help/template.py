#!/usr/bin/python
# TODO

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
if i!=TODO:
    print "Usage: TODO.py TODO"
    exit(1)
