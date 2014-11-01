#!/usr/bin/python
# compile set of png images into a movie
# top left:     dm density
# top right:    star density
# middle left:  gas density
# middle right: gas temperature
# bottom left:  metallicity
# bottom right: empty

import os
import threading

from sys import *
from numpy import *
from math import *
#from numarray import *

from pylab import *
import matplotlib
matplotlib.use('Agg') # use png output

import Image
import fortranfile    # for fortran data import

# check syntax
i=len(sys.argv)
if i!=1:
    print "Usage: vid_compile.py"
    sys.exit(1)
    
nproc = 48
os.nice(8)

semaphore = threading.Semaphore(nproc)
def run_command(cmd):
    with semaphore:
        os.system(cmd)

# TODO:step through all snapshots, execute in parallel

# create image array
ov = zeros((3,2))
#im=Image.open("mov/dmdens/map_00154.png")
#ivals=array(im.getdata())
ivals = imread("mov/dmdens/map_00154.png")
print ivals
ov[0][0]=ivals
ov[0,1]=ivals
ov[1,0]=ivals
ov[1,1]=ivals

print ov
