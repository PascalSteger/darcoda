#!/usr/bin/python
# run AHFstep from python, test script, will be included in master.py 2

import os
import threading

import sys
import numpy
import math

import Image
import fortranfile    # for fortran data import

# check syntax
i=len(sys.argv)
if i!=1:
    print "Usage: runahf.py"

nproc = 48
os.nice(8)

semaphore = threading.Semaphore(nproc)
def run_command(cmd):
    with semaphore:
        os.system(cmd)

simdir = "/scratch/psteger/sim_aboley/"
nc=260
num5 = str(nc).zfill(5)
d  = simdir + "/output_"+num5+"/"
os.system("cat "+d+"a.par")
os.system("cd "+d+" && AHFstep a.par")
