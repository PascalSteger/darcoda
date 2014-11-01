#!/usr/bin/python
# find_min_rhos.py
# invoke find_min_rho.py for all halos

import os
import sys
import threading
import initialize as my

os.nice(10)
nproc=32

d = sys.argv[1]+"/"

os.nice(10)
halo = open(d+"halo","r")
os.system("rm "+d+"halosp")
i = 0
for line in halo:
    i = i + 1
    si = str(i)
    cmd = "find_min_rho.py "+d+"dm/rho_"+si+".dat >>"+d+"halosp"
    print cmd
    os.system(cmd)
    # this is not done for the stars
    # we want to center only on dm particles
