#!/usr/bin/python
# shrink_spheres.py
# invoke shrink_sphere.py for all halos on stellar part

import os
import sys
import threading
import initialize as my

d = sys.argv[1]+"/"

os.nice(10)
halo = open(d+"halo","r")
os.system("rm "+d+"halosp")
i = 0
for line in halo:
    i = i + 1
    si = str(i)
    cmd = "shrink_sphere.py "+d+"stars/stars_"+si+".dat >>"+d+"halosp"
    print cmd
    my.run(cmd)
    # this is not done for the stars
    # we want to center only on dm particles
