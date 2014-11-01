#!/usr/bin/python
# shrink_spheres.py
# invoke shrink_sphere.py for all halos

import os
import sys
import threading
import initialize as my

os.nice(10)

d = sys.argv[1]+"/"

for i in range(mys.get_nhalo(snap)):
    j = i + 1
    si = str(j)
    cmd = "shrink_sphere.py "+str(snap)
    print cmd
    os.system(cmd)
    # this is not done for the stars
    # we want to center only on dm particles
