#!/usr/bin/python
# shrink_spheres.py
# invoke shrink_sphere.py for all halos

import os
import sys
import threading

nproc = 24

semaphore = threading.Semaphore(nproc)
def run_command(cmd):
    with semaphore:
        os.system(cmd)

os.nice(10)
pos = open("mt/cen","r")
i = 0; x=[]; y=[]; z=[]; r=[]
for line in pos:
    val = line.split()
    x.append(float(val[0]))
    y.append(float(val[1]))
    z.append(float(val[2]))
    r.append(float(val[3]))

for i in range(len(x)):
    if(x[-i-1]==0):
        # no valid position found, approximate by last known position
        x[-i-1]=x[-i]
        y[-i-1]=y[-i]
        z[-i-1]=z[-i]
        r[-i-1]=r[-i]

for j in range(i):
    print x[j],y[j],z[j],r[j]
