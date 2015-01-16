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
pos = open("mt/pos","r")
i = 0
for line in pos:
    i = i + 1
    #threading.Thread(target=run_command,args=("shrink_sphere.py vid/part_"+str(i).zfill(5), )).start()
    os.system("shrink_sphere.py vid/part_"+str(i).zfill(5))
    # we want to center on dm particles only
