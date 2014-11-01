#!/usr/bin/python
# given a simulation and a file with halo centers, do gen_prof

import sys
import os

i=len(sys.argv)
if i != 2:
    print "usage: ./gen_profs.py centerfile"
    print "example: ./gen_profs.py mt/cen"
    sys.exit()

cen = sys.argv[1]

f=open(cen,"r")
i=0
for line in f:
    i = i+1
    val=line.split()
    x = val[2]; y = val[3]; z = val[4]; r = val[5]
    cmd = "gen_prof.py "+str(i)+" "+x+" "+y+" "+z+" "+r;
    print cmd;
    thread(cmd)
    #os.system(cmd)
