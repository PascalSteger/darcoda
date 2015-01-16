#!/usr/bin/python

import os
import sys
import ../lib/mysql as mys
import fileinput

print "8"
nsnap=int(sys.argv[1])

d=mys.d(nsnap)
print d
exit(0)
# reading center, rvir, make comoving again

if(not os.path.exists(d+"halos")):
    exit
    
halos = open(d+"halos","r")
hid=0
for line in halos:
    if(line[0]=="#"):
        continue
    hid=hid+1
    val = line.split()
    mys.fill_halo(nsnap,hid,val)
    # if less than 200 particles in halo
    # if less than 90% of mass given by high res particles
    # (more than 10% contamination from low res particles)
    # exclude halo from further analysis
    # if outside refined zone (0.25 from edge): exclude
halos.close()
mys.physical_xcm(nsnap,0.719)

mys.exclude(nsnap,100,0.5,0.0)
    
mys.fill_snapshot(nsnap)
