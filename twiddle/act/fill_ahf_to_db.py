#!/usr/bin/env python2

## \file
# after AHF has finished, fill found halo values into MySQL database

# (c) 2015 ETHZ, Pascal Steger, pascal@steger.aero

import os, sys
import lib.mysql as mys
import lib.initialize as my

nsnap = int(sys.argv[1])
d = mys.d(nsnap)

# reading center, rvir, make comoving again
if(not os.path.exists(d+"halos")):
    exit

co, cu = my.sqlstart()

halos = open(d+"halos","r")
hid=0
for line in halos:
    if(line[0]=="#"):
        continue
    hid = hid+1
    val = line.split()
    mys.fill_halo(nsnap, hid, val, co, cu)

    # if less than 100 particles in halo
    # if nvpart (mass in internal units)/npart = fmhires < 5 (arbitrary! but excludes halos with many lowres particles)
    # attention: does not work, as star particles are less massive than DM particles, have lower fmhires like that, so we exclude halos with no stars more often
    # workaround: set cut in fmhires to large value
    # if outside refined zone (0.25 from edge): exclude
    # exclude halo from further analysis

halos.close()
co.commit()
my.sqlstop(co,cu)

#mys.physical_xcm(nsnap,0.719)
#mys.physical_xcm(nsnap,0.702*0.702) #convert all values from AHFstep output to physical units
mys.physical_xcm(nsnap, 1000.0) # convert all values from AHFstep output (physical) to ramses code units (0,1)

#mys.exclude(nsnap,300,100,0.3)       #excluding the outer 30%, and halos with < 300 particles
mys.exclude(nsnap, 30, 0, 0.0) # out to very edge. TODO: check/implement overwrapping
mys.fill_snapshot(nsnap)
