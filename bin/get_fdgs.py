#!/usr/bin/python

import sys

i = len(sys.argv)
if i!=4:
    print "usage: get_fdgs.py dm.dat star.dat mtot"
    exit(1)

mtot= float(sys.argv[3])/0.719

fdm = open(sys.argv[1])
mdm = 0.0
for line in fdm:
    val = line.split()
    mdm += float(val[0])
fdm.close()

fs  = open(sys.argv[2])
ms  = 0.0
for line in fs:
    val = line.split()
    ms  += float(val[0])
fs.close()

fdm = mdm/mtot
fs  = ms/mtot
fg  = (mtot-mdm-ms)/mtot

print "f_dm = ",fdm
print "f_s  = ",fs
print "f_g  = ",fg
