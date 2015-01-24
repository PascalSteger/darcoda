#!/usr/bin/env python2

## \file
# find shared particles between two snapshots
# can be used with twice the same snapshot to get subhalo tree

# (c) 2014 ETHZ, Pascal Steger, pascal@steger.aero

import sys
from numpy import zeros
import lib.initialize as my
import lib.mysql as mys

if(len(sys.argv)!=2):
    print("usage: substructure.py snap")
    print("assumption: output_snap contains particles_dm in AHF format")
    print("wrong number of arguments")
    exit(1)

snap = int(sys.argv[1])
# read in halos repeatedly
inf = my.open_file(mys.d(snap)+"particles_dm","r")
alm = int(inf.readline())
size = []; split = []
for i in range(alm):
    siz = int(inf.readline())
    size.append(siz)
    part = []
    for c in range(siz):
        part.append(int(inf.readline()))
    split.append(part)
inf.close()
print("read in successful")


def count_shared(m1,m2):
    cc = 0
    for i in range(len(m1)):
        if m1[i] in m2:
            cc = cc + 1
    return cc
## \fn count_shared(m1,m2)
# count common elements in arrays m1 and m2
# @param m1 array1
# @param m2 array2
                
print("find particles shared by halo_in_1, halo_in_2")
# create sparse matrix of particles shared by each halo with each other
M = zeros([alm,alm],int)
for i in range(alm):
    for j in range(alm):
        M[i,j]=count_shared(split[i],split[j])
        if(M[i,j]!=0):
            print(i,j,M[i,j])

print("search host (most particles shared)")
for i in range(alm):
    maxi = 0; maxj = i
    for j in range(alm):
        if(M[i,j]>maxi):
            maxi = M[i,j]
            maxj = j
    my.sql("update halo set hosthid="+str(maxj)+" where snap="+str(snap)
	   +" and hid="+str(i))
    print(i, maxj, maxi)
