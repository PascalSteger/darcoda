#!/usr/bin/env python2

## \file
# find shared particles between two snapshots
# can be used with twice the same snapshot to get subhalo tree

# (c) 2014 ETHZ, Pascal Steger, pascal@steger.aero

import sys
import lib.initialize as my
import lib.mysql as mys
from numpy import zeros

if(len(sys.argv)!=3):
    print("usage: merger_tree.py snap1 snap2")
    print("assumption: output_00snap1,2 contain particles_dm in AHF format")
    exit(1)

snap1 = int(sys.argv[1]); snap2 = int(sys.argv[2])

# read in halos repeatedly
fname=mys.d(snap1)+"particles_dm"
print('file = ', fname)
in1 = my.open_file(fname,"r")
alm1 = int(in1.readline())
size1 = []; split1 = []
for i in range(alm1):
    size = int(in1.readline())
    size1.append(size)
    part = []
    for c in range(size):
        part.append(int(in1.readline()))
    split1.append(part)
in1.close()

in2 = my.open_file(mys.d(snap2)+"particles_dm","r")
alm2 = int(in2.readline())
size2 = []; split2 = []
for i in range(alm2):
    size = int(in2.readline())
    size2.append(size)
    part = []
    for c in range(size):
        part.append(int(in2.readline()))
    split2.append(part)
in2.close()

print("read in successful")

def count_shared(m1,m2):
    cc = 0
#    print("count")
    for i in range(len(m1)):
        if m1[i] in m2:
            cc += 1
    return cc
                
print("find particles shared by halo_in_1, halo_in_2")
# create sparse matrix of particles shared by each halo with each other
M = zeros([alm1,alm2],int)
for i in range(alm1):
    for j in range(alm2):
        M[i,j]=count_shared(split1[i],split2[j])
        if(M[i,j]>0):
            print(i,j,M[i,j])

print("search progenitor (most particles shared)")
# find progenitor (most particles shared)
for i in range(alm1):
    maxi = 0; maxj = i
    for j in range(alm2):
        if(M[i,j]>maxi):
            maxi = M[i,j]
            maxj = j
    my.sql("update halo set proghid="+str(maxj)+" where snap="+str(snap1)
	   +" and hid="+str(i))
    print(i, maxj, maxi)
