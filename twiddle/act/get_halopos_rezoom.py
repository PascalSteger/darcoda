#!/usr/bin/env python2

## \file
# read halo props from DB, find three most massive halos with no other
# halo of mvir/2 inside 5rvir

# (c) 2014 ETHZ, Pascal Steger, pascal@steger.aero

import lib.mysql as mys
from math import sqrt


def dist_mod(x1,y1,z1,x2,y2,z2):
    dx = abs(x1-x2)
    if(dx>0.5):
        dx = 1-dx
    dy = abs(y1-y2)
    if(dy>0.5):
        dy = 1-dy
    dz = abs(z1-z2)
    if(dz>0.5):
        dz = 1-dz
    return sqrt(dx**2.+dy**2.+dz**2.)
## \fn dist_mod(x1,y1,z1,x2,y2,z2)
# calc_dist modulo Lbox
# @param x1 pos1
# @param y1 pos1
# @param z1 pos1
# @param x2 pos2
# @param y2 pos2
# @param z2 pos2



def snap(snapstop):
    print("find 3 most massive halos")
    # read in pos
    x,y,z,m,r=mys.getxyzmr(snapstop, 1)
    print('mass ', m)
    print('radius ',r)
    
    count=0
    for i in range(len(m)):
        print("halo ", i, m[i], r[i])
        hassub = False
        for j in range(i+1,len(m)):
            # if too low mass halo, continue
            if(m[j]<0.5*m[i]):
                continue
            # have massive halo, possibly subhalo
            if(dist_mod(x[i],y[i],z[i],x[j],y[j],z[j])<10.*r[i]):
                hassub = True
                print(" lies near ", j, m[j], r[j])
                break
        if(count>=3):
            break
        if(not hassub):
            count+=1
            print(x[i],y[i],z[i],r[i],m[i])
## \fn snap(snapstop)
# output three most massive halos
# @param snapstop snapshot ID

            
if __name__ == "__main__":
    snap(10)
