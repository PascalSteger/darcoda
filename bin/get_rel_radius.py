#!/usr/bin/python
import sys
import math
from math import *
import frange as fr

NBIN = 1e4

print(" constants [in SI]")
h     = 0.719
pc    = 3.08568025e16
Mpch  = 1e6*pc/h
yr    = 3.1556926e7
G     = 6.67300e-11
Msunh = 1.98892e30/h

print("read in")
fn   = sys.argv[1]
xc   = float(sys.argv[2])*Mpch
yc   = float(sys.argv[3])*Mpch
zc   = float(sys.argv[4])*Mpch
rvir = float(sys.argv[5])*Mpch

tsim = 0.466e9*yr    # z=10
dr   = rvir / NBIN
bmax = rvir
bmin = 8e-6*Mpch

print("read in all particles of the halo")
f = file(fn,"r");
x=[];y=[];z=[];m=[];r=[]; Ntot = 0
for line in f:
    val = line.split()
    m.append(float(val[0])*Msunh)
    x.append(float(val[1])*Mpch);
    y.append(float(val[2])*Mpch)
    z.append(float(val[3])*Mpch)
    r.append(sqrt((x[Ntot]-xc)**2+(y[Ntot]-yc)**2+(z[Ntot]-zc)**2))
    Ntot = Ntot + 1
    
print("increase radius until relaxation time is simulation time")
for rtmp in fr.frange(0.0,rvir,dr):
    #print "rtmp       = ",rtmp
    #determine N,M
    M = 0.0; N = 0;
    for i in range(0,Ntot):
        if(r[i]<rtmp):
            M = M + m[i];
            N = N + 1;
    if(M<=0.0):
        continue
        
    CL     = log(rtmp/bmin)
    trelax = N*rtmp**1.5/(8*sqrt(G)*sqrt(M)*CL)

    if(trelax>tsim):
        print "r_relax      = ",rtmp/Mpch*1e6,"pc/h"
        exit(0)

print "r_relax > rvir"
