#!/usr/bin/env python3

## \file
# plot scatter plots from halos file

import sys
i=len(sys.argv)
if i<4:
    print("Usage: scatter.py infile outfile columnx(start from 1) columny [log]")
    print("Example: scatter.py halos scatter.png 12 4 1")
    exit(1)
    
infile  = sys.argv[1]; outfile = sys.argv[2]
colx     = [int(sys.argv[3])-1]
coly     = [int(sys.argv[4])-1]
boolog  = False
if(i>5):
    boolog  = int(sys.argv[5])>0

from matplotlib import pyplot as PLT
fig = PLT.figure()
ax1 = fig.add_subplot(111)

import numpy as NP
with open(infile) as f:
    vx = NP.loadtxt(f, dtype='float', comments="#", skiprows=0, usecols=colx)
    vx = NP.ravel(vx)   # 'flatten' v, cumulative sum
with open(infile) as g:
    vy = NP.loadtxt(g, dtype='float', comments="#", skiprows=0, usecols=coly)
    vy = NP.ravel(vy)   # 'flatten' v, cumulative sum
with open(infile) as h:
    #38 for fmhires
    fhi = NP.loadtxt(h, dtype='float', comments="#", skiprows=0, usecols=[2-1])
    fhi = NP.ravel(fhi)   # 'flatten' v, cumulative sum

thresh = 0.0
vxhi= vx[fhi>thresh];  vyhi= vy[fhi>thresh]
vxlo= vx[fhi<=thresh]; vylo= vy[fhi<=thresh]

if(boolog):
    vxhi = NP.log10(vxhi);    vyhi = NP.log10(vyhi)
    vxlo = NP.log10(vxlo);    vylo = NP.log10(vylo)

import math
vx = [ x for x in vx if (not math.isnan(x)) ]
vy = [ y for y in vy if (not math.isnan(y)) ]

#PLT.scatter(vxlo,vylo,color='red',alpha=0.5)
PLT.scatter(vxhi,vyhi,alpha=0.5)

PLT.xlabel(r'${\rm log}_{10}r_{\rm vir}/kpc$')
PLT.ylabel(r'${\rm log}_{10}M_{\rm vir}/M_\odot$')


#PLT.xlim(3.5,8.5)
#PLT.ylim(0.5,10000)

PLT.savefig(outfile)
