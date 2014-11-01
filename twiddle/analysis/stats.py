#!/usr/bin/env python3

## \file
# plot histograms from halos file

import sys
i=len(sys.argv)
if i<4:
    print("Usage: stats.py infile outfile column(start from 1) [log]")
    print("Example: stats.py halos stat.png 4 1")
    exit(1)
    
infile  = sys.argv[1]; outfile = sys.argv[2]
col     = [int(sys.argv[3])-1]
boolog  = False
if(i>4):
    boolog  = int(sys.argv[4])>0

from matplotlib import pyplot as PLT
fig = PLT.figure()
ax1 = fig.add_subplot(111)

import numpy as NP
with open(infile) as f:
    v = NP.loadtxt(f, dtype='float', comments="#", skiprows=0, usecols=col)#delimiter=","
    v = NP.ravel(v)   # 'flatten' v, cumulative sum

if(boolog):
    v = NP.log10(v)
import math
v = [ x for x in v if not math.isnan(x) ]
#print('v = ',v)

import scipy
from scipy.stats import scoreatpercentile
q1    = scoreatpercentile(v,25)
q3    = scoreatpercentile(v,75)
iqd   = q3 - q1
dv    = 2*iqd*len(v)**(-1/3.)
nbins = (max(v)-min(v))/dv

nbins = 21
counts,bins,rects = PLT.hist(v,normed=False,log=boolog,histtype='step',bins=nbins)
#PLT.xlim(3.5,8.5)
#PLT.ylim(0.5,10000)
PLT.xlabel(r'${\rm log}_{10}M/M_\odot$')
PLT.ylabel(r'count')

print(counts)
print("sum(counts) = ",sum(counts))
width=0.7*(bins[1]-bins[0])
center=(bins[:-1]+bins[1:])/2
print("bin centers: ",center)

hist=1.0*counts/sum(counts)
print("normalized counts: ",hist)

#PLT.clf()
#PLT.step(center,NP.log10(hist))
#print("area: ", sum(NP.diff(r[1])*r[0]))


#counts,bins,patches=PLT.hist(v, density=True, bins=10, log=True)
# hist         = counts/sum(counts * NP.diff(bins))
# print(sum(counts/sum(counts * NP.diff(bins))))

# PLT.clf()
# PLT.bar(center,counts/sum(counts),align='center')

PLT.savefig(outfile)
