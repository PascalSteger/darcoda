#!/usr/bin/env python3

## \file
# plot histograms from halos file

import sys
i=len(sys.argv)
if i<4:
    print("Usage: stats_fmhi.py infile outfile column(start from 1) fmsplit")
    print("Example: stats_fmhi.py halos stat.png 4 0.8")
    exit(1)
    
infile  = sys.argv[1]; outfile = sys.argv[2]
col     = [int(sys.argv[3])-1]
boolog  = True
if(i>4):
    fmsplit=float(sys.argv[4])
else:
    fmsplit=0.0
print("fmsplit = ",fmsplit)


import numpy as NP
with open(infile) as f:
    v = NP.loadtxt(f, dtype='float', comments="#", skiprows=0, usecols=col)#delimiter=","
    v = NP.ravel(v)   # 'flatten' v, cumulative sum
with open(infile) as g:
    w = NP.loadtxt(g, dtype='float', comments="#", skiprows=0, usecols=[37])
    w = NP.ravel(w)

if(boolog):
    v = NP.log10(v)
print('len(v) = ',len(v))
vlo=v[w<=fmsplit]
vhi=v[w>fmsplit]

print("lengths:",len(vlo),len(vhi))

# import scipy
# from scipy.stats import scoreatpercentile
#q1    = scoreatpercentile(vlo,25); q3    = scoreatpercentile(vlo,75);iqd   = q3 - q1
#dv    = 2*iqd*len(vlo)**(-1/3.);   nbins = (max(vlo)-min(vlo))/dv

nbins=12

from pylab import *
ion(); subplot(111)

#
#hist(vlo,log=False,histtype='step',bins=nbins,color="red")#normed=normed
counts,bins,rects=hist(vhi,bins=nbins,color="blue")#normed=normed

xlabel(r'${\rm log}_{10}M/M_\odot$')
ylabel(r'count')

print(counts)
print(bins)
# width=0.7*(bins[1]-bins[0])
# center=(bins[:-1]+bins[1:])/2
# hist=1.0*counts/sum(counts)

#bar(center,counts/sum(counts),align='center')

savefig(outfile)
ioff();show()
