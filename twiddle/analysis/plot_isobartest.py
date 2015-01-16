#!/usr/bin/env python3

## \file
# plot abundances from cat rectest.log|grep "Y:"|cut -d":" -f2|pr -s -t -l9|tee rectest.col

import sys
infile  = sys.argv[1]; outfile = sys.argv[2]

from matplotlib import pyplot as PLT
fig = PLT.figure()
ax1 = fig.add_subplot(111)

import numpy as NP
with open(infile) as f:
    v = NP.loadtxt(f, dtype='float', comments="#", skiprows=0, unpack=True)#delimiter=",", usecols=[col]
    print(v)

import math
import scipy

z = 1/v[0]-1

PLT.plot(z,v[4],c='black',label='e')
PLT.plot(z,v[5],c='red',label='HI')
PLT.plot(z,v[6],c='orange',label='HII')
PLT.plot(z,v[7],c='green',label='HeI')
PLT.plot(z,v[8],c='cyan',label='HeII')
PLT.plot(z,v[9],c='violet',label='HeIII')
PLT.plot(z,v[10],c='black',label='H-')
PLT.plot(z,v[11],c='blue',label='H2')
PLT.plot(z,v[12],c='brown',label='H2*')

PLT.xscale('log'); PLT.yscale('log')
#PLT.xlim(0.5,512)
#PLT.ylim(10**-24,10**0)
PLT.xlabel(r'z')
#PLT.xticks(NP.logspace(0,2,3),['1','10','100'])
PLT.ylabel(r'$Y_i$')
#PLT.yticks(NP.logspace(-7,3,6))
PLT.legend(loc=2)

PLT.savefig(outfile)
