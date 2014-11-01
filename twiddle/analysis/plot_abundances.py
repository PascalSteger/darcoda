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

PLT.plot(z,v[1],c='black',label='e')
PLT.plot(z,v[2],c='red',label='HI')
PLT.plot(z,v[3],c='orange',label='HII')
PLT.plot(z,v[4],c='green',label='HeI')
PLT.plot(z,v[5],c='cyan',label='HeII')
PLT.plot(z,v[6],c='violet',label='HeIII')
PLT.plot(z,v[7],c='black',label='H-')
PLT.plot(z,v[8],c='blue',label='H2')
PLT.plot(z,v[9],c='brown',label='H2*')

PLT.xscale('log'); PLT.yscale('log')
PLT.xlim(8e1,10**4)
PLT.ylim(10**-15,10**1)
PLT.xlabel(r'z')
#PLT.xticks(NP.logspace(0,2,3),['1','10','100'])
PLT.ylabel(r'$n_i/n_{(HI+HII)}$')
#PLT.yticks(NP.logspace(-7,3,6))
PLT.legend(loc=2)

PLT.savefig(outfile)
