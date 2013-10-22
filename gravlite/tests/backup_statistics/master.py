#!/usr/bin/python

import pymc
import pymcgau
import pylab

S = pymc.MCMC(pymcgau, db='pickle')
S.sample(iter=10000, burn=5000, thin=2)
pylab.ion()
pymc.Matplot.plot(S)
pylab.show()
