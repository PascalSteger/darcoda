#!/usr/bin/python
import numpy as np
import sys
import pdb
from pymc import *

from pylab import *
import gl_params as gp
import gr_params as gpr
from gl_helper import expDtofloat
from gl_class_files import *


size = 100
p = Uniform( "p", 0 , 1) #this is the fraction that come from mean1 vs mean2
ber = Bernoulli( "ber", p = p, size = size) # produces 1 with proportion p.

precision = Gamma('precision', alpha=0.1, beta=0.1)

mean1 = Normal( "mean1", 0, 0.001 ) #better to use normals versus Uniforms
                                        # (unless you are certain the value is truncated
mean2 = Normal( "mean2", 0, 0.001 )

@deterministic
def mean( ber = ber, mean1 = mean1, mean2 = mean2):
    return ber*mean1 + (1-ber)*mean2


# generate some artificial data (addition of two Gaussians)
# v = np.random.randint( 0, 2, size)
# data = v*(0+ np.random.randn(size) ) + (1-v)*(3 + np.random.randn(size ) )



print 'input:'
print gpr.fil
x0,y0,z0,vz0,vb0,Mg0,PM0,comp0=np.genfromtxt(gpr.fil,skiprows=0,unpack=True,\
                                             usecols=(0,1,2,11,12,13,19,20),\
                                             dtype="d17",\
                                             converters={0:expDtofloat,  # x0  in pc \
                                                         1:expDtofloat,  # y0  in pc \
                                                         2:expDtofloat,  # z0  in pc \
                                                         11:expDtofloat, # vz0 in km/s\
                                                         12:expDtofloat, # vb0(LOS due binary), km/s\
                                                         13:expDtofloat, # Mg0 in Angstrom\
                                                         19:expDtofloat, # PM0 [1]\
                                                         20:expDtofloat}) # comp0 1,2,3(background)
data = Mg0

obs = Normal( "obs", mean, precision, value = data, observed = True)
model = Model( {"p":p, "precision": precision, "mean1": mean1, "mean2":mean2, "obs":obs} )

if __name__=="__main__":
    from pymc import MCMC, Matplot

    S = MCMC(locals(), db='pickle')
    S.sample(iter=10000, burn=8000, thin=2)
    import pylab
    Matplot.plot(S)
    pylab.ion()
    pylab.hist(data,100)
    pylab.savefig('data.png')
    pylab.ioff()
    pylab.show()
