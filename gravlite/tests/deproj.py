#!/usr/bin/python3.2
# (c) 2013 Pascal Steger, ETH Zurich, psteger@phys.ethz.ch

import numpy as np
# import gl_int as gi
# import physics_sphere as phys


r0 = np.arange(0.,3.,0.01)
rho = np.exp(-r0)

def rho(r):
    return np.exp(-1.*r)

from scipy.special import kv
def th_Sigma(R):
    return 2.*R*kv(1,R)

def igra(r,R):
    return 2.*r*rho(r)/np.sqrt(r**2-R**2)

from scipy.integrate import quad

def intSigma(R):
    return quad(igra,R,np.inf,args=(R),limit=100,full_output=0)

from pylab import *
ion()
figure()
plot(r0, th_Sigma(r0))
draw()

import pdb
pdb.set_trace()

ioff()
show()

