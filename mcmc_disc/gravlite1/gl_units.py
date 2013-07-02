#!/usr/bin/python
# (c) 2013 Pascal S.P. Steger
'''work in pc, Msun, km/s'''

import gl_params as gp
import numpy as np

#    x, nu1, sig1 = units.get_physical(gp.xipol, gp.nu1_x, gp.sig1_x)
def get_physical(x, nu1, sig1, ntracer):
    return x * gp.rcore[ntracer], nu1*gp.dens0pc_2D[ntracer], sig1*gp.maxvlos[ntracer]
