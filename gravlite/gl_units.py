#!/usr/bin/env python3

##
# @file
# work in pc, Msun, km/s

# (c) 2013 Pascal S.P. Steger, psteger@phys.ethz.ch


import gl_params as gp
import numpy as np

## x, nu1, sig1 = units.get_physical(gp.xipol, gp.nu1_x, gp.sig1_x)
# @param x radius in [rcore]
# @param nu1 2D tracer density in [dens0pc_2D]
# @param sig1 velocity dispersion profile in [maxvlos]
# @param ntracer population
def get_physical(x, nu1, sig1, ntracer):
    return x * gp.rcore[ntracer], nu1*gp.dens0pc_2D[ntracer], sig1*gp.maxvlos[ntracer]
