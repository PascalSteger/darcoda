#!/usr/bin/env ipython3

##
# @file
# all integrals from gl_physics

# (c) 2013 Pascal S.P. Steger

import numpy as np
import pdb, scipy
from scipy.integrate import simps,trapz,quad
from scipy.interpolate import splrep, splev, splint
import gl_helper as gh
import gl_physics as phys
from pylab import *
