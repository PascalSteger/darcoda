#!/usr/bin/env python3

##
# @file
# all integrals from gl_physics

# (c) GPL v3 2014 Pascal S.P. Steger

import numpy as np
import pdb, scipy
from scipy.integrate import simps,trapz,quad
from scipy.interpolate import splrep, splev, splint
import gl_helper as gh
import gl_physics as phys
from pylab import *
ion()
