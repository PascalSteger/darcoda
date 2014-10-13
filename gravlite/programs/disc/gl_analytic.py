#!/usr/bin/env ipython3

##
# @file
# @ingroup gravlite
# all analytic profiles from gl_physics

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import pdb
import gl_project as glp

asech = lambda x: np.arccosh(1./x)
asec  = lambda x: np.arccos(1./x)
