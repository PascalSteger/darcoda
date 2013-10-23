#!/usr/bin/env python3

##
# @file
# test rho_INT_Rho projection with mock data

# (c) 2013 psteger@phys.ethz.ch

import numpy as np
import pdb

import gl_params
import gl_plot as gpl
from gl_project import *
from gl_analytic import *

r0 = np.arange(0.1,3.0,0.1)
rho0 = rho_anf(r0,1.,1.e6)
Rho0 = Sigma_anf(r0,1.,1.e6)
gpl.ion()
gpl.start()
# gpl.yscale('linear')

gpl.plot(r0, Rho0, color='blue', lw=2)
# gpl.plot(r0, rho_INTDIRECT_Rho(r0, rho0), color='blue')
try:
    Rhonu = rho_INT_Rho(r0,rho0)
except Exception as detail:
    print('NaN in rho_INT_Rho')
gpl.plot(r0, Rhonu, '.', color='red')
# pdb.set_trace()
print('error fraction determined/true = ',Rhonu/Rho0)
print('gp.ascale = ',gp.ascale)
print('gp.Mscale = ',gp.Mscale)
gpl.ioff(); gpl.show()
