#!/usr/bin/env python3

import pymultinest
import math
import os
if not os.path.exists("chains"): os.mkdir("chains")

# our probability functions
# Taken from the eggbox problem.

def myprior(cube, ndim, nparams):
	for i in range(ndim):
		cube[i] = cube[i] * 10 * math.pi

def myloglike(cube, ndim, nparams):
	chi = 1.
	for i in range(ndim):
		chi *= math.cos(cube[i] / 2.)
	return math.pow(2. + chi, 5)

# number of dimensions our problem has
parameters = ["x", "y"]
n_params = len(parameters)

# run MultiNest
pymultinest.run(LogLikelihood=myloglike,
              Prior=myprior,
              nest_ndims = n_params,
              nest_fb = True,
              nest_resume = False)


# run
# $ multinest_marginals.py chains/1-
# which will produce pretty marginal pdf plots

# for code to analyse the results, and make plots see full demo
