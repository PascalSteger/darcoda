#!/usr/bin/env ipython3
from __future__ import absolute_import, unicode_literals, print_function

import pymultinest
import math
import os
import subprocess

if not os.path.exists("chains"): os.mkdir("chains")
def show(filepath):
    """ open the output (pdf) file for the user """
    if os.name == 'mac': subprocess.call(('open', filepath))
    elif os.name == 'nt': os.startfile(filepath)
    elif os.name == 'posix': subprocess.call(('xdg-open', filepath))
## \fn show(filepath)
# call xdg-open
# @param filepath string



def myprior(cube, ndim, nparams):
    #print "cube before", [cube[i] for i in range(ndim)]
    for i in range(ndim):
        cube[i] = cube[i] * 10 * math.pi
    #print "python cube after", [cube[i] for i in range(ndim)]
## \fn myprior(cube, ndim, nparams)
# our probability functions
# Taken from the eggbox problem.
# @param cube [0,1]^ndim
# @param ndim int
# @param nparams int

def myloglike(cube, ndim, nparams):
    chi = 1.
    #print "cube", [cube[i] for i in range(ndim)], cube
    for i in range(ndim):
        chi *= math.cos(cube[i] / 2.)
    #print "returning", math.pow(2. + chi, 5)
    return math.pow(2. + chi, 5)
## \fn myloglike(cube, ndim, nparams)
# callback function for MultiNest
# @param cube physical parameter space
# @param ndim int
# @param nparams

if __name__=="__main__":
    # number of dimensions our problem has
    parameters = ["x", "y"]
    n_params = len(parameters)

    # run MultiNest
    pymultinest.run(myloglike,\
                    myprior,\
                    n_params,\
                    importance_nested_sampling = False,\
                    resume = True,\
                    verbose = True,\
                    sampling_efficiency = 'model',\
                    n_live_points = 1000,\
                    outputfiles_basename='chains/2-',\
                    init_MPI = False)
