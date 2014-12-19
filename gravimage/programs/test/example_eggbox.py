#!/usr/bin/env python3
import math
import sys

try:
    from mpi4py import MPI

    myrank = MPI.COMM_WORLD.Get_rank()
    nprocs = MPI.COMM_WORLD.Get_size()
    procnm = MPI.Get_processor_name()
except:
    myrank = 0
    nprocs = 1
    procnm = 'localhost'

sys.stdout.write("Hello, World!! I am process %d of %d on %s.\n" % (myrank, nprocs, procnm))
sys.stdout.flush()

# print os.environ['SHELL']
# print os.environ['LD_LIBRARY_PATH']

import pymultinest

# this function needs to do two things:
# - get the [0,1] parameters in cube and scale them in place to their physical
#   values (implementing priors if needed)
# - return the likelihood
#
# note that len(cube) may be larger than ndim; this allows returning
# "nonphysical" parameters (not used in the search, but reported in output)

evals = 0

def LogLike(cube, ndim):
    global evals
    evals = evals + 1

    chi = 1.0

    for i in range(ndim):
        cube[i] = 10.0 * math.pi * cube[i]
        chi *= math.cos(0.5 * cube[i])

    return (chi + 2.0)**5


mmodal = True       # search for multiple modes
ceff = 0            # SMBH tuning
nlive = 1000        # number of "live" points
tol = 0.5           # final tolerance in computing evidence
efr = 0.9           # sampling efficiency (0.8 and 0.3 are recommended for parameter estimation & evidence evaluation)
ndims = 2           # number of search parameters
nPar = 2            # number of reported parameters (see above)
nClsPar = 2         # number of parameters on which to attempt mode separation (first nClsPar parameters used)
maxModes = 50       # maximum number of modes
updInt = 100        # interval of updates
Ztol = -1.e90       # threshold for reporting evidence about modes
root = "chains/1-"  # prefix for output files
seed = 13           # seed for random numbers
periodic = [1,1]    # period conditions
fb = True           # output status updates?
resume = False      # resume from previous run?

pymultinest.nested.nestRun(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, periodic, fb, resume, LogLike, 0)

sys.stdout.write("Hello, World!! I am process %d of %d on %s, and I ran %d evaluations.\n" % (myrank, nprocs, procnm, evals))
sys.stdout.flush()
