#!/usr/bin/env python3
from numpy cimport ndarray
from numpy import empty

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

cdef extern:
    dumper()


cdef extern:
    loglike()

cdef extern:
    void wrap_nestrun(bint* nest_IS, bint* nest_mmodal, bint* nest_ceff,
                      int* nest_nlive, double* nest_tol, double* nest_ef,
                      int* nest_ndims, int* nest_totPar, int* nest_nCdims,
                      int* maxClst, int* nest_updInt, double* nest_Ztol,
                      char* str, int* seed, int* nest_pWrap, bint* nest_fb,
                      bint* nest_resume, bint* nest_outfile, bint* initMPI,
                      double* nest_logZero, int* nest_maxIter, int* context)

def run(bint nest_IS, bint nest_mmodal, bint nest_ceff,
                      int nest_nlive, double nest_tol, double nest_ef,
                      int nest_ndims, int nest_totPar, int nest_nCdims,
                      int maxClst, int nest_updInt, double nest_Ztol,
                      char str, int seed, int nest_pWrap, bint nest_fb,
                      bint nest_resume, bint nest_outfile, bint initMPI,
                      double nest_logZero, int nest_maxIter, int context):
    wrap_nestrun(&nest_IS, &nest_mmodal, &nest_ceff,
                 &nest_nlive, &nest_tol, &nest_ef,
                 &nest_ndims, &nest_totPar, &nest_nCdims,
                 &maxClst, &nest_updInt, &nest_Ztol,
                 &str, &seed, &nest_pWrap, &nest_fb,
                 &nest_resume, &nest_outfile, &initMPI,
                 &nest_logZero, &nest_maxIter, &context   )

#def mesh_exp(double r_min, double r_max, double a, int N):
#    cdef ndarray[double, mode="c"] mesh = empty(N, dtype=double)
#    c_mesh_exp(&r_min, &r_max, &a, &N, &mesh[0])
#    return mesh


import math, os
import multinest


# number of dimensions our problem has
parameters = ["x", "y"]
n_dims = len(parameters)

print('reached state before run')
# run MultiNest
multinest.run(myloglike,\
              myprior,\
              n_dims,\
              nest_resume = True, \
              nest_fb = True,\
              initMPI = False)

print('run successfully finished')
