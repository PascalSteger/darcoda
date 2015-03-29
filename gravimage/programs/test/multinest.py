#!/usr/bin/env ipython3
from __future__ import absolute_import, unicode_literals, print_function
from ctypes import *

try:
    lib = CDLL('libnest3.so')
except OSError as e:
    print('error loading library')
    print('problem:', e)
    import sys
    sys.exit(1)

from numpy.ctypeslib import as_array
import signal, sys
import inspect

def interrupt_handler(signal, frame):
    sys.stderr.write('ERROR: Interrupt received: Terminating\n')
    sys.exit(1)

def run(LogLikelihood,
        Prior,
        nest_ndims,
        nest_totPar = None,
        nest_nCdims = None,
        wrapped_params = None,
        nest_IS = True,
        nest_mmodal = True,
        nest_ceff = False,
        nest_nlive = 400,
        nest_tol = 0.5,
        nest_ef = 0.8,
        nest_updInt = 100,
        null_log_evidence = -1e90,
        maxClst = 100,
        nest_Ztol = -1e30,
        outputfiles_basename = "chains/1-",
        seed = -1,
        nest_fb = False,
        nest_resume = True,
        context = 0,
        nest_outfile = True,
        nest_logZero = -1e100,
        nest_maxIter = 0,
        initMPI = True,
        dump_callback = None):

    print('got into run() function')
    if nest_totPar == None:
        nest_totPar = nest_ndims
    if nest_nCdims == None:
        nest_nCdims = nest_ndims
    if wrapped_params == None:
        wrapped_params = [0] * nest_ndims

    WrappedType = c_int * len(wrapped_params)
    nest_pWrap = WrappedType(*wrapped_params)

    if nest_ef == 'parameter':
        nest_ef = 0.8
    if nest_ef == 'model':
        nest_ef = 0.3

    loglike_type = CFUNCTYPE(c_double, POINTER(c_double),
        c_int, c_int, c_double, c_void_p)

    dumper_type  = CFUNCTYPE(c_void_p,         # nSamples
                             c_int,            # nlive
                             c_int,            # nPar
                             c_int,            # physLive
                             POINTER(c_double),# posterior
                             POINTER(c_double),# paramConstr
                             POINTER(c_double),# maxLogLike
                             c_double,         # logZ
                             c_double,         # INSlogZ
                             c_double,         # logZerr
                             c_void_p)         # context_pass

    # check if lnew is supported by user function
    nargs = 3
    try:
        nargs = len(inspect.getargspec(LogLikelihood).args)
    except:
        print('lnew not supported by user function')
        pass

    if nargs == 4:
        def loglike(cube, ndim, nparams, lnew, nullcontext):
            if Prior:
                Prior(cube, ndim, nparams)
            return LogLikelihood(cube, ndim, nparams, lnew)
    else:
        def loglike(cube, ndim, nparams, lnew, nullcontext):
            if Prior:
                Prior(cube, ndim, nparams)
            return LogLikelihood(cube, ndim, nparams)

    print('before dumper')
    def dumper(nSamples,nlive,nPar,
               physLive,posterior,paramConstr,
               maxLogLike,logZ,logZerr,nullcontext):
        if dump_callback:
            # It's not clear what the desired MultiNest dumper callback
            # syntax is... but this should pass back the right numpy arrays,
            # without copies. Untested!
            pc =  as_array(paramConstr,shape=(nPar,4))

            dump_callback(nSamples,nlive,nPar,
                as_array(physLive,shape=(nPar+1,nlive)).T,
                as_array(posterior,shape=(nPar+2,nSamples)).T,
                (pc[0,:],pc[1,:],pc[2,:],pc[3,:]), # (mean,std,bestfit,map)
                maxLogLike,logZ,logZerr)
    prev_handler = signal.signal(signal.SIGINT, interrupt_handler)

    print('nest_nlive = '+str(nest_nlive))

    print('before __nested_MOD_nestrun, here')

    lib.__nested_MOD_nestrun(c_bool(nest_IS),\
                             c_bool(nest_mmodal), \
                             c_bool(nest_ceff),\
                             c_int(nest_nlive), \
                             c_double(nest_tol),\
                             c_double(nest_ef),\
                             c_int(nest_ndims),\
                             c_int(nest_totPar),\
                             c_int(nest_nCdims),\
                             c_int(maxClst),\
                             c_int(nest_updInt),\
                             c_double(nest_Ztol),\
                             create_string_buffer(outputfiles_basename.encode(),100),\
                             c_int(seed), \
                             nest_pWrap,\
                             c_bool(nest_fb),\
                             c_bool(nest_resume),\
                             c_bool(nest_outfile),\
                             c_bool(initMPI),\
                             c_double(nest_logZero),\
                             c_int(nest_maxIter),\
                             loglike_type(loglike),\
                             dumper_type(dumper),\
                             c_int(context))
    print('started __nested_MOD_nest')
    signal.signal(signal.SIGINT, prev_handler)

## \fn run(LogLikelihood, Prior, nest_ndims, nest_totPar = None, nest_nCdims = None, wrapped_params = None, nest_IS = True, nest_mmodal = True, nest_ceff = False, nest_nlive = 400, nest_tol = 0.5, nest_ef = 0.8, nest_updInt = 100, null_log_evidence = -1e90, maxClst = 100, nest_Ztol = -1e30, outputfiles_basename = "chains/1-", seed = -1, nest_fb = False, nest_resume = True, context = 0, nest_outfile = True, nest_logZero = -1e100, nest_maxIter = 0, initMPI = True, dump_callback = None)
# Runs MultiNest
# The most important parameters are the two log-probability functions Prior
# and LogLikelihood. They are called by MultiNest.
# Prior should transform the unit cube into the parameter cube. Here
# is an example for a uniform prior::
#    def Prior(cube, ndim, nparams):
#        for i in range(ndim):
#            cube[i] = cube[i] * 10 * math.pi
# The LogLikelihood function gets this parameter cube and should
# return the logarithm of the likelihood.
# Here is the example for the eggbox problem::
#    def Loglike(cube, ndim, nparams, lnew):
#        chi = 1.
#        for i in range(ndim):
#            chi *= math.cos(cube[i] / 2.)
#        return math.pow(2. + chi, 5)
# Some of the parameters are explained below. Otherwise consult the
# MultiNest documentation.
# @param nest_IS:
# If True, Multinest will use Importance Nested Sampling (INS). Read http://arxiv.org/abs/1306.2144
# for more details on INS. Please read the MultiNest README file before using the INS in MultiNest v3.0.
# @param nest_totPar:
#    Total no. of parameters, should be equal to ndims in most cases
#    but if you need to store some additional
#    parameters with the actual parameters then you need to pass
#    them through the likelihood routine.
# @param maxClst ?
# @param nest_ef:
#    defines the sampling efficiency. 0.8 and 0.3 are recommended
#    for parameter estimation & evidence evalutation
#    respectively.
#    use 'parameter' or 'model' to select the respective default
#    values
# @param nest_Ztol:
#    MultiNest can find multiple modes & also specify which samples belong to which mode. It might be
#    desirable to have separate samples & mode statistics for modes with local log-evidence value greater than a
#    particular value in which case Ztol should be set to that value. If there isn't any particularly interesting
#    Ztol value, then Ztol should be set to a very large negative number (e.g. -1e90).
# @param nest_tol:
#    A value of 0.5 should give good enough accuracy.
# @param nest_nCdims:
#    If mmodal is T, MultiNest will attempt to separate out the
#    modes. Mode separation is done through a clustering
#    algorithm. Mode separation can be done on all the parameters
#    (in which case nCdims should be set to ndims) & it
#    can also be done on a subset of parameters (in which case
#    nCdims < ndims) which might be advantageous as
#    clustering is less accurate as the dimensionality increases.
#    If nCdims < ndims then mode separation is done on
#    the first nCdims parameters.
# @param null_log_evidence:
#    If mmodal is T, MultiNest can find multiple modes & also specify
#    which samples belong to which mode. It might be
# desirable to have separate samples & mode statistics for modes
# with local log-evidence value greater than a
# particular value in which case nullZ should be set to that
# value. If there isn't any particulrly interesting
# nullZ value, then nullZ should be set to a very large negative
# number (e.g. -1.d90).
# @param initMPI:
#     initialize MPI routines?, relevant only if compiling with MPI
# @param nest_logZero:
#     points with loglike < logZero will be ignored by MultiNest
# @param nest_maxIter:
#     maximum number of iterations. 0 is unlimited.
# @param nest_outfile:
#     write output files? This is required for analysis.
# @param dump_callback:
#     a callback function for dumping the current status
