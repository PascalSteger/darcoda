#!/usr/bin/env python3

##
# @file
# pymultinest run of gravlite integrals needs pymultinest from
# http://johannesbuchner.github.io/PyMultiNest/
# http://johannesbuchner.github.io/PyMultiNest/install.html#install-on-linux
# needs Multinest from https://github.com/JohannesBuchner/MultiNest

# TODO: run with mpirun -np <N> gravlite.py
# where <N> is an integer <= number of processors
# (c) 2013 ETHZ Pascal Steger, psteger@phys.ethz.ch

### imports
# from __future__ import absolute_import, unicode_literals, print_function
#from mpi4py import MPI
import subprocess
import numpy as np
import pymultinest
import pickle
import pdb
# increment NICEness of process by 1, CPU usage shall not block others
# import os
# os.nice(1)

# optionally start with -i and -c switches, to batch start gaia and walk runs
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--investigation", dest="investigation",
                      default="", help="investigation to run: gaia, walk, hern, triax, discmock, simplenu")
parser.add_option("-c", "--case", dest="case",
                      default=-1, help="case: 1, 2, ..")
(options, args) = parser.parse_args()
print('gravlite.py '+str(options.investigation)+' '+str(options.case))
import gl_params
import warnings
warnings.simplefilter('ignore') # set to 'error' when debugging
ts = '' # empty timestamp means: create new timestamp with folder
gp = gl_params.Params(ts, options.investigation, int(options.case))
import gl_file as gf

def show(filepath):
    subprocess.call(('xdg-open', filepath))
    return
## \fn show(filepath) open the output (pdf) file for the user @param
# filepath filename with full path

def myprior(cube, ndim, nparams):
    mycube = Cube(gp)
    mycube.copy(cube)
    try:
        cube = mycube.convert_to_parameter_space(gp)
    except Exception:
        gh.LOG(1, 'parameters not fulfilling prior requirements')

    return
## \fn myprior(cube, ndim, nparams) priors
# @param cube [0,1]^ndim cube, array of dimension ndim
# @param ndim number of dimensions, 2*npop*nipol + nipol
# @param nparams = ndim + additional parameters
# stored with actual parameters

def myloglike(cube, ndim, nparams):
    if min(cube) == -9999:  # parameters not fulfilling prior requirements,
        return -1e300       #      return very large chi2
    tmp_profs = geom_loglike(cube, ndim, nparams, gp)
    # store tmp_prof by appending it to pc2.save
    # TODO: with parallel version, need to append to CPU-based output name
    # we only store models after the initial Sigma burn-in
    if gp.chi2_nu_converged:
        tmp_profs.x0 = gp.xipol
        tmp_profs.xbins = np.hstack([gp.binmin, gp.binmax[-1]])
        with open(gp.files.outdir+'pc2.save', 'ab') as fi:
            pickle.dump(tmp_profs, fi)
            # convention: use chi^2 directly, not log likelihood
    # for output:
    # from   likelihood L = exp(-\chi^2/2), want log of that
    return -tmp_profs.chi2/2.
## \fn myloglike(cube, ndim, nparams) calculate probability function
# @param cube ndim cube of physical parameter space (nr)
# @param ndim number of dimensions, 2*npop*nipol + nipol
# @param nparams = ndim + additional parameters
# stored with actual parameters

def prepare_data(gp):
    if gp.getnewdata:
        if gp.getnewpos:
            gf.read_data(gp)
        gf.bin_data(gp)
    gf.get_binned_data_noscale(gp)    #H Silverwood 20/11/14
    gp.files.populate_output_dir(gp)
    gf.get_rhohalfs(gp)
## \fn prepare_data(gp)
# prepare everything for multinest(.MPI) run
# @param gp global parameters
def run(gp):
    pymultinest.run(myloglike,   myprior,
                    gp.ndim, n_params = gp.ndim+1, # None beforehands
                    n_clustering_params = gp.ndim,# separate modes on
                                                  # the rho parameters
                                                  # only: gp.nrho
                    wrapped_params = [ gp.pops, gp.nipol, gp.nrho], # do
                                                                     #not
                                                                     #wrap-around
                                                                     #parameters
                    importance_nested_sampling = True, # INS enabled
                    multimodal = True,           # separate modes
                    const_efficiency_mode = True, # use const sampling efficiency
                    n_live_points = gp.nlive,
                    evidence_tolerance = 0.0, # set to 0 to keep
                                              # algorithm working
                                              # indefinitely
                    sampling_efficiency = 0.95,
                    n_iter_before_update = 100, # output after this many iterations
                    null_log_evidence = -1e100,
                    max_modes = gp.nlive,   # preallocation of modes:
                                            #max. = number of live
                                            #points
                    mode_tolerance = -1.e100,   # mode tolerance in the
                                               #case where no special
                                               #value exists: highly
                                               #negative
                    outputfiles_basename = gp.files.outdir,
                    seed = -1,
                    verbose = True,
                    resume = False,
                    context = 0,
                    write_output = True,
                    log_zero = -1e500,    # points with log likelihood
                                          #< log_zero will be
                                          #neglected
                    max_iter = 0,         # set to 0 for never
                                          #reaching max_iter (no
                                          #stopping criterium based on
                                          #number of iterations)
                    init_MPI = True,     # use MPI
                    dump_callback = None)

if __name__=="__main__":
    global Cube, geom_loglike
    from gl_class_cube import Cube
    from gl_loglike import geom_loglike

    # hwmess = "Hello, World!! I am process %d of %d on %s.\n"
    # myrank = MPI.COMM_WORLD.Get_rank()
    # nprocs = MPI.COMM_WORLD.Get_size()
    # procnm = MPI.Get_processor_name()
    # import sys
    # sys.stdout.write(hwmess % (myrank, nprocs, procnm))

    # TODO: wait for prepare_data to finish
    # if MPI.COMM_WORLD.Get_rank() == 0:
    #     # TODO: wrong: have 0 take part in sampling as well
    prepare_data(gp) # run once
    # else:
    #     # run with full MPI
    run(gp)
