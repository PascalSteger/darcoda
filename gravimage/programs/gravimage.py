#!/usr/bin/env python3

##
# @file
# pymultinest run of gravimage integrals needs pymultinest from
# http://johannesbuchner.github.io/PyMultiNest/
# http://johannesbuchner.github.io/PyMultiNest/install.html#install-on-linux
# needs Multinest from https://github.com/JohannesBuchner/MultiNest

# TODO: run with mpirun -np <N> gravimage.py
# where <N> is an integer <= number of processors
# (c) 2013 ETHZ Pascal Steger, psteger@phys.ethz.ch

### imports
# from __future__ import absolute_import, unicode_literals, print_function
from mpi4py import MPI
import subprocess
import numpy as np
import pymultinest
import pickle
import gl_helper as gh
import pdb
import numpy.random as npr
import time
import gl_multinest_helper as glmh
# increment NICEness of process by 1, CPU usage shall not block others
# import os
# os.nice(1)

# optionally start with -i and -c switches, to batch start gaia and walk runs
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--investigation", dest="investigation",
                      default="", help="investigation to run: simplenu, disc_nbody")
parser.add_option("-c", "--case", dest="case",
                      default=-1, help="case: 1 (Hunt & Kawata GCD), 2 (Garbari N-body), ..")
(options, args) = parser.parse_args()
print('gravimage.py '+str(options.investigation)+' '+str(options.case))
import gl_params
import warnings
warnings.simplefilter('ignore') # set to 'error' when debugging
ts = '' # empty timestamp means: create new timestamp with folder
#gp = gl_params.Params(ts, options.investigation, int(options.case))
#import gl_file as gf

def show(filepath):
    subprocess.call(('xdg-open', filepath))
    return
## \fn show(filepath) open the output (pdf) file for the user @param
# filepath filename with full pathThen

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
    myrank = MPI.COMM_WORLD.Get_rank()

    if min(cube[0:ndim]) == -9999:  # parameters not fulfilling prior requirements,
        return -1e300       #      return very large chi2
    try:
        tmp_profs = geom_loglike(cube, ndim, nparams, gp)
    except ValueError:
        return -1e100  # SS: think about sign...
    # store tmp_prof by appending it to pc2.save
    # TODO: with parallel version, need to append to CPU-based output name
    # we only store models after the initial Sigma burn-in
    if npr.random() < gp.save_fraction: #save only a fraction of models
        tmp_profs.x0 = gp.z_bincenters
        tmp_profs.xbins = np.hstack([gp.z_binmins, gp.z_binmaxs[-1]])
        with open(gp.files.outdir+'pc2.save', 'ab') as fi:
            pickle.dump(tmp_profs, fi)
            # convention: use chi^2 directly, not log likelihood
    # for output:
    # from   likelihood L = exp(-\chi^2/2), want log of that
    #print('P', myrank, ': tmp_profs.chi2 = ', tmp_profs.chi2)
    return -tmp_profs.chi2/2.
## \fn myloglike(cube, ndim, nparams) calculate probability function
# @param cube ndim cube of physical parameter space (nr)
# @param ndim number of dimensions, 2*npop*nipol + nipol
# @param nparams = ndim + additional parameters
# stored with actual parameters

def prepare_data(gp):
    hwmess = "prepare_data, process %d of %d on %s.\n"
    myrank = MPI.COMM_WORLD.Get_rank()
    nprocs = MPI.COMM_WORLD.Get_size()
    procnm = MPI.Get_processor_name()
    import sys
    sys.stdout.write(hwmess % (myrank, nprocs, procnm))

    if gp.getnewdata:
        if gp.getnewpos:
            gf.read_data(gp)
        gf.bin_data(gp)
    gf.get_binned_data_noscale(gp)    #H Silverwood 20/11/14
    gp.files.populate_output_dir(gp)
    gf.get_rhohalfs(gp)
    print('P', myrank, ': prepare_data complete')




if __name__=="__main__":
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    nprocs = comm.Get_size()
    procnm = MPI.Get_processor_name()

    #Instantiate the gp class on P0, hold other processes until completed
    if myrank == 0:
        print('P', myrank, ':, creating gp.')
        gp = gl_params.Params(ts, options.investigation, int(options.case))
        print('P', myrank, ':, created gp.')
        glmh.write_mn_info(gp)
        comm.Barrier()

    elif myrank != 0:
        print('P', myrank, ':, On hold while gp class is instantiated')
        gp = None
        comm.Barrier()

    #Broadcase gp class to other processors
    gp = comm.bcast(gp, root=0)

    #Set import paths for P1 thru PN (this is done for P0 when instantiating gp,
    # but does not propagate to other processes). TO DO fix this.
    import import_path as ip
    ip.set_geometry(gp.geom)

    import gl_file as gf
    global Cube, geom_loglike
    from gl_class_cube import Cube
    from gl_loglike import geom_loglike

    #Prepare data on P0, hold other processes until completed
    if myrank == 0:
        print('P', myrank, ': Preparing Data')
        prepare_data(gp)
        print('P', myrank, ': Data prepared\n')
        time.sleep(5)
        comm.Barrier()

    elif myrank != 0:
        print('P', myrank, ': On hold while data prepared\n')
        comm.Barrier()

    #Broadcase gp again
    gp = comm.bcast(gp, root=0)

    print('\n P', myrank, ': Running pymultinest')

    start_time = time.time()
    pymultinest.run(myloglike,   myprior,
                    gp.ndim, n_params = gp.ndim+1, # None beforehands
                    n_clustering_params = gp.ndim,# separate modes on
                                                  # the rho parameters
                                                  # only: gp.nrho
                    wrapped_params = [gp.ntracer_pops, gp.nbins[0], gp.nrhonu], # do
                                                                     #not
                                                                     #wrap-around
                                                                     #parameters
                    importance_nested_sampling = False, # INS enabled
                    multimodal = True,           # separate modes
                    const_efficiency_mode = False, # use const sampling efficiency
                    n_live_points = gp.nlive,
                    evidence_tolerance = 0.5, # set to 0 to keep
                                              # algorithm working
                                              # indefinitely
                    sampling_efficiency = 0.8,#0.8 or 'parameter' for parameter estimation,
                                                #0.3 or 'evidence' for evidence calculation
                    n_iter_before_update = 100, # output after this many iterations
                    null_log_evidence = -1e100,
                    max_modes = gp.nlive,   # preallocation of modes:
                                            #max. = number of live
                                            #points
                    mode_tolerance = -1.e100,   # mode tolerance in the
                                               #case where no special
                                               #value exists: highly
                                               #negative
                    outputfiles_basename = gp.files.outdir + '/output',
                    seed = -1, # -1 for seed from clock
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
                    init_MPI = False,     # use MPI
                    dump_callback = None)


    end_time = time.time()
    run_time = end_time - start_time
    print('\n P', myrank, ': On ', nprocs, ' cores, run time is ',  run_time, ' sec')
