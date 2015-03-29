#!/usr/bin/env python3

##
# @file
# needs Multinest from https://github.com/JohannesBuchner/MultiNest

# Copyright 2015 GPLv3 ETHZ, Pascal Steger, pascal@steger.aero

# get latest version of PyMultiNest from
# svn checkout http://lisasolve.googlecode.com/svn/trunk/ lisasolve-read-only

import subprocess
import pymultinest
import pickle
import warnings
import numpy as np
import pdb
#from multiprocessing import Pool

# increment NICEness of process by 1, if CPU usage shall not block others
# import os
# os.nice(1)

# optionally start with -i and -c switches, to batch start gaia and walk runs
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--investigation", dest="investigation", default="", help="investigation to run: gaia, walk, hern, triax, discmock")
parser.add_option("-c", "--case", dest="case", default='-1', help="case: 1, 2, ..")
parser.add_option("-t", "--timestamp", dest="timestamp", default='-1', help="timestamp: 201501221224")
(options, args) = parser.parse_args()
print('gravimage.py '+options.investigation+' '+str(options.case)+' '+str(options.timestamp))
if options.timestamp != '-1':
    import gi_base as gb
    basepath = gb.get_basepath()
    import import_path as ip
    ip.insert_sys_path(basepath+"DT"+options.investigation+"/"+options.case+"/"+options.timestamp+"/programs/")
import gi_params
warnings.simplefilter('ignore') # set to 'error' when debugging
gp = gi_params.Params(options.timestamp, options.investigation, int(options.case))
if options.timestamp != '-1':
    ip.remove_third()
    #import os
    #os.system('cd '+basepath+'DT'+options.investigation+'/'+options.case+'/'+options.timestamp)
    gp.restart = True
    gp.chi2_Sig_converged = 0
import gi_file as gf

def show(filepath):
    subprocess.call(('xdg-open', filepath))
    return
## \fn show(filepath) open the output (pdf) file for the user @param
# filepath filename with full path

def myprior(cube, ndim, nparams):
    mycube = Cube(gp)
    mycube.copy(cube)
    cube = mycube.convert_to_parameter_space(gp)
    return
## \fn myprior(cube, ndim, nparams)
# priors
# @param cube [0,1]^ndim cube, array of dimension ndim
# @param ndim number of dimensions, 2*npop*nipol + nipol
# @param nparams = ndim + additional parameters
# stored with actual parameters

def myloglike(cube, ndim, nparams):
    tmp_profs = geom_loglike(cube, ndim, nparams, gp)
    # store tmp_prof by appending it to pc2.save
    # we only store models after the initial Sigma burn-in
    if gp.chi2_Sig_converged <= 0:
        tmp_profs.x0 = gp.xepol
        tmp_profs.xbins = np.hstack([gp.dat.binmin, gp.dat.binmax[-1]])
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
        gf.get_pos_and_COM(gp)
        gf.bin_data(gp)
    if gp.getSigdata:
        # if Sig convergence finished already
        gf.read_Sigdata(gp)
    gf.get_binned_data(gp)
    if not gp.restart:
        gp.files.populate_output_dir(gp)
    gf.get_rhohalfs(gp)
## \fn prepare_data(gp)
# prepare everything for multinest(.MPI) run
# @param gp global parameters

def run(gp):
    pymultinest.run(myloglike, myprior, gp.ndim, n_params = gp.ndim+1,
                    n_clustering_params = gp.nrho, # gp.ndim, or separate modes on the rho parameters only: gp.nrho
                    wrapped_params = [ gp.pops, gp.nipol, gp.nrho],
                    importance_nested_sampling = False, # INS enabled
                    multimodal = False,           # separate modes
                    const_efficiency_mode = True, # use const sampling efficiency
                    n_live_points = gp.nlive,
                    evidence_tolerance = 0.0,   # 0 to keep algorithm working indefinitely
                    sampling_efficiency = 0.05, # 0.05, MultiNest README for >30 params
                    n_iter_before_update = 2,  # output after this many iterations
                    null_log_evidence = -1e100,
                    max_modes = gp.nlive, # preallocation of modes: max=number of live points
                    mode_tolerance = -1.e100,   # mode tolerance in the case where no special value exists: highly negative
                    outputfiles_basename = gp.files.outdir,
                    seed = -1, verbose = True,
                    resume = gp.restart,
                    context = 0, write_output = True,
                    log_zero = -1e500, # points with log likelihood<log_zero will be neglected
                    max_iter = 0, # set to 0 for never reaching max_iter (no stopping criterium based on number of iterations)
                    init_MPI = False, dump_callback = None)

if __name__=="__main__":
    global Cube, geom_loglike
    from gi_class_cube import Cube
    from gi_loglike import geom_loglike
    prepare_data(gp)
    run(gp)
