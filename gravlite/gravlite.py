#!/usr/bin/env python3

##
# @file
# pymultinest run of gravlite integrals
# needs pymultinest from http://johannesbuchner.github.io/PyMultiNest/
# http://johannesbuchner.github.io/PyMultiNest/install.html#install-on-linux
# needs Multinest from https://github.com/JohannesBuchner/MultiNest

# (c) 2013 ETHZ Pascal Steger, psteger@phys.ethz.ch

from __future__ import absolute_import, unicode_literals, print_function
import numpy as np
import pdb
import threading, subprocess
import pymultinest


def show(filepath):
    subprocess.call(('xdg-open', filepath))
    return
## \fn show(filepath)
# open the output (pdf) file for the user
# @param filepath filename with full path


def myprior(cube, ndim, nparams):
    mycube = Cube(gp)
    mycube.copy(cube)
    cube = mycube.convert_to_parameter_space(gp)
    return
## \fn myprior(cube, ndim, nparams)
# priors
# @param cube [0,1]^ndim cube, array of dimension ndim
# @param ndim number of dimensions, 2*npop*nipol + nipol
# @param nparams = ndim + additional parameters stored with actual parameters


def myloglike(cube, ndim, nparams):
    return geom_loglike(cube, ndim, nparams, gp)
## \fn myloglike(cube, ndim, nparams)
# calculate probability function
# @param cube [0,1]^ndim cube
# @param ndim number of dimensions, 2*npop*nipol + nipol
# @param nparams = ndim + additional parameters stored with actual parameters


def run(gp):
    gp.files.create_output_dir(gp)
    
    import gl_file as gfile
    if gp.getnewdata:
        gfile.bin_data(gp)
    gfile.get_data(gp)

    gp.files.populate_output_dir(gp)
    
    ## number of dimensions
    n_dims = gp.nepol + gp.pops*gp.nepol + gp.pops*gp.nbeta #rho, (nu, beta)_i
    
    # show live progress
    # progress = pymultinest.ProgressPlotter(n_params = n_dims)
    # progress.start()
    # threading.Timer(2, show, [gp.files.outdir+'/phys_live.points.pdf']).start() 

    # print(str(len(gp.files.outdir))+': len of gp.files.outdir')
    pymultinest.run(myloglike,   myprior,
                    n_dims,      n_params = n_dims+1, # None beforehands
                    n_clustering_params = n_dims, # separate modes on
                                                  #the rho parameters
                                                  #only (gp.nepol in
                                                  #this case)
                    wrapped_params = [ gp.pops, gp.nipol, gp.nepol], # do
                                                                     #not
                                                                     #wrap-around
                                                                     #parameters
                    importance_nested_sampling = True, # INS enabled
                    multimodal = True,            # separate modes
                    const_efficiency_mode = True, # use const sampling efficiency
                    n_live_points = gp.nlive,
                    evidence_tolerance = 0.0, # set to 0 to keep
                                              #algorithm working
                                              #indefinitely
                    sampling_efficiency = 0.03, # very low eff. in case of const efficiency mode, README
                    n_iter_before_update = 1, # output after this many iterations
                    null_log_evidence = 1., # separate modes if
                                            #logevidence > this param.
                    max_modes = gp.nlive,   # preallocation of modes:
                                            #max. = number of live
                                            #points
                    mode_tolerance = -1.e30,   # mode tolerance in the
                                               #case where no special
                                               #value exists: highly
                                               #negative
                    outputfiles_basename = gp.files.outdir,
                    seed = -1,
                    verbose = True,
                    resume = False,
                    context = 0,
                    write_output = True,
                    log_zero = -1e6,      # points with log likelihood
                                          #< log_zero will be
                                          #neglected
                    max_iter = 0,         # set to 0 for never
                                          #reaching max_iter (no
                                          #stopping criterium based on
                                          #number of iterations)
                    init_MPI = True,      # use MPI (TODO: debug)
                    dump_callback = None)
    # progress.stop()    # ok, done. Stop our progress watcher

    # exit(0)

    # # store names of parameters, always useful
    # f = open('%sparams.json' % a.outputfiles_basename, 'w')
    # import json
    # json.dump(parameters, f, indent=2)
    # f.close()

    # # lets analyse the results
    # a = pymultinest.Analyzer(n_params = n_dims)
    # s = a.get_stats()
    
    # # store derived stats
    # f = open('%sstats.json' % a.outputfiles_basename, 'w')
    # json.dump(s, f, indent=2)
    # f.close()
    
    # print();     print('-' * 30, 'ANALYSIS', '-' * 30)
    # print('Global Evidence:\n\t%.15e +- %.15e'\
    #       % ( s['nested sampling global log-evidence'],\
    #           s['nested sampling global log-evidence error'] ))


if __name__=="__main__":
    import gl_params
    gp = gl_params.Params()
    
    from gl_class_cube import Cube
    from gl_loglike import geom_loglike
            

    run(gp)
