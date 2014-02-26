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

import gl_priors as gprio
import gl_chi as gc
import gl_physics as phys
from gl_class_cube import Cube
from gl_class_profiles import Profiles
from gl_project import rho_SUM_Mr
import gl_helper as gh

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
    tmp_profs = Profiles(gp.pops, gp.nipol)
    off = 0

    rho_param = np.array(cube[off:off+gp.nepol])
    # TODO: enable bound on n(r) wiggling
    # better approach: use other representation of n(r) from [0,1] priors
    # if gprio.check_nr(rho_param[2:-1]):
    #     print('dn/dr too big!')
    #     return gh.err(0.7, gp)

    tmp_rho = phys.rho(gp.xepol, rho_param, gp)
    if(gprio.check_rho(tmp_rho, gp.rprior, gp.rhotol)):
        print('rho slope error') # should be ensured by [0,1] priors and maxslope
        return gh.err(1.)
    tmp_profs.set_rho(tmp_rho[:gp.nipol])
    off += gp.nepol

    nuparstore = []
    for pop in np.arange(gp.pops)+1:
        nu_param = cube[off:off+gp.nepol]
        nuparstore.append(nu_param)
        
        tmp_nu = phys.rho(gp.xepol, nu_param, gp) #  [1], [pc]
        # if gprio.check_nu(tmp_nu):
        #     print('nu error')
        #     return err/2.
        # TODO: turn on baryonic prior check
        # if gp.bprior and gprio.check_bprior(tmp_rho, tmp_nu):
        #      print('bprior error')
        #      return gh.err(1.5, gp)
        tmp_profs.set_nu(pop, tmp_nu[:gp.nipol]) # [munit/pc^3]
        off += gp.nepol

        beta_param = np.array(cube[off:off+gp.nbeta])
        tmp_beta = phys.beta(gp.xipol, beta_param, gp)
        if gprio.check_beta(tmp_beta, gp):
            print('beta error')
            return gh.err(2., gp)
        tmp_profs.set_beta(pop, tmp_beta)
        off += gp.nbeta

        try:
            sig, kap = phys.sig_kap_los(gp.xepol, pop, rho_param, nu_param, beta_param, gp)
            # sig and kap already are on data radii only, so no extension by 3 bins here
        except Exception as detail:
            return gh.err(3., gp)
        tmp_profs.set_sig_kap(pop, sig, kap)
    
    # determine log likelihood (*not* reduced chi2)
    chi2 = gc.calc_chi2(tmp_profs, nuparstore, gp)
    print('found log likelihood = ', -chi2/2.)
    return -chi2/2.   # from   likelihood L = exp(-\chi^2/2), want log of that

## \fn myloglike(cube, ndim, nparams)
# calculate probability function
# @param cube [0,1]^ndim cube
# @param ndim number of dimensions, 2*npop*nipol + nipol
# @param nparams = ndim + additional parameters stored with actual parameters


def stringlist(pops, nipol):
    tmp = []
    for i in range(nipol):
        tmp.append('rho_%d' % i)
    for j in range(pops):
        for i in range(nipol):
            tmp.append('nu_%d,%d' % (j,i))
        for i in range(gp.nbeta):
            tmp.append('betastar_%d,%d' % (j,i))
    return tmp
## \fn stringlist(pops, nipol)
# show parameter names in array
# @param pops int, number of populations
# @param nipol int, number of radial bins


def run(gp):
    import gl_file as gfile
    if gp.getnewdata:
        gfile.bin_data(gp)
    gfile.get_data(gp)
    
    ## number of dimensions
    n_dims = gp.nepol + gp.pops*gp.nepol + gp.pops*gp.nbeta #rho, (nu, beta)_i
    parameters = stringlist(gp.pops, gp.nepol)
    
    # show live progress
    # progress = pymultinest.ProgressPlotter(n_params = n_dims)
    # progress.start()
    # threading.Timer(2, show, [gp.files.outdir+'/phys_live.points.pdf']).start() 

    # print(str(len(gp.files.outdir))+': len of gp.files.outdir')
    pymultinest.run(myloglike,   myprior,
                    n_dims,      n_params = n_dims+1, # None beforehands
                    n_clustering_params = gp.nepol, # separate modes on the rho parameters only (gp.nepol in this case)
                    wrapped_params = [ gp.pops, gp.nipol, gp.nepol],          # do not wrap-around parameters
                    importance_nested_sampling = False, # INS enabled
                    multimodal = False,  # separate modes
                    const_efficiency_mode = True, # use const sampling efficiency
                    n_live_points = gp.nlive,
                    evidence_tolerance = 0.0, # set to 0 to keep algorithm working indefinitely
                    sampling_efficiency = 0.5,
                    n_iter_before_update = gp.nlive, # output after this many iterations
                    null_log_evidence = 1e30, # separate modes if logevidence > this param.
                    max_modes = gp.nlive,   # preallocation of modes: max. = number of live points
                    mode_tolerance = -1.,
                    outputfiles_basename = gp.files.outdir,
                    seed = -1,
                    verbose = True,
                    resume = False,
                    context = 0,
                    write_output = True,
                    log_zero = -1e6,
                    max_iter = -1,
                    init_MPI = False,
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
    
    gp.files.makedir(gp)
    run(gp)

