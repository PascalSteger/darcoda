#!/usr/bin/env python3

## @file
# split populations in observed dwarf galaxies
# based on metallicity, and Mg index, but no v_LOS nor position
# convention:       1, 2 is for first, second component

# (c) 2014 Pascal S.P. Steger

import pdb
import numpy as np
import multinest
import gl_helper as gh

gh.DEBUGLEVEL = 1
DEBUG = True

def myprior(cube, ndim, nparams):
    # convert to physical space
    off = 0
    cube[off] = cube[off] # fraction of particles in part 1
    off +=1
    for pop in range(gp.pops): # no. of pops goes in here, first MW, then 1,2,..
        # TODO: add up differences from 0
        cube[off] = cube[off]*(Mg_max-Mg_min)+Mg_min # Mg_mu
        off += 1
        cube[off] = cube[off]*(Mg_max-Mg_min) # Mg_sig
        off += 1
    if off != ndim:
        gh.LOG(1, 'wrong number of parameters in myprior.cube')
        pdb.set_trace()
    return cube
## \fn myprior(cube, ndim, nparams) priors
# @param cube [0,1]^ndim cube, array of dimension ndim
# @param ndim number of dimensions
# @param nparams = ndim + additional parameters
# stored with actual parameters

def w(Rk):
    gh.sanitize_vector(Rk, Nsample, 0, 1e30, DEBUG)
    w_ipol = np.zeros(Nsample)
    for k in range(Nsample):
        w_ipol[k] = wpt[np.where(abs(Rk[k]-Rpt) == min(abs(Rk[k]-Rpt)))]
    return w_ipol
## \fn w(Rk)
# return selection function as function of radius
# take vector as input
# @param Rk radius [pc]

def myloglike(cube, ndim, nparams):
    off = 0
    Mg_mu = []; Mg_sig = []
    frac = cube[off]
    off += 1
    for pop in range(gp.pops):
        Mg_mu.append(cube[off])
        off += 1
        Mg_sig.append(cube[off])
        off += 1
    if off != ndim:
        gh.LOG(1, 'wrong number of parameters in myloglike.cube')
        pdb.set_trace()
    gh.LOG(2,'starting logev evaluation')
    p1_Mg= 1/np.sqrt(2*np.pi*(Mg_sig[0]**2+Mg_err**2))*np.exp(-(Mg-Mg_mu[0])**2/(2*np.sqrt(Mg_sig[0]**2+Mg_err**2)))
    p2_Mg= 1/np.sqrt(2*np.pi*(Mg_sig[1]**2+Mg_err**2))*np.exp(-(Mg-Mg_mu[1])**2/(2*np.sqrt(Mg_sig[1]**2+Mg_err**2)))
    p1 = frac*PM*p1_Mg
    p2 = (1-frac)*PM*p2_Mg
    pcom = p1+p2
    print('pcom (min, max) = ', min(pcom), max(pcom))
    print('fraction of pcom == 0 : ', sum(pcom==0)/len(pcom))
    lpcom = np.log(pcom)
    logev = np.sum(lpcom)
    gh.LOG(1, 'logL:',logev)
    #if logev < -1e300:
    #    pdb.set_trace()
    return logev
## \fn myloglike(cube, ndim, nparams) calculate probability function
# @param cube ndim cube of physical parameter space (nr)
# @param ndim number of dimensions, 2*npop*nipol + nipol
# @param nparams = ndim + additional parameters
# stored with actual parameters

def run(gp):
    n_dims = 1+gp.pops*2
    multinest.run(myloglike,   myprior,
                    n_dims,      n_params = n_dims,
                    n_clustering_params = n_dims, # separate modes on
                                                  # the rho parameters
                                                  # only (gp.nrho in
                                                  # this case)
                    wrapped_params = [ gp.pops, gp.nipol, gp.nrho], # do
                                                                     #not
                                                                     #wrap-around
                                                                     #parameters
                    importance_nested_sampling = True, # INS enabled
                    multimodal = True,            # separate modes
                    const_efficiency_mode = True, # use const sampling efficiency
                    n_live_points = Nsample,
                    evidence_tolerance = 0.0,   # 0 to keep working infinitely
                    sampling_efficiency = 0.95,
                    n_iter_before_update = 100, # output after this many iterations
                    null_log_evidence = 1., # separate modes if
                                            #logevidence > this param.
                    max_modes = Nsample,
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
                    log_zero = -999999, # points with log L < log_zero will be
                                          # neglected
                    max_iter = 0,
                    init_MPI = True,     # use MPI
                    dump_callback = None)

if __name__=="__main__":
    from mpi4py import MPI
    import numpy

    comm = MPI.Comm.Get_parent()
    size = comm.Get_size()
    rank = comm.Get_rank()

    import gl_params
    gp = gl_params.Params()

    globs = None
    globs = comm.bcast(globs, root=0)
    global Nsample, wpt, Rpt, Mg, Mg_err, PM, Mg_min, Mg_max
    gp, Nsample, wpt, Rpt, Mg, Mg_err, PM, Mg_min, Mg_max = globs
    run(gp)
    finish = 1
    comm.Reduce(finish, None,
                op=MPI.SUM, root=0)
    comm.Disconnect()