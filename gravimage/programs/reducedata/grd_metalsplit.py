#!/usr/bin/env ipython3

## @file
# split populations in observed dwarf galaxies
# based on metallicity, and Mg index, but no v_LOS nor position
# convention:       1, 2 is for first, second component

# (c) 2014 Pascal S.P. Steger, pascal@steger.aero

import pdb
import sys
import numpy as np
import pymultinest
import gi_helper as gh

gh.DEBUGLEVEL = 1
DEBUG = True

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

def bufcount(filename):
    f = open(filename)
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read # loop optimization
    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)
    return lines
## \fn bufcount(filename)
# determine no. lines optimally
# @param filename filename

def myprior(cube, ndim, nparams):
    # convert to physical space
    off = 0
    cube[off] = cube[off]*0.8+0.1 # fraction of particles in part 1, with min 0.1, max 0.9
    # such that each population has at least 10% of the total no. stars
    off +=1
    for pop in range(gp.pops): # no. of pops goes in here, first MW, then 1,2,..
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

# def w(Rk):
#     gh.sanitize_vector(Rk, Nsample, 0, 1e30, DEBUG)
#     w_ipol = np.zeros(Nsample)
#     for k in range(Nsample):
#         w_ipol[k] = wpt[np.where(abs(Rk[k]-Rpt) == min(abs(Rk[k]-Rpt)))]
#     return w_ipol
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
    gh.sanitize_vector(Mg_mu, 2, -10, 10, True)
    if off != ndim:
        gh.LOG(1, 'wrong number of parameters in myloglike.cube')
        pdb.set_trace()
    gh.LOG(2,'starting logev evaluation')
    p1_Mg= 1/np.sqrt(2*np.pi*(Mg_sig[0]**2+Mg_err**2))*np.exp(-(Mg-Mg_mu[0])**2/(2*np.sqrt(Mg_sig[0]**2+Mg_err**2)))
    p2_Mg= 1/np.sqrt(2*np.pi*(Mg_sig[1]**2+Mg_err**2))*np.exp(-(Mg-Mg_mu[1])**2/(2*np.sqrt(Mg_sig[1]**2+Mg_err**2)))
    p1 = frac*PM*p1_Mg
    for i in range(0,len(p1)):
        if p1[i] == 0.0:
            p1[i] = 1e-30
    p2 = (1-frac)*PM*p2_Mg
    for i in range(0, len(p2)):
        if p2[i] == 0.0:
            p2[i] = 1e-30
    pcom = p1+p2
    #print('pcom (min, max) = ', min(pcom), max(pcom))
    #print('fraction of pcom == 0 : ', sum(pcom==0)/len(pcom))
    lpcom = np.log(pcom)
    logev = np.sum(lpcom)
    #print(logev)
    #gh.LOG(1, 'logL:',logev)
    if logev < -1e300:
        logev = -1e300
    #    pdb.set_trace()
    return logev
## \fn myloglike(cube, ndim, nparams) calculate probability function
# @param cube ndim cube of physical parameter space (nr)
# @param ndim number of dimensions, 2*npop*nipol + nipol
# @param nparams = ndim + additional parameters
# stored with actual parameters

def run(gp):
    n_dims = 1+gp.pops*2
    pymultinest.run(myloglike,
                  myprior,
                  n_dims, # nest_ndims
                  n_dims, # nest_totPar
                  n_dims, # separate modes on nest_nCdims
                  # the rho parameters only (gp.nrho in this case)
                  [ gp.pops, gp.nipol, gp.nrho], # do
                  #not wrap-around parameters wrapped_params =
                  True, # nest_IS = INS enabled
                  True, #nest_mmodal =            # separate modes
                  True, # nest_ceff = use const sampling efficiency
                  Nsample, # nest_nlive =
                  0.0,   # nest_tol = 0 to keep working infinitely
                  0.25, # nest_ef =
                  4000, # nest_updInt = output after this many iterations
                  1., # null_log_evidence separate modes if
                  #logevidence > this param.
                  Nsample, # maxClst =
                  -1.e30,   # nest_Ztol = mode tolerance in the
                  #case where no special value exists: highly negative
                  gp.files.outdir, # outputfiles_basename =
                  -1, # seed =
                  True, # nest_fb =
                  False, # nest_resume =
                  0, # context =
                  True, # nest_outfile =
                  -999999, # nest_logZero = points with log L < log_zero will be
                  # neglected
                  0, # nest_maxIter =
                  True,     # initMPI =  use MPI
                  None) #dump_callback =

## \fn run(gp)
# run MultiNest
# @param gp global parameters defined in gi_params.py


if __name__=="__main__":
#    from mpi4py import MPI
#    comm = MPI.Comm.Get_parent()
#    size = comm.Get_size()
#    rank = comm.Get_rank()

    import gi_params
    global gp
    gp = gi_params.Params()
    if gp.pops < 2:
        gh.LOG(1, " population splitting needs 2 or more populations, corrected")
        gp.pops = 2

    #globs = comm.bcast(globs, root=0)
    global Nsample, Mg, Mg_err, PM, Mg_min, Mg_max

    Mg_min = -3 # according to WalkerPenarrubia2011
    Mg_max = 3

    import gr_params
    gpr = gr_params.grParams(gp)
    # convention: directory names have ending "/"
    gpr.fil = gpr.dir+"table_merged.bin"

    # number of measured tracer stars
    Nsample = bufcount(gpr.fil)


    delim = [0,22,3,3,6,4,3,5,6,6,7,5,6,5,6,5,6]
    ID = np.genfromtxt(gpr.fil, skiprows=29, unpack=True,\
                       usecols=(0,1),delimiter=delim)

    RAh,RAm,RAs,DEd,DEm,DEs,Vmag,VI,\
      VHel,e_VHel,SigFe,e_SigFe,\
      Mg,Mg_err,PM = np.genfromtxt(gpr.fil, skiprows=29, unpack=True, \
                                   usecols=tuple(range(2,17)), delimiter=delim, filling_values=-1)

    run(gp)
    finish = 1
#    comm.Reduce(finish, None,
#                op=MPI.SUM, root=0)
#    comm.Disconnect()
