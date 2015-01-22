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


#import matplotlib
#matplotlib.use('pdf')
#from pylab import *
#ion()

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

def introduce_points_in_between(r0, gp):
    rmin = np.log10(min(r0))
    rmax = np.log10(max(r0))
    return np.logspace(rmin, rmax, gp.nfine)
## \fn introduce_points_in_between(r0, gp)
# get gp.fine points logarithmically spaced points
# @param r0 [pc] gp.xipol
# @param gp global parameter

def myprior(cube, ndim, nparams):
    # convert to physical space
    off = 0
    cube[off] = cube[off]*0.8+0.1 # fraction of particles in part 1, with min 0.1, max 0.9
    # such that each population has at least 10% of the total no. stars
    off +=1
    for pop in range(2): # 2 pops
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
    for pop in range(2):
        Mg_mu.append(cube[off])
        off += 1
        Mg_sig.append(cube[off])
        off += 1
    gh.sanitize_vector(Mg_mu, 2, -10, 10, True)
    if off != ndim:
        gh.LOG(1, 'wrong number of parameters in myloglike.cube')
        pdb.set_trace()
    gh.LOG(2,'starting logev evaluation')
    p1_Mg= 1/np.sqrt(2*np.pi*(Mg_sig[0]**2+Mg_err**2))*\
           np.exp(-(Mg-Mg_mu[0])**2/(2*(Mg_sig[0]**2+Mg_err**2)))
    p2_Mg= 1/np.sqrt(2*np.pi*(Mg_sig[1]**2+Mg_err**2))*\
           np.exp(-(Mg-Mg_mu[1])**2/(2*(Mg_sig[1]**2+Mg_err**2)))
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
    import gr_params
    gpr = gr_params.grParams(gp)

    #globs = comm.bcast(globs, root=0)
    global Nsample, Mg, Mg_err, PM, Mg_min, Mg_max
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


    sel = (Mg>-1)  # exclude missing data on Mg
    RAh = RAh[sel]
    RAm = RAm[sel]
    RAs = RAs[sel]
    DEd = DEd[sel]
    DEm = DEm[sel]
    DEs = DEs[sel]
    Vmag = Vmag[sel]
    VI  = VI[sel]
    VHel = VHel[sel]
    e_VHel = e_VHel[sel]
    Mg = Mg[sel]
    Mg_err = Mg_err[sel]
    PM = PM[sel]

    Mg_min = min(Mg) # -3, 3 if according to WalkerPenarrubia2011
    Mg_max = max(Mg)

    # easiest way for visualization: use histogram to show data
    #hist(Mg, np.sqrt(len(Mg))/2, normed=True)

    # but: it's not as easy as that
    # we have datapoints with errors and probability of membership weighting
    # thus, we need to smear the values out using a Gaussian of width = Mg_err
    # and add them up afterwards after scaling with probability PM
    x = np.array(np.linspace(min(Mg), max(Mg), 100))
    Mgdf = np.zeros(100)
    for i in range(len(Mg)):
        Mgdf += PM[i]*gh.gauss(x, Mg[i], Mg_err[i])
    Mgdf /= sum(PM)

    #plot(x, Mgdf, 'g', lw=2)
    # only then we want to compare to Gaussians

    n_dims = 1+gp.pops*2
    #Nsample = 10*n_dims
    pymultinest.run(myloglike,
                  myprior,
                  n_dims, # nest_ndims
                  n_dims+1, # nest_totPar
                  n_dims, # separate modes on nest_nCdims
                  # the rho parameters only (gp.nrho in this case)
                  [ gp.pops, gp.nipol, gp.nrho],
                  True, # nest_IS = INS enabled
                  True, #nest_mmodal =            # separate modes
                  True, # nest_ceff = use const sampling efficiency
                  Nsample, # nest_nlive =
                  0.0,   # nest_tol = 0 to keep working infinitely
                  0.8, # nest_ef =
                  10000, # nest_updInt = output after this many iterations
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
                  1000, # nest_maxIter =
                  True,     # initMPI =  use MPI
                  None) #dump_callback =


    #params after xmas
    # TODO: get most likely parameters from stat.dat
    cubeML= np.array([0.475319236624166197E+00, 0.621662395675444568E+00, 0.798401723057411417E-01, 0.550211197376269112E+00, 0.158468782949331616E+00])

    cubeMLphys=myprior(cubeML, 1+gp.pops*2, 1+gp.pops*2)
    #myloglike(cubeMLphys, 1+gp.pops*2, 1+gp.pops*2)
    pML, mu1ML, sig1ML, mu2ML, sig2ML = cubeMLphys
    #pdb.set_trace()
    g1 = pML*gh.gauss(x, mu1ML, sig1ML)
    g2 = (1-pML)*gh.gauss(x, mu2ML, sig2ML)
    gtot = g1+g2

    #plot(x, pML*g1, 'white')
    #plot(x, (1-pML)*g2, 'white')
    #plot(x, gtot, 'r')
    #xlabel('Mg')
    #ylabel('pdf')

    # get center of photometric measurements by deBoer
    # for Fornax, we have
    com_x = 96203.736358393697
    com_y = -83114.080684733024
        # only use stars which are members of the dwarf
    sig = abs(RAh[0])/RAh[0]
    RAh = RAh/sig
    xs = 15*(RAh*3600+RAm*60+RAs)*sig       # [arcsec/15]
    sig = abs(DEd[0])/DEd[0]
    DEd = DEd/sig
    ys = (DEd*3600+DEm*60+DEs)*sig          # [arcsec]
    arcsec = 2.*np.pi/(360.*60.*60) # [pc]
    kpc = 1000 # [pc]
    DL = {0: lambda x: x * (138),#+/- 8 for Fornax
          1: lambda x: x * (101),#+/- 5 for Carina
          2: lambda x: x * (79),  #+/- 4 for Sculptor
          3: lambda x: x * (86) #+/- 4 for Sextans
      }[gp.case](kpc)

    xs *= (arcsec*DL) # [pc]
    ys *= (arcsec*DL) # [pc]


    # instantiate different samplings, store half-light radii (2D)
    R1half = []
    R2half = []

    # get a sample assignment:
    import numpy.random as npr
    popass = []
    for i in range(sum(sel)):
        mg = Mg[i]
        ppop1 = pML*gh.gauss(mg, mu1ML, sig1ML)
        ppop2 = (1-pML)*gh.gauss(mg, mu2ML, sig2ML)
        if npr.rand() <= ppop1/(ppop1+ppop2):
            popass.append(1)
        else:
            popass.append(2)

    popass = np.array(popass)
    sel1 = (popass==1)
    sel2 = (popass==2)
    #radii of all stellar tracers from pop 1 and 2
    R1 = np.sqrt((xs[sel1]-com_x)**2 + (ys[sel1]-com_y)**2)
    R2 = np.sqrt((xs[sel2]-com_x)**2 + (ys[sel2]-com_y)**2)
    R1.sort()
    R2.sort()

    import gi_project as gip
    for pop in np.arange(2)+1:
        if pop == 1:
            R0 = R1 # [pc]
            Rhalf = R1[len(R1)/2]
            R1half.append(Rhalf)
            co = 'blue'
        else:
            R0 = R2 # [pc]
            Rhalf = R2[len(R2)/2]
            R2half.append(Rhalf)
            co = 'red'
        Rdiff = np.abs(R1[len(R1)/2] - R2[len(R2)/2])
        if Rdiff < 95 or Rdiff > 105:
            continue

        # calculate 2D Sig profiles
        Rmin = min(R0) # [pc]
        Rmax = max(R0) # [pc]
        Binmin, Binmax, Rbin = gh.determine_radius(R0, Rmin, Rmax, gp) # [pc]
        gp.xipol = Rbin # [pc]
        minr = min(Rbin)                           # [pc]
        maxr = max(Rbin)                           # [pc]
        Vol = gh.volume_circular_ring(Binmin, Binmax, gp) # [pc^2]
        totmass_tracers = float(len(x))
        Rsi   = gh.add_errors(R0,   gpr.Rerr)   # [pc], gpr.Rerr was in
        tpb = np.zeros(gp.nipol)
        Sig_phot = np.zeros(gp.nipol)
        for i in range(gp.nipol):
            ind1 = np.argwhere(np.logical_and(Rsi >= Binmin[i], \
                                              Rsi <  Binmax[i])).flatten() # [1]
            tpb[i] = float(len(ind1)) # [1]
            Sig_phot[i] = float(len(ind1))*totmass_tracers/Vol[i] # [Munit/pc^2]

        # loglog(gp.xipol, Sig_phot, co)
        # axvline(Rhalf, color=co)
        # xlim([min(gp.xipol), max(gp.xipol)])
        # xlabel(r'$R$')
        # ylabel(r'$\Sigma(R)$')

        # deproject to get 3D nu profiles
        gp.xipol = Rbin
        minr = min(Rbin)                           # [pc]
        maxr = max(Rbin)                           # [pc]
        gp.xepol = np.hstack([minr/8., minr/4., minr/2.,\
                              Rbin, \
                              2*maxr, 4*maxr, 8*maxr]) # [pc]
        gp.xfine = introduce_points_in_between(gp.xepol, gp)
        Sigdatnu, Sigerrnu = gh.complete_nu(Rbin, Sig_phot, Sig_phot/10., gp.xfine)
        dummyx, nudatnu, nuerrnu, Mrnu = gip.Sig_NORM_rho(gp.xfine, \
                                                          Sigdatnu, Sigerrnu,\
                                                          gp)
        #nudat = gh.linipollog(gp.xfine, nudatnu, gp.xipol)
        #nuerr = gh.linipollog(gp.xfine, nuerrnu, gp.xipol)
        #loglog(gp.xipol, nudat, co)
        #axvline(Rhalf, color=co)
        #xlim([min(gp.xipol), max(gp.xipol)])
        #xlabel(r'$R$')
        #ylabel(r'$\nu(R)$')
        #plum = 100*gh.plummer(gp.xipol, Rhalf, len(R0))
        #loglog(gp.xipol, plum, color=co, linestyle='--')
        #ylim([min(plum), max(plum)])
        #pdb.set_trace()

    R1half = np.array(R1half)
    R2half = np.array(R2half)
    Rdiffhalf = np.abs(R1half-R2half)

    #clf()
    #hist(Rdiffhalf, np.sqrt(len(Rdiffhalf))/2)
    #xlabel(r'\Delta R')
    #ylabel('pdf')
    #pdb.set_trace()
    np.savetxt(gp.files.dir+'popass', popass)
    return
## \fn run(gp)
# run MultiNest
# @param gp global parameters defined in gi_params.py

if __name__=="__main__":
    import gi_params
    global gp
    gp = gi_params.Params()
    if gp.pops < 2:
        gh.LOG(1, " population splitting needs 2 or more populations, corrected")
        gp.pops = 2
    # convention: directory names have ending "/"
    run(gp)
