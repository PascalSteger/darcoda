#!/usr/bin/env python3

## @file
# split populations in observed dwarf galaxies
# based on metallicity, and Mg index, but no v_LOS nor position
# convention:       1, 2 is for first, second component

# (c) 2014 Pascal S.P. Steger

import pdb
import numpy as np
#import matplotlib as plt
from pylab import *
import pymultinest
import gl_helper as gh
from gl_centering import com_shrinkcircle_v_2D
import gl_units as gu

gh.DEBUGLEVEL = 3
DEBUG = True

def myprior(cube, ndim, nparams):
    # convert to physical space
    off = 0
    cube[off] = cube[off] # fraction of particles in part 1
    off +=1
    for pop in range(0, gp.pops+1): # no. of pops goes in here, first MW, then 1,2,..
        cube[off] = cube[off] # Fe_mu
        off += 1
        cube[off] = cube[off] # Fe_sig
        off += 1
        cube[off] = cube[off] # Mg_mu
        off += 1
        cube[off] = cube[off] # Mg_sig
        off += 1
    if off != ndim:
        gh.LOG(1, 'wrong number of parameters in myprior.cube')
        pdb.set_trace()
    return
## \fn myprior(cube, ndim, nparams) priors
# @param cube [0,1]^ndim cube, array of dimension ndim
# @param ndim number of dimensions, 2*npop*nipol + nipol
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
    Fe_mu = []; Fe_sig = []; Mg_mu = []; Mg_sig = []
    frac = cube[off]
    off += 1
    for pop in range(gp.pops):
        Fe_mu.append(cube[off])
        off += 1
        Fe_sig.append(cube[off])
        off += 1
        Mg_mu.append(cube[off])
        off += 1
        Mg_sig.append(cube[off])
        off += 1
    if off != ndim:
        gh.LOG(1, 'wrong number of parameters in myloglike.cube')
        pdb.set_trace()
    gh.LOG(2,'starting logev evaluation')
    p1_Fe= 1/np.sqrt(2*np.pi*(Fe_sig[0]**2+Fe_err**2))*np.exp(-(Fe-Fe_mu[0])**2/(2*np.sqrt(Fe_sig[0]**2+Fe_err**2)))
    p1_Mg= 1/np.sqrt(2*np.pi*(Mg_sig[0]**2+Mg_err**2))*np.exp(-(Mg-Mg_mu[0])**2/(2*np.sqrt(Mg_sig[0]**2+Mg_err**2)))
    p2_Fe= 1/np.sqrt(2*np.pi*(Fe_sig[1]**2+Fe_err**2))*np.exp(-(Fe-Fe_mu[1])**2/(2*np.sqrt(Fe_sig[1]**2+Fe_err**2)))
    p2_Mg= 1/np.sqrt(2*np.pi*(Mg_sig[1]**2+Mg_err**2))*np.exp(-(Mg-Mg_mu[1])**2/(2*np.sqrt(Mg_sig[1]**2+Mg_err**2)))
    logev = np.sum(np.log(frac*PM*p1_Fe*p1_Mg + (1-frac)*PM*p2_Fe*p2_Mg))
    gh.LOG(1, 'logL:',logev)
    return logev
## \fn myloglike(cube, ndim, nparams) calculate probability function
# @param cube ndim cube of physical parameter space (nr)
# @param ndim number of dimensions, 2*npop*nipol + nipol
# @param nparams = ndim + additional parameters
# stored with actual parameters


def run(gp):
    import gr_params
    gpr = gr_params.Params(gp)
    global DL
    DL = {0: lambda x: x * (138),#+/- 8 for Fornax
          1: lambda x: x * (101),#+/- 5 for Carina
          2: lambda x: x * (79), #+/- 4 for Sculptor
          3: lambda x: x * (86)  #+/- 4 for Sextans
      }[gp.case](gu.kpc__pc)
    global k2 # overall half-light radius in [pc] from Irwin,Hatzidimitriou1995
    k2 = {0: lambda x: x * (339),#+/-36 for Fornax
          1: lambda x: x * (137),#+/-22 for Carina
          2: lambda x: x * (94), #+/-26 for Sculptor
          3: lambda x: x * (294) #+/-38 for Sextans
      }[gp.case](1)

    global wpt, Rpt, V0, Ve0, Mg0, Mge0, PM0
    gpr.fil = gpr.dir+"/table_merged.bin"
    delim = [0,22,3,3,6,4,3,5,6,6,7,5,6,5,6,5,6]
    #ID = np.genfromtxt(gpr.fil, skiprows=29, unpack=True,\
    # usecols=(0,1),delimiter=delim)
    RAh,RAm,RAs,DEd,DEm,DEs,Vmag,VI,\
      VHel,e_VHel,Fe,Fe_err,\
      Mg,Mg_err,PM = np.genfromtxt(gpr.fil, skiprows=29, unpack=True, \
                                       usecols=tuple(range(2,17)), \
                                       delimiter=delim, filling_values=-1)
    # attention, we do not have Mg measurements for 501 stars in Fornax,
    #  visible by missing SigMg values, set to -1
    #   we exclude them from all further analysis
    global Nsample
    Nsample = len(PM)

    # where no Fe or Mg is measured, set the corresponding error to infinity
    for i in range(Nsample):
        if Fe[i] < 0:
            Fe_err[i] = np.inf
        elif Mg[i] < 0:
            Mg_err[i] = np.inf

    # use all stellar tracer particles from now on, independent on their probability of membership
    scatter(Fe, Mg)
    xlabel('Fe')
    ylabel('Mg')
    show()
    pdb.set_trace()

    global alpha_s, delta_s
    sig = abs(RAh[0])/RAh[0]
    RAh = RAh/sig
    # stellar position alpha_s, delta_s
    # 15degrees in 1 hour right ascension
    alpha_s = 15*(RAh*3600+RAm*60+RAs)*sig       # [arcsec]
    sig = abs(DEd[0])/DEd[0]                # +/-
    DEd = DEd/sig
    delta_s = (DEd*3600+DEm*60+DEs)*sig          # [arcsec]
    # unit conversion into a set of [pc], [km/s]
    #arcsec = 2.*np.pi/(360.*60.*60)      # [rad/arcsec]
    alpha_s /= gu.rad__arcsec  # [rad]
    delta_s /= gu.rad__arcsec  # [rad]

    # instead of using other datasets for individual dSph,
    # determine Heliocentric-Rest-Frame Line-Of-Sight velocity of dwarf Vd
    # and position of dwarf (ad, dd)
    # from probability-of-membership weighted means directly
    global Vd, ad, dd
    Vd = np.sum(VHel * PM)/np.sum(PM)    # [km/s]
    ad = np.sum(alpha_s * PM)/np.sum(PM) # [arcsec]
    dd = np.sum(delta_s * PM)/np.sum(PM) # [arcsec]
    # determine distance to dwarf TODO reference
    xs = alpha_s*DL # [pc]
    ys = delta_s*DL # [pc]
    PM0 = 1.*np.copy(PM)
    x0 = 1.*np.copy(xs)
    y0 = 1.*np.copy(ys) # [pc]

    # remove center displacement, already change x0
    com_x, com_y, com_vz = com_shrinkcircle_v_2D(x0, y0, VHel, PM) # [pc], [km/s]
    global R0
    R0 = np.sqrt(x0**2+y0**2)
    # sort by R0 so integral makes sense later
    order = np.argsort(R0)
    R0 = R0[order]; PM0 = PM0[order]
    x0 = x0[order]; y0  = y0[order]
    xs = xs[order]; ys  = ys[order]
    alpha_s = alpha_s[order]; delta_s = delta_s[order]
    V0 = V0[order]; Ve0 = Ve0[order]
    Fe = Fe[order]; Fe_err = Fe_err[order]
    Mg = Mg[order]; Mg_err = Mg_err[order]

    A = np.loadtxt(gp.files.dir+'w_2.0.dat')
    Rpt, wpt = A.T # [arcmin], [1]
    arcmin__pc = 2.*np.pi* DL / (360 * 60) # [pc] at distance 138 kpc for Fornax
    Rpt *= arcmin__pc # [pc]
    global gw
    gw = w(R0)

    gh.LOG(1,'starting MultiNest run:')
    n_dims = 1+gp.pops*4
    pymultinest.run(myloglike,   myprior,
                    n_dims,      n_params = n_dims+1, # None beforehands
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
                    n_live_points = gp.nlive,
                    evidence_tolerance = 0.0,   # 0 to keep working infinitely
                    sampling_efficiency = 0.05, # very low eff. in
                                                #case of const efficiency mode,
                                                #README
                    n_iter_before_update = 10, # output after this many iterations
                    null_log_evidence = 1., # separate modes if
                                            #logevidence > this param.
                    max_modes = gp.nlive,
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


if __name__=='__main__':
    import gl_params
    gp = gl_params.Params()

    run(gp)

# works with investigation = 'obs', pops = 2, metalpop = True
# profile with python3 -m cProfile grd_split.py
### TODO output:
