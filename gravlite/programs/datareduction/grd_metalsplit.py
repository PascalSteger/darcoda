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

gh.DEBUGLEVEL = 1
DEBUG = True

def myprior(cube, ndim, nparams):
    # convert to physical space
    off = 0
    cube[off] = cube[off] # fraction of particles in part 1
    off +=1
    for pop in range(gp.pops): # no. of pops goes in here, first MW, then 1,2,..
        # TODO: add up differences from 0
        #cube[off] = cube[off]*(Fe_max-Fe_min)+Fe_min # Fe_mu
        #off += 1
        #cube[off] = cube[off]*(Fe_max-Fe_min) # Fe_sig
        #off += 1
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
    Fe_mu = []; Fe_sig = []; Mg_mu = []; Mg_sig = []
    frac = cube[off]
    off += 1
    for pop in range(gp.pops):
        #Fe_mu.append(cube[off])
        #off += 1
        #Fe_sig.append(cube[off])
        #off += 1
        Mg_mu.append(cube[off])
        off += 1
        Mg_sig.append(cube[off])
        off += 1
    if off != ndim:
        gh.LOG(1, 'wrong number of parameters in myloglike.cube')
        pdb.set_trace()
    gh.LOG(2,'starting logev evaluation')
    #p1_Fe= 1/np.sqrt(2*np.pi*(Fe_sig[0]**2+Fe_err**2))*np.exp(-(Fe-Fe_mu[0])**2/(2*np.sqrt(Fe_sig[0]**2+Fe_err**2)))
    p1_Mg= 1/np.sqrt(2*np.pi*(Mg_sig[0]**2+Mg_err**2))*np.exp(-(Mg-Mg_mu[0])**2/(2*np.sqrt(Mg_sig[0]**2+Mg_err**2)))
    #p2_Fe= 1/np.sqrt(2*np.pi*(Fe_sig[1]**2+Fe_err**2))*np.exp(-(Fe-Fe_mu[1])**2/(2*np.sqrt(Fe_sig[1]**2+Fe_err**2)))
    p2_Mg= 1/np.sqrt(2*np.pi*(Mg_sig[1]**2+Mg_err**2))*np.exp(-(Mg-Mg_mu[1])**2/(2*np.sqrt(Mg_sig[1]**2+Mg_err**2)))
    p1 = frac*PM*p1_Mg # *p1_Fe
    p2 = (1-frac)*PM*p2_Mg # *p2_Fe
    pcom = p1+p2
    lpcom = np.log(pcom)
    logev = np.sum(lpcom)
    gh.LOG(1, 'logL:',logev)
    if logev < -1e300:
        pdb.set_trace()
    return logev
## \fn myloglike(cube, ndim, nparams) calculate probability function
# @param cube ndim cube of physical parameter space (nr)
# @param ndim number of dimensions, 2*npop*nipol + nipol
# @param nparams = ndim + additional parameters
# stored with actual parameters

def show_Fe_Mg_2D(Fe, Fe_err, Mg, Mg_err):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import kde
    np.random.seed(1977)

    # Generate 200 correlated x,y points
    data = np.vstack([Fe, Mg]).T
    #data = np.random.multivariate_normal([0, 0], [[1, 0.5], [0.5, 3]], 200)
    x, y = data.T
    nbins = 40

    fig, axes = plt.subplots(ncols=2, nrows=2, sharex=True, sharey=True)

    axes[0, 0].set_title('Scatterplot')
    axes[0, 0].plot(x, y, 'ko')
    axes[0, 0].set_aspect('equal')

    axes[0, 1].set_title('Hexbin plot')
    axes[0, 1].hexbin(x, y, gridsize=nbins)
    axes[0, 1].set_aspect('equal')

    axes[1, 0].set_title('2D Histogram')
    axes[1, 0].hist2d(x, y, bins=nbins)
    axes[1, 0].set_aspect('equal')

    # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
    k = kde.gaussian_kde(data.T)
    xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    axes[1, 1].set_title('Gaussian KDE')
    axes[1, 1].pcolormesh(xi, yi, zi.reshape(xi.shape))
    axes[1, 1].set_aspect('equal')

    fig.tight_layout()
    plt.show()
    return

    from matplotlib.patches import Ellipse
    NUM = len(Fe)
    ells = [Ellipse(xy=np.array([Fe[i], Mg[i]]), width=Fe_err[i], height=Mg_err[i], angle=0,\
                    color='blue', linewidth=0, alpha=0.01*PM[i])
            for i in range(NUM)]

    fig = figure()
    ax = fig.add_subplot(111, aspect='equal')
    for e in ells:
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)
    ax.set_xlabel('Fe [Ang]')
    ax.set_ylabel('Mg [Ang]')
    ax.set_xlim(-1, 1.5)
    ax.set_ylim(-1, 1.5)
    show()
## \fn show_Fe_Mg_2D(Fe, Fe_err, Mg, Mg_err)
# show ellipses with error bars for each star's Fe and Mg
# @param Fe iron abundance
# @param Fe_err error on it
# @param Mg Magnesium abundance
# @param Mg_err error on it

def show_Mg(Mg, Mg_err):
    hist(Mg, np.sqrt(Nsample))
    return
## \fn show_Mg(Mg, Mg_err)
# show histogram of Magnesium indices
# @param Mg vector, indices
# @prama Mg_err vector, corresponding errors

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

    global wpt, Rpt, V0, Ve0, Mg, Mg_err, Fe, Fe_err,PM
    gpr.fil = gpr.dir+"/table_merged.bin"
    delim = [0,22,3,3,6,4,3,5,6,6,7,5,6,5,6,5,6]
    #ID = np.genfromtxt(gpr.fil, skiprows=29, unpack=True,\
    # usecols=(0,1),delimiter=delim)
    RAh,RAm,RAs,DEd,DEm,DEs,Vmag,VI,\
      Vhel,Vhel_err,Fe,Fe_err,\
      Mg,Mg_err,PM = np.genfromtxt(gpr.fil, skiprows=29, unpack=True, \
                                       usecols=tuple(range(2,17)), \
                                       delimiter=delim, filling_values=-1)
    # exclude 0 probability of memory models, so we can weight by PM later on
    sel = (PM>0)
    RAh=RAh[sel]; RAm=RAm[sel]; RAs=RAs[sel]; DEd=DEd[sel]; DEm=DEm[sel]; DEs=DEs[sel]
    Vmag=Vmag[sel]; VI=VI[sel]; Vhel=Vhel[sel]; Vhel_err=Vhel_err[sel]
    Fe=Fe[sel]; Fe_err=Fe_err[sel]; Mg=Mg[sel]; Mg_err=Mg_err[sel]; PM=PM[sel]
    global Fe_min, Fe_max, Mg_min, Mg_max
    Fe_min = min(Fe); Fe_max = max(Fe); Mg_min = min(Mg); Mg_max = max(Mg)

    # attention, we do not have Mg measurements for 501 stars in Fornax,
    #  visible by missing SigMg values, set to -1
    #   we exclude them from all further analysis
    global Nsample
    Nsample = len(PM)

    # where no Fe or Mg is measured, set the corresponding error to infinity
    Fe_mean = np.mean(Fe[Fe>-1])
    Mg_mean = np.mean(Mg[Mg>-1])
    penalty_err = np.inf
    for i in range(Nsample):
        if Fe[i] == -1:
            Fe[i] = Fe_mean
            Fe_err[i] = penalty_err
        if Mg[i] == -1:
            Mg[i] = Mg_mean
            Mg_err[i] = penalty_err
        if Fe_err[i] == -1:
            Fe_err[i] = penalty_err
        if Mg_err[i] == -1:
            Mg_err[i] = penalty_err

    # use all stellar tracer particles from now on, independent on their probability of membership
    #scatter(Fe, Mg)
    #show_Fe_Mg_2D(Fe, Fe_err, Mg, Mg_err)
    # TODO show ellipses with errors, and alpha
    #xlabel('Fe')
    #ylabel('Mg')
    #show()
    #pdb.set_trace()

    show_Mg(Mg, Mg_err)
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
    Vd = np.sum(Vhel * PM)/np.sum(PM)    # [km/s]
    ad = np.sum(alpha_s * PM)/np.sum(PM) # [arcsec]
    dd = np.sum(delta_s * PM)/np.sum(PM) # [arcsec]
    # determine distance to dwarf TODO reference
    xs = alpha_s*DL # [pc]
    ys = delta_s*DL # [pc]
    PM0 = 1.*np.copy(PM)
    x0 = 1.*np.copy(xs)
    y0 = 1.*np.copy(ys) # [pc]

    # remove center displacement, already change x0
    com_x, com_y, com_vz = com_shrinkcircle_v_2D(x0, y0, Vhel, PM) # [pc], [km/s]
    global R0
    R0 = np.sqrt(x0**2+y0**2)
    # sort by R0 so integral makes sense later
    order = np.argsort(R0)
    R0 = R0[order]; PM0 = PM0[order]
    x0 = x0[order]; y0  = y0[order]
    xs = xs[order]; ys  = ys[order]
    alpha_s = alpha_s[order]; delta_s = delta_s[order]
    Vhel = Vhel[order]; Vhel_err = Vhel_err[order]
    Fe = Fe[order]; Fe_err = Fe_err[order]
    Mg = Mg[order]; Mg_err = Mg_err[order]

    A = np.loadtxt(gp.files.dir+'w_2.0.dat')
    Rpt, wpt = A.T # [arcmin], [1]
    arcmin__pc = 2.*np.pi* DL / (360 * 60) # [pc] at distance 138 kpc for Fornax
    Rpt *= arcmin__pc # [pc]

    gh.LOG(1,'starting MultiNest run:')
    n_dims = 1+gp.pops*2 # 4 if using Fe too
    pymultinest.run(myloglike,   myprior,
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
                    sampling_efficiency = 0.05, # very low eff. in
                                                #case of const efficiency mode,
                                                #README
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


if __name__=='__main__':
    import gl_params
    gp = gl_params.Params()

    run(gp)

# works with investigation = 'obs', pops = 2, metalpop = True
# profile with python3 -m cProfile grd_split.py
### TODO output:
