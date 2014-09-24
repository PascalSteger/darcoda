#!/usr/bin/env python3

## @file
# split populations in observed dwarf galaxies
# based on metallicity, Mg index, v_LOS, position

# (c) 2013 Pascal S.P. Steger

import sys, pdb, ipdb
import numpy as np
from numpy import sin, cos

import pymultinest

import gr_params as gpr
import gl_helper as gh
from gl_centering import com_shrinkcircle_v_2D
import BiWeight as BW

def p_plummer(R, rs):
    return np.log(2.*R/rs**2/(1.+R**2/rs**2)**2)
## \fn p_plummer(R, rh)
# eq. 8 Walker 2011, likelihood that a tracer star is member of Plummer sphere
# @param R projected radius from center, [pc]
# @param rs scale radius, [pc]


def p_gauss(X, Xmean, sigmaX, errorX):
    prefactor = 1./np.sqrt(2.*np.pi*(sigmaX**2+errorX**2))
    exponent = -0.5*(X-Xmean)**2/(sigmaX**2+errorX**2)
    print('prefactor, exponent =',prefactor,exponent)
    return np.log(prefactor*np.exp(exponent))
## \fn p_gauss(X, Xmean, sigmaX, errorX)
# eq. 9, 11 Walker 2011, likelihood based on generic Gauss function
# @param X variable, property of stellar tracer
# @param Xmean mean of all stars in that population
# @param sigmaX spread of Gaussian
# @param errorX observation error


def p_MW(X, Xi, Xierror, PMi):
    nom = 0.
    denom = 0.
    for k in range(len(Xi)):
        nom += (1.-PMi[k])/np.sqrt(2.*np.pi*Xierror[k]**2)*np.exp(-0.5*(Xi[k]-X)**2/Xierror[k]**2)
        denom += (1.-PMi[k])
    return np.log(nom/denom)
## \fn p_MW(X, Xi, Xierror, PMi, error)
# eq. 12 generic Walker 2011, likelihood that star is MW foreground
# @param X variable, property of one star
# @param Xi properties of stellar tracers
# @param Xierror errors on these
# @param PMi probability of membership
# @param error observation error


def pR(Rk, pop):
    print('pR')
    return p_plummer(Rk, Rhalf_i[pop])
## \fn pR(Rk, pop)
# eq. 8 Walker 2011, probability distribution of radii
# @param Rk radius of stellar tracer k, [pc]
# @param pop int for population (0: MW, 1: 1, ...)


def pV(Vk, k, pop):
    print('pV')
    # mean LOS velocity
    Vmean = calc_Vmean(alpha_s[k], delta_s[k])

    return p_gauss(Vk, Vmean, sigmaV[pop], ve0[k])
## \fn pV(Vk, pop)
# eq. 9 Walker 2011, probability distribution of LOS velocities
# @param Vk velocity of stellar tracer k, [km/s]
# @param k ID of tracer star under investigation
# @param pop int for population (0: MW, 1, 2..)

def pW(Wk, k, pop):
    print('pW')
    return p_gauss(Wk, Wmean[pop], sigmaW[pop], We0[k])
## \fn pW(Wk, k, pop)
# eq. 11 Walker 2011, probability distribution of reduced magnesium index
# @param Wk magnesium index, [Ang]
# @param k ID for stellar tracer
# @param pop int for population (0: MW, 1, 2..)


def pjoint(R, k2, V, Verror, W, Werror, PM, k, pop):
    print('pjoint for tracer star ',k,' of ',len(PM))
    if pop == 0:
        pr = p_MW(R[k], R, k2, PM) # k2 is measured half-light radius of overall distro
        pv = p_MW(V[k], V, Verror, PM)
        pw = p_MW(W[k], W, Werror, PM)
        return pr*pv*pw
    else:
        return pR(R[k], pop)*pV(V[k], k, pop)*pW(W[k], k, pop)
## \fn pjoint(R, k2, V, W, PM, k, pop)
# eq. 13 Walker 2011, joint probability distributions
# @param R projected radius, [pc]
# @param k2 smoothing scale = half-light radius of overall distro
# @param V LOS velocity [km/s]
# @param Verror error on V
# @param W reduced magnesium index [Ang]
# @param Werror error on W
# @param PM probability of membership
# @param k int ID of star under investigation
# @param pop int for population


def calc_Vmean(als, des):
    # TODO debug

    term1 = cos(des)*sin(als)*(Vd*cos(dd)*sin(ad)+\
                           DL*mu_alpha*cos(dd)*cos(ad)-\
                           DL*mu_delta*sin(dd)*sin(ad))
    term2 = cos(des)*cos(als)*(Vd*cos(dd)*cos(ad)\
                           -DL*mu_delta*sin(dd)*cos(ad)\
                           -DL*mu_alpha*cos(dd)*sin(ad))
    term3 = sin(des)*(Vd*sin(dd)+DL*mu_delta*cos(dd))
    Vm = term1 + term2 + term3
    pdb.set_trace()
    print(' Vmean =',Vm)
    return Vm
## \fn calc_Vmean(als, des)
# eq. 10 Walker 2011, dwarf spheroidal systemic HRF, [km/s]
# @param als right ascension of star in [arcsec]
# @param des declination of star in [arcsec]


def w(Rk):
    w_ipol = wpt[np.where(abs(Rk-Rpt) == min(abs(Rk-Rpt)))]
    return w_ipol
## \fn w(Rk)
# return selection function as function of radius
# @param Rk radius [pc]

def show(filepath):
    subprocess.call(('xdg-open', filepath))
    return
## \fn show(filepath) open the output (pdf) file for the user @param
# filepath filename with full path


def myprior(cube, ndim, nparams):
    # convert to physical space
    cube[0] = cube[0] # fmem
    cube[1] = cube[1] # fsub
    cube[2] = 10**(4.5*cube[2]) # r_half [pc]

    cube[3] = cube[3]*2000.-1000. # proper motion in x [km/s], mu_alpha
    cube[4] = cube[4]*2000.-1000. # proper motion in y [km/s], mu_delta
    off = 5 # offset from common parameters
    for pop in range(gp.pops+1): # no. of pops goes in here, first MW, then 1,2,..
        cube[off] = cube[off] # Rhalf_i / Rhalf
        off += 1
        cube[off] = cube[off]*6.-3. # Wmean
        off += 1
        cube[off] = np.sqrt(10.**(cube[off]*6.-5.))  # sigmaW [Ang]
        off += 1
        cube[off] = np.sqrt(10.**(cube[off]*10.-5.)) # sigmaV [km/s]
        off += 1
    return
## \fn myprior(cube, ndim, nparams) priors
# @param cube [0,1]^ndim cube, array of dimension ndim
# @param ndim number of dimensions, 2*npop*nipol + nipol
# @param nparams = ndim + additional parameters
# stored with actual parameters


def myloglike(cube, ndim, nparams):
    ev = 0.
    fmem = cube[0]
    fsub = cube[1]
    f1 = fmem*fsub
    f2 = fmem*(1-fsub)
    ftot = [1-f1-f2, f1, f2]
    # TODO: calculate log likelihood of assigning stars to each population
    r_half = cube[2]

    global mu_alpha, mu_delta
    mu_alpha = cube[3]
    mu_delta = cube[4]
    off = 5

    global Rhalf_i, sigmaW, sigmaV, Wmean
    Rhalf_i = []; Wmean = []; sigmaW = []; sigmaV = []
    for pop in range(gp.pops+1):
        Rhalf_i.append(cube[off]*r_half)
        off += 1
        Wmean.append(cube[off])
        off += 1
        sigmaW.append(cube[off])
        off += 1
        sigmaV.append(cube[off])
        off += 1


    print('starting logev evaluation')
    logev = 0.0
    for k in range(Nsample):
        term = 0.0
        for pop in range(gp.pops+1):
            term += ftot[pop] * w(R0[k])*\
                    pjoint(R0, Rhalf_i[pop]*np.ones(len(v0)), v0, ve0,\
                           W0, We0, PM0, k, pop)
            # TODO divide by integral
        logev += np.log(term)

        print(' found log(likelihood) = ',logev)
    return logev
## \fn myloglike(cube, ndim, nparams) calculate probability function
# @param cube ndim cube of physical parameter space (nr)
# @param ndim number of dimensions, 2*npop*nipol + nipol
# @param nparams = ndim + additional parameters
# stored with actual parameters


def run(gp):
    global Nsample, R0, wpt, Rpt, v0, ve0, W0, We0, PM0
    global Vmean, alpha_s, delta_s
    gpr.fil = gpr.dir+"/table_merged.bin"
    delim = [0,22,3,3,6,4,3,5,6,6,7,5,6,5,6,5,6]
    ID = np.genfromtxt(gpr.fil, skiprows=29, unpack=True,\
                       usecols=(0,1),delimiter=delim)


    RAh,RAm,RAs,DEd,DEm,DEs,Vmag,VI,\
      VHel,e_VHel,SigFe,e_SigFe,\
      SigMg,e_SigMg,PM = np.genfromtxt(gpr.fil, skiprows=29, unpack=True, \
                                       usecols=tuple(range(2,17)), delimiter=delim, filling_values=-1)
    # use all stellar tracer particles from now on, independent on their probability of membership
    W0 = SigMg
    We0 = e_SigMg
    sig = abs(RAh[0])/RAh[0]
    print('RAh: signum = ',sig)
    RAh = RAh/sig
    # stellar position alpha_s, delta_s
    # 15degrees in 1 hour right ascension

    alpha_s = 15*(RAh*3600+RAm*60+RAs)*sig       # [arcsec/15]

    sig = abs(DEd[0])/DEd[0]                # +/-
    print('DEd: signum = ',sig)
    DEd = DEd/sig
    delta_s = (DEd*3600+DEm*60+DEs)*sig          # [arcsec]



    # unit conversion into a set of [pc], [km/s]
    arcsec = 2.*np.pi/(360.*60.*60)      # [rad/arcsec]
    kpc = 1000 # [pc]

    alpha_s *= arcsec  # [rad]
    delta_s *= arcsec  # [rad]


    # instead of using other datasets for individual dSph,
    # determine Heliocentric-Rest-Frame Line-Of-Sight velocity of dwarf Vd
    # and position of dwarf (ad, dd)
    # from probability-of-membership weighted means directly
    global Vd, ad, dd
    Vd = np.sum(VHel * PM)/np.sum(PM)    # [km/s]
    ad = np.sum(alpha_s * PM)/np.sum(PM) # [arcsec]
    dd = np.sum(delta_s * PM)/np.sum(PM) # [arcsec]


    # determine distance to dwarf (TODO: reference)
    global DL
    DL = {0: lambda x: x * (138),#+/- 8 for Fornax
          1: lambda x: x * (101),#+/- 5 for Carina
          2: lambda x: x * (79),  #+/- 4 for Sculptor
          3: lambda x: x * (86) #+/- 4 for Sextans
      }[gp.case](kpc)

    xs = alpha_s*DL # [pc]
    ys = delta_s*DL # [pc]

    PM0 = 1.*np.copy(PM)

    x0 = 1.*np.copy(xs)
    y0 = 1.*np.copy(ys) # [pc]
    R0 = np.sqrt(x0*x0+y0*y0)
    v0 = 1.*np.copy(VHel) # [km/s]
    ve0 = 1.*e_VHel # velocity error

    Nsample = len(PM)


    A = np.loadtxt(gp.files.dir+'w_2.0.dat')
    Rpt, wpt = A.T # [arcmin], [1]
    arcmin = 2.*np.pi* DL / (360 * 60) # [pc] at distance 138 kpc for Fornax
    Rpt *= arcmin # [pc]


    n_dims = 17

    print('starting MultiNest run:')
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
                    multimodal = False,            # separate modes
                    const_efficiency_mode = True, # use const sampling efficiency
                    n_live_points = gp.nlive,
                    evidence_tolerance = 0.0, # set to 0 to keep
                                              #algorithm working
                                              #indefinitely
                    sampling_efficiency = 0.05, # very low eff. in
                                                #case of const efficiency mode,
                                                #README
                    n_iter_before_update = 1, # output after this many iterations
                    null_log_evidence = 1., # separate modes if
                                            #logevidence > this param.
                    max_modes = gp.nlive,   # preallocation of modes:
                                            #max. = number of live
                                            #points
                    mode_tolerance = -1.e60,   # mode tolerance in the
                                               #case where no special
                                               #value exists: highly
                                               #negative
                    outputfiles_basename = gp.files.outdir,
                    seed = -1,
                    verbose = True,
                    resume = False,
                    context = 0,
                    write_output = True,
                    log_zero = 0,      # points with log likelihood
                                          #< log_zero will be
                                          #neglected
                    max_iter = 0,         # set to 0 for never
                                          #reaching max_iter (no
                                          #stopping criterium based on
                                          #number of iterations)
                    init_MPI = True,     # use MPI
                    dump_callback = None)


if __name__=='__main__':
    import gl_params
    gp = gl_params.Params()

    run(gp)

# works with investigation = 'obs', pops = 2, metalpop = True
# profile with python3 -m cProfile grd_split.py
# TODO output:
