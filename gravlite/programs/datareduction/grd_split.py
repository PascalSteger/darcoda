#!/usr/bin/env python3

## @file
# split populations in observed dwarf galaxies
# based on metallicity, Mg index, v_LOS, position
# convention: pop = 0 is for MW population
#                   1, 2 is for first, second component

# (c) 2013 Pascal S.P. Steger

import sys, ipdb
import numpy as np

from scipy.integrate import simps
import pymultinest

import gr_params as gpr
import gl_helper as gh
from gl_centering import com_shrinkcircle_v_2D
import BiWeight as BW
import gl_units as gu

gh.DEBUGLEVEL = 1
DEBUG = True


def p_plummer(R, rs):
    logev_plummer = np.log(2.*R/rs**2/(1.+R**2/rs**2)**2)
    gh.sanitize_vector(logev_plummer, -1, -1e30, 1e6, DEBUG)
    return logev_plummer
## \fn p_plummer(R, rh)
# eq. 8 Walker 2011, likelihood that a tracer star is member of Plummer sphere
# @param R projected radius from center, [pc]
# @param rs scale radius, [pc]


def p_gauss(X, Xmean, sigmaX, errorX):
    prefactor = 1./np.sqrt(2.*np.pi*(sigmaX**2+errorX**2))
    exponent = -0.5*(X-Xmean)**2/(sigmaX**2+errorX**2)
    logev_gauss = np.log(prefactor)+exponent
    gh.LOG(3,'prefactor = ',prefactor)
    gh.LOG(3,'exponent = ',exponent)
    gh.sanitize_scalar(logev_gauss, -1e30, 1e6, DEBUG)
    return logev_gauss
## \fn p_gauss(X, Xmean, sigmaX, errorX)
# eq. 9, 11 Walker 2011, log likelihood based on generic Gauss function
# @param X variable, property of stellar tracer
# @param Xmean mean of all stars in that population
# @param sigmaX spread of Gaussian
# @param errorX observation error


def lp_MW(X, Xi, Xierror, PMi):
    nom = 0.
    denom = 0.
    for k in range(len(Xi)):
        nom += (1.-PMi[k])/np.sqrt(2.*np.pi*Xierror[k]**2)*np.exp(-0.5*(Xi[k]-X)**2/Xierror[k]**2)
        denom += (1.-PMi[k])
    prob_MW =  nom/denom
    gh.sanitize_scalar(np.log(prob_MW), -1e30, 1e6, DEBUG)
    return np.log(prob_MW)
## \fn lp_MW(X, Xi, Xierror, PMi, error)
# eq. 12 generic Walker 2011, likelihood that star is MW foreground
# @param X variable, property of one star
# @param Xi properties of stellar tracers
# @param Xierror errors on these
# @param PMi probability of membership
# @param error observation error


def lpR(Rk, pop):
    gh.LOG(3,'pR')
    log_prob_R = p_plummer(Rk, Rhalf_i[pop])
    return log_prob_R
## \fn lpR(Rk, pop)
# eq. 8 Walker 2011, probability distribution of radii
# @param Rk radius of stellar tracer k, [pc]
# @param pop int for population (0: MW, 1: 1, ...)


def lpV(Vk, k, pop):
    gh.LOG(3,'pV')
    # mean LOS velocity
    Vmean = calc_Vmean(alpha_s[k], delta_s[k])
    log_prob_V = p_gauss(Vk, Vmean, sigmaV[pop], Ve0[k])
    return log_prob_V
## \fn lpV(Vk, pop)
# eq. 9 Walker 2011, probability distribution of LOS velocities
# @param Vk velocity of stellar tracer k, [km/s]
# @param k ID of tracer star under investigation
# @param pop int for population (0: MW, 1, 2..)


def lpW(Wk, k, pop):
    gh.LOG(3,'pW')
    log_prob_W = p_gauss(Wk, Wmean[pop], sigmaW[pop], We0[k])
    return log_prob_W
## \fn lpW(Wk, k, pop)
# eq. 11 Walker 2011, probability distribution of reduced magnesium index
# @param Wk magnesium index, [Ang]
# @param k ID for stellar tracer
# @param pop int for population (0: MW, 1, 2..)


def pjoint(R, k2, V, Verror, W, Werror, PM, k, pop):
    gh.LOG(2,'pjoint for tracer star ', k)
    gh.LOG(2,' of ', len(PM))
    gh.LOG(2,' and pop ', pop)
    if pop == 0:
        lpr = lp_MW(R[k], R, k2, PM) # k2 is measured half-light radius of overall distro
        lpv = lp_MW(V[k], V, Verror, PM)
        lpw = lp_MW(W[k], W, Werror, PM)
        log_p_joint = lpr+lpv+lpw
    else:
        lpr = lpR(R[k], pop)
        lpv = lpV(V[k], k, pop)
        lpw = lpW(W[k], k, pop)
        log_p_joint = lpr+lpv+lpw
    gh.sanitize_scalar(log_p_joint, -1e30, 0, DEBUG)
    return np.exp(log_p_joint)
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



def myprior(cube, ndim, nparams):
    # convert to physical space
    off = 0
    cube[off] = cube[off] # fmem
    off += 1
    cube[off] = cube[off] # fsub
    off += 1
    cube[off] = 10**(4.5*cube[off]) # r_half [pc]
    off += 1
    cube[off] = cube[off]*2000.-1000. # proper motion in x [mas/century], mu_alpha
    off += 1
    cube[off] = cube[off]*2000.-1000. # proper motion in y [mas/century], mu_delta
    off += 1
    for pop in range(0, gp.pops+1): # no. of pops goes in here, first MW, then 1,2,..
        cube[off] = cube[off] # Rhalf_i / Rhalf
        off += 1
        cube[off] = cube[off]*6.-3. # Wmean
        off += 1
        cube[off] = np.sqrt(10.**(cube[off]*6.-5.))  # sigmaW [Ang]
        off += 1
        cube[off] = np.sqrt(10.**(cube[off]*10.-5.)) # sigmaV [km/s]
        off += 1
    if off != ndim:
        gh.LOG(1, 'wrong number of parameters in myprior.cube')
        ipdb.set_trace()
    return
## \fn myprior(cube, ndim, nparams) priors
# @param cube [0,1]^ndim cube, array of dimension ndim
# @param ndim number of dimensions, 2*npop*nipol + nipol
# @param nparams = ndim + additional parameters
# stored with actual parameters



def w(Rk):
    gh.sanitize_vector(Rk, Nsample, 0, 1e30, True)
    w_ipol = np.zeros(Nsample)
    for k in range(Nsample):
        w_ipol[k] = wpt[np.where(abs(Rk[k]-Rpt) == min(abs(Rk[k]-Rpt)))]
    return w_ipol
## \fn w(Rk)
# return selection function as function of radius
# take vector as input
# @param Rk radius [pc]

glob_w = w_pt(R0)


def int_wp(Rhalf):
    integral = simps(w(R0)*p_plummer(R0, Rhalf), R0)
    gh.sanitize_scalar(integral, 0, 1e30, True)
    return integral
## \fn int_wp(pop, Rhalf)
# calculate denominator integral for population pop=1,2 in eq. 14 Walker 2011
# @param Rhalf half-light radius for Plummer spheres of pop 1 and pop 2
# @param gp global parameters, for radii


gh.LOG(2, 'calculate integrals in denominator')
iiinteg = []
integ.append(1)             # for pop 0, MW, want to divide by 1
integ.append(int_wp(Rhalf_i[1]))     # for pop 1
integ.append(int_wp(Rhalf_i[2]))   # for pop 2

def calc_Vmean(als, des):
    gh.sanitize_vector(als, Nsample, -1e6, 1e6)
    gh.sanitize_vector(des, Nsample, -1e6, 1e6)
    mu_alpha__ascentury_1 = mu_alpha / gu.arcsec__mas
    mu_alpha__ass_1 = mu_alpha__ascentury_1 / gu.century__s
    mu_alpha__rad_per_s = mu_alpha__ass_1 / gu.rad__arcsec
    mu_delta__ascentury_1 = mu_delta / gu.arcsec__mas
    mu_delta__ass_1 = mu_delta__ascentury_1 / gu.century__s
    mu_delta__rad_per_s = mu_delta__ass_1 / gu.rad__arcsec
    DL__km = DL * gu.pc__m / gu.km__m # [km]
    term1 = np.cos(des)*np.sin(als)*(Vd*np.cos(dd)*np.sin(ad)+\
                           DL__km*mu_alpha__rad_per_s*np.cos(dd)*np.cos(ad)-\
                           DL__km*mu_delta__rad_per_s*np.sin(dd)*np.sin(ad))
    term2 = np.cos(des)*np.cos(als)*(Vd*np.cos(dd)*np.cos(ad)\
                           -DL__km*mu_delta__rad_per_s*np.sin(dd)*np.cos(ad)\
                           -DL__km*mu_alpha__rad_per_s*np.cos(dd)*np.sin(ad))
    term3 = np.sin(des)*(Vd*np.sin(dd)+DL__km*mu_delta__rad_per_s*np.cos(dd))
    Vm = term1 + term2 + term3
    gh.LOG(3,' Vmean =',Vm)
    gh.sanitize_scalar(Vm, -1000, 1000, DEBUG)
    return Vm
## \fn calc_Vmean(als, des)
# eq. 10 Walker 2011, dwarf spheroidal systemic HRF, [km/s]
# @param als right ascension of star in [arcsec]
# @param des declination of star in [arcsec]

glob_Vm = calc_Vmean(alpha_s, delta_s)

glob_sum_1_PM = np.sum(1-PM)
glob_M_r = np.zeros((Nsample, Nsample)); glob_phat_r = np.zeros(Nsample)
glob_M_v = np.zeros((Nsample, Nsample)); glob_phat_v = np.zeros(Nsample)
glob_M_w = np.zeros((Nsample, Nsample)); glob_phat_w = np.zeros(Nsample)

for i in range(Nsample):
    prefac = (1-PM[i])/np.sqrt(2*np.pi)
    glob_M_r[i,:] = prefac/k2*np.exp(-(R0[i]-R0)**2/(2*k2*k2))
    glob_M_v[i,:] = prefac/epsv*np.exp(-(V0[i]-V0)**2/(2*Ve0*Ve0))
    glob_M_w[i,:] = prefac/epsw*np.exp(-(W0[i]-W0)**2/(2*We0*We0))
glob_phat_r = np.sum(glob_M_r, 0)/glob_sum_1_PM
glob_phat_v = np.sum(glob_M_v, 0)/glob_sum_1_PM
glob_phat_w = np.sum(glob_M_w, 0)/glob_sum_1_PM

def myloglike(cube, ndim, nparams):
    fmem = cube[0]
    fsub = cube[1]
    f1 = fmem*fsub
    f2 = fmem*(1-fsub)
    ftot = [1-f1-f2, f1, f2]
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
    if off != ndim:
        gh.LOG(1, 'wrong number of parameters in myloglike.cube')
        ipdb.set_trace()


    p_R_1 = TODO

    gh.LOG(2,'starting logev evaluation')
    logev = 0.0

        single0 = ftot[pop]*glob_w
        single2 = pjoint(R0, Rhalf_i[pop]*np.ones(len(v0)), v0, Ve0,\
                             W0, We0, PM0, k, pop)


        if np.isnan(single):
            ipdb.set_trace()
        if single < -1e30:
            ipdb.set_trace()
        term += single


        logev += np.log(term)
        gh.LOG(1,' found log(likelihood) = ',logev)
    return logev
## \fn myloglike(cube, ndim, nparams) calculate probability function
# @param cube ndim cube of physical parameter space (nr)
# @param ndim number of dimensions, 2*npop*nipol + nipol
# @param nparams = ndim + additional parameters
# stored with actual parameters


def run(gp):
    global Nsample, wpt, Rpt, V0, Ve0, W0, We0, PM0
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
    gh.LOG(3,'RAh: signum = ',sig)
    RAh = RAh/sig
    # stellar position alpha_s, delta_s
    # 15degrees in 1 hour right ascension

    alpha_s = 15*(RAh*3600+RAm*60+RAs)*sig       # [arcsec]

    sig = abs(DEd[0])/DEd[0]                # +/-
    gh.LOG(3,'DEd: signum = ',sig)
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


    # determine distance to dwarf
### TODO reference
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

    # remove center displacement, already change x0
    com_x, com_y, com_vz = com_shrinkcircle_v_2D(x0, y0, VHel, PM) # [pc], [km/s]

    R0 = np.sqrt(x0**2+y0**2)
    global R0, Rfine
    Rfine = np.logspace(np.log10(min(R0)), np.log10(max(R0)), 100)
    V0 = 1.*np.copy(VHel) # [km/s] not necessary to remove center LOS velocity
    Ve0 = 1.*e_VHel # velocity error

    Nsample = len(PM)

    A = np.loadtxt(gp.files.dir+'w_2.0.dat')
    Rpt, wpt = A.T # [arcmin], [1]
    arcmin = 2.*np.pi* DL / (360 * 60) # [pc] at distance 138 kpc for Fornax
    Rpt *= arcmin # [pc]


    n_dims = 5+(gp.pops+1)*4

    gh.LOG(1,'starting MultiNest run:')
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
                    n_live_points = 3,
                    ### TODO gp.nlive,
                    evidence_tolerance = 0.0, # set to 0 to keep
                                              #algorithm working
                                              #indefinitely
                    sampling_efficiency = 0.05, # very low eff. in
                                                #case of const efficiency mode,
                                                #README
                    n_iter_before_update = 1, # output after this many iterations
                    null_log_evidence = 1., # separate modes if
                                            #logevidence > this param.
                    max_modes = 3,
                    ### TODO gp.nlive,   # preallocation of modes:
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
                    log_zero = -999999,      # points with log likelihood
                                          #< log_zero will be
                                          #neglected
                    max_iter = 1,
                    ### TODO set to 0 for never
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
### TODO output:
