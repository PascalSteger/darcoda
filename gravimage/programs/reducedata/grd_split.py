#!/usr/bin/env ipython3

## @file
# split populations in observed dwarf galaxies
# based on metallicity, Mg index, v_LOS, position
# convention: pop = 0 is for MW population
#                   1, 2 is for first, second component

# (c) 2013 Pascal S.P. Steger

import ipdb
import numpy as np

from scipy.integrate import simps
import pymultinest

import gi_helper as gh
from gi_centering import com_shrinkcircle_v_2D
#import BiWeight as BW
import gi_units as gu

gh.DEBUGLEVEL = 3
DEBUG = True

def lp_plummer(Rad, rs):
    logev_plummer = np.log(2.*Rad/rs**2/(1.+Rad**2/rs**2)**2)
    gh.sanitize_vector(logev_plummer, -1, -1e30, 1e6, DEBUG)
    return logev_plummer
## \fn lp_plummer(Rad, rs)
# eq. 8 Walker 2011, likelihood that a tracer star is member of Plummer sphere
# @param Rad projected radius from center, [pc]
# @param rs scale radius, [pc]


def lp_gauss(X, Xmean, sigmaX, errorX):
    prefactor = 1./np.sqrt(2.*np.pi*(sigmaX**2+errorX**2))
    exponent = -0.5*(X-Xmean)**2/(sigmaX**2+errorX**2)
    logev_gauss = np.log(prefactor)+exponent
    gh.LOG(3,'prefactor = ',prefactor)
    gh.LOG(3,'exponent = ',exponent)
    gh.sanitize_vector(logev_gauss, Nsample, -1e30, 1e6, DEBUG)
    return logev_gauss
## \fn lp_gauss(X, Xmean, sigmaX, errorX)
# eq. 9, 11 Walker 2011, log likelihood based on generic Gauss function
# @param X variable, property of stellar tracer
# @param Xmean mean of all stars in that population
# @param sigmaX spread of Gaussian
# @param errorX observation error


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
    gh.sanitize_vector(Rk, Nsample, 0, 1e30, DEBUG)
    w_ipol = np.zeros(Nsample)
    for k in range(Nsample):
        w_ipol[k] = wpt[np.where(abs(Rk[k]-Rpt) == min(abs(Rk[k]-Rpt)))]
    return w_ipol
## \fn w(Rk)
# return selection function as function of radius
# take vector as input
# @param Rk radius [pc]


def int_wp(Rhalf):
    integral = simps(gw*np.exp(lp_plummer(R0, Rhalf)), R0)
    gh.sanitize_scalar(integral, 0, 1e30, DEBUG)
    return integral
## \fn int_wp(pop, Rhalf)
# calculate denominator integral for population pop=1,2 in eq. 14 Walker 2011
# @param Rhalf half-light radius for Plummer spheres of pop 1 and pop 2
# @param gp global parameters, for radii


def calc_Vmean(als, des):
    gh.sanitize_vector(als, Nsample, -1e6, 1e6, DEBUG)
    gh.sanitize_vector(des, Nsample, -1e6, 1e6, DEBUG)
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
    gh.sanitize_vector(Vm, Nsample, -500, 500, DEBUG)
    return Vm
## \fn calc_Vmean(als, des)
# eq. 10 Walker 2011, dwarf spheroidal systemic HRF, [km/s]
# @param als right ascension of star in [arcsec]
# @param des declination of star in [arcsec]


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

    gh.LOG(2, 'calculate integrals in denominator')

    lpR1 = lp_plummer(R0, Rhalf_i[1])
    lpR2 = lp_plummer(R0, Rhalf_i[2])
    gVm = calc_Vmean(alpha_s, delta_s)
    lpV1 = lp_gauss(V0, gVm, sigmaV[1], Ve0)
    lpV2 = lp_gauss(V0, gVm, sigmaV[2], Ve0)
    lpW1 = lp_gauss(W0, Wmean[1], sigmaW[1], We0)
    lpW2 = lp_gauss(W0, Wmean[2], sigmaW[2], We0)

    gh.LOG(2,'starting logev evaluation')
    term_MW = ftot[0]*gphat_r*gphat_v*gphat_w
    gintw1 = int_wp(Rhalf_i[1])     # for pop 1
    gintw2 = int_wp(Rhalf_i[2])   # for pop 2
    term_pop1 = ftot[1]*gw*np.exp(lpR1+lpV1+lpW1)/gintw1
    term_pop2 = ftot[2]*gw*np.exp(lpR2+lpV2+lpW2)/gintw2

    sumterms = term_MW+term_pop1+term_pop2
    print('f_MW,pop1,pop2=',np.median(term_MW/sumterms),\
          np.median(term_pop1/sumterms), np.median(term_pop2/sumterms))
    ipdb.set_trace()
    logterm14sum = np.log(sumterms)
    logev = np.sum(logterm14sum)
    gh.sanitize_scalar(logev, -1e30, 1e6, DEBUG)

    # TODO write out probability of membership in MW, pop1, pop2


    gh.LOG(1, 'logL:',logev)
    return logev
## \fn myloglike(cube, ndim, nparams) calculate probability function
# @param cube ndim cube of physical parameter space (nr)
# @param ndim number of dimensions, 2*npop*nipol + nipol
# @param nparams = ndim + additional parameters
# stored with actual parameters


def run(gp):
    import gr_params
    gpr = gr_params.grParams(gp)

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


    global wpt, Rpt, V0, Ve0, W0, We0, PM0
    gpr.fil = gpr.dir+"/table_merged.bin"
    delim = [0,22,3,3,6,4,3,5,6,6,7,5,6,5,6,5,6]
    #ID = np.genfromtxt(gpr.fil, skiprows=29, unpack=True,\
    #                   usecols=(0,1),delimiter=delim)
    RAh,RAm,RAs,DEd,DEm,DEs,Vmag,VI,\
      VHel,e_VHel,SigFe,e_SigFe,\
      SigMg,e_SigMg,PM = np.genfromtxt(gpr.fil, skiprows=29, unpack=True, \
                                       usecols=tuple(range(2,17)), \
                                       delimiter=delim, filling_values=-1)
    # attention, we do not have Mg measurements for 501 stars in Fornax,
    #  visible by missing SigMg values, set to -1
    #   we exclude them from all further analysis
    sel = (SigMg>-1)
    RAh=RAh[sel]; RAm=RAm[sel]; RAs=RAs[sel]; DEd=DEd[sel]; DEm=DEm[sel]; DEs=DEs[sel]
    Vmag=Vmag[sel]; VI=VI[sel]; VHel=VHel[sel]; e_VHel=e_VHel[sel]; PM=PM[sel]
    SigMg=SigMg[sel]; e_SigMg=e_SigMg[sel]; SigFe=SigFe[sel]; e_SigFe=e_SigFe[sel]

    global Nsample
    Nsample = len(PM)

    # use all stellar tracer particles from now on, independent on their probability of membership
    V0 = 1.*np.copy(VHel)  # [km/s] not necessary to remove center LOS velocity
    Ve0 = 1.*e_VHel        # velocity error
    W0 = 1.*np.copy(SigMg)
    We0 = 1.*np.copy(e_SigMg)

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
    W0 = W0[order]; We0 = We0[order]
    # Rfine = np.logspace(np.log10(min(R0)), np.log10(max(R0)), 100)

    A = np.loadtxt(gp.files.dir+'w_2.0.dat')
    Rpt, wpt = A.T # [arcmin], [1]
    arcmin__pc = 2.*np.pi* DL / (360 * 60) # [pc] at distance 138 kpc for Fornax
    Rpt *= arcmin__pc # [pc]
    global gw
    gw = w(R0)

    sum_1_PM = np.sum(1-PM)
    gM_r = np.zeros((Nsample, Nsample))
    gM_v = np.zeros((Nsample, Nsample))
    gM_w = np.zeros((Nsample, Nsample))

    # TODO check this is right
    for i in range(Nsample):
        prefac = (1-PM[i])/np.sqrt(2*np.pi)
        gh.sanitize_scalar(prefac, 0, 1/np.sqrt(2*np.pi), DEBUG)
        gM_r[i,:] = prefac/np.abs(k2)*np.exp(-((R0[i]-R0)**2)/(2*k2*k2))
        gM_v[i,:] = prefac/np.abs(Ve0)*np.exp(-((V0[i]-V0)**2)/(2*Ve0*Ve0))
        gM_w[i,:] = prefac/np.abs(We0)*np.exp(-((W0[i]-W0)**2)/(2*We0*We0))
    global gphat_r, gphat_v, gphat_w
    gphat_r = np.sum(gM_r, 0)/sum_1_PM
    gh.sanitize_vector(gphat_r, Nsample, 0, 1e30, DEBUG)
    gphat_v = np.sum(gM_v, 0)/sum_1_PM
    gh.sanitize_vector(gphat_v, Nsample, 0, 1e30, DEBUG)
    gphat_w = np.sum(gM_w, 0)/sum_1_PM
    gh.sanitize_vector(gphat_w, Nsample, 0, 1e30, DEBUG)

    gh.LOG(1,'starting MultiNest run:')
    n_dims = 5+(gp.pops+1)*4
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
                    n_iter_before_update = 1, # output after this many iterations
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
    import gi_params
    gp = gi_params.Params()

    run(gp)

# works with investigation = 'obs', pops = 2, metalpop = True
# profile with python3 -m cProfile grd_split.py
