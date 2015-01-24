#!/usr/bin/env ipython3

##
# @file
# calculations for velocity dispersion
# disc version

# (c) GPL v3 2014 Pascal Steger, pascal@steger.aero

import pdb
import numpy as np
from scipy.interpolate import splrep, splint
import gi_helper as gh


def nr(z0, dlr, pop, gp):
    # extend asymptotes to 0, and high radius
    znu = np.hstack([z0[0]/gp.rinfty, z0, gp.rinfty*z0[-1]])
    logznu = np.log(znu/gp.Xscale[pop])
    dlr = -1.*dlr

    # use linear spline interpolation in r
    spline_n = splrep(logznu, dlr, k=1)

    # evaluate spline at any points in between
    return spline_n #, splev(z0, spline_n)
## \fn nr(z0, dlr, pop, gp)
# calculate n(r) at any given radius, as linear interpolation with two asymptotes
# @param z0 bin centers, [pc]
# @param dlr : asymptote at 0, n(r) for all bins, asymptote at infinity
# @param pop int for population
# @param gp global parameters


def rho(z0, rhopar, pop, gp):
    vec = rhopar
    rho_at_rhalf = vec[0]
    vec = vec[1:]
    # get spline representation on gp.xepol, where rhopar are defined on
    spline_n = nr(gp.xepol, vec, pop, gp)

    # and apply it to these radii, which may be anything in between
    zs =  np.log(z0/gp.Xscale[pop]) # have to integrate in d log(r)
    logrright = []; logrleft = []
    if np.rank(zs) == 0:
        if zs>0:
            logrright.append(zs)
        else:
            logrleft.append(zs)
    else:
        logrright = zs[(zs>=0.)]
        logrleft  = zs[(zs<0.)]
        logrleft  = logrleft[::-1] # inverse order

    # integrate to left and right of halflight radius
    logrhoright = []
    for i in np.arange(0, len(logrright)):
        logrhoright.append(np.log(rho_at_rhalf) + \
                           splint(0., logrright[i], spline_n))
                           # integration along dlog(r) instead of dr

    logrholeft = []
    for i in np.arange(0, len(logrleft)):
        logrholeft.append(np.log(rho_at_rhalf) + \
                          splint(0., logrleft[i], spline_n))

    tmp = np.exp(np.hstack([logrholeft[::-1], logrhoright])) # still defined on log(r)
    gh.checkpositive(tmp, 'rho()')
    return tmp
## \fn rho(z0, vec, pop, gp)
# calculate density, from interpolated n(r) = -log(rho(r))
# using interpolation to left and right of r=r_{*, 1/2}
# @param vec: rho(rstarhalf), asymptote_0, nr(gp.xipol), asymptote_infty
# @param z0: radii to calculate density for, in physical units (pc)
# @param pop int for which population
# @param gp global parameters


def beta(zipol, param, gp):
    betatmp = 0.
    for i in range(gp.nbeta):
        betatmp += param[i]*(zipol/max(gp.xipol))**i
    return betatmp, betatmp/(2.+betatmp)
## \fn beta(zipol, param, gp)
# return sum of polynomials for tilt as fct of radius
# @param zipol [pc]
# @param param n_beta parameters
# @param gp global parameters


def tilt(zipol, param, gp):
    tilttmp = 0.
    for i in range(gp.nbeta):
        tilttmp += param[i]*(zipol/max(zipol))**i
    return tilttmp
## \fn tilt(zipol, param, gp)
# return sum of polynomials for tilt as fct of radius
# @param zipol [pc]
# @param param n_beta parameters
# @param gp global parameters


def kappa(xipol, Kz):
    z0 = xipol
    kzpar = -np.hstack([Kz[0]/z0[0], (Kz[1:]-Kz[:-1])/(z0[1:]-z0[:-1])])
    return kzpar
## \fn kappa(xipol, Kz)
# calculate vertical force from Kz
# @param xipol height above midplane [pc]
# @param Kz acceleration in z direction


def nu_decrease(zpars, pars, gp):
    parsu = abs(pars)                        # Mirror prior
    if gp.monotonic:
        rnuz_z = np.zeros(len(zpars))
        rnuz_z[0] = parsu[0]
        for i in range(1,len(rnuz_z)):
            rnuz_z[i] = rnuz_z[i-1] + parsu[i]
        fun = rnuz_z[::-1]
    else:
        # Alternative here --> don't assume monotonic!
        fun = parsu

    # Normalise:
    if not gp.quadratic:
        # Exact linear interpolation integral:
        norm_nu = 0.
        for i in range(len(zpars)-1):
            zl = zpars[i]; zr = zpars[i+1]; zz = zr-zl
            b = fun[i]; a = fun[i+1]; norm_nu = norm_nu + (a-b)/2./zz*(zr**2.-zl**2.) + b*zr - a*zl
    else:
        # Exact quadratic interpolation:
        norm_nu = 0.
        for i in range(1,len(zpars)-1):
            z0 = zpars[i-1]; z1 = zpars[i]; z2 = zpars[i+1]
            f0 = fun[i-1];   f1 = fun[i];   f2 = fun[i+1]
            h = z2-z1
            a = f0; b = (f1-a)/h; c = (f2-2.*b*h-a)/2./h**2.
            z1d = z1-z0; z2d = z1**2.-z0**2.; z3d = z1**3.-z0**3.

            if i == len(zpars)-2:
                # Last bin integrate from z0 --> z2:
                z1d = z2-z0; z2d = z2**2.-z0**2.; z3d = z2**3.-z0**3.

            intquad = (a - b*z0 + c*z0*z1)*z1d + (b/2. - c/2.*(z0+z1))*z2d + c/3.*z3d
            norm_nu = norm_nu + intquad

            #sel  = (z > z0 and z < z2)
            #zcut = z[sel]
            #tcut = testy[sel]

    fun /= norm_nu

    # Check for negative density:
    sel = (fun>0)
    small = min(fun[sel])
    for jj in range(len(fun)):
        if (fun[jj] < 0):
            fun[jj] = small
    return fun
## \fn nu_decrease(zpars, pars, gp)
# Fully non-parametric monotonic *decreasing* function [c.f. Kz function].
# Note, here we explicitly avoid rnuz_z(0) = 0 since this would correspond
# to a constraint rnu_z(zmax) = 0 which is only true if zmax = infty.
# @param zpars [pc]
# @param pars nu parametrization in bins
# @param gp global parameters


def kz(zpars, kzpar, gp):
    # Solve interpolated integral for Kz:
    if not gp.quadratic:
        # Linear interpolation here:
        kz_z = np.zeros(len(zpars))
        for i in range(1, len(zpars)):
            zl = zpars[i-1]
            zr = zpars[i]
            zz = zr-zl
            b = kzpar[i-1]
            a = kzpar[i]
            kz_z[i] = kz_z[i-1]+(a-b)/2./zz*(zr**2.-zl**2.)+b*zr-a*zl
    else:
        # Quadratic interpolation here:
        kz_z = np.zeros(gp.nrho)
        for i in range(1,gp.nrho-1):
            z0 = zpars[i-1]
            z1 = zpars[i]
            z2 = zpars[i+1]
            f0 = kzpar[i-1]
            f1 = kzpar[i]
            f2 = kzpar[i+1]
            h = z2-z1
            a = f0
            b = (f1-a)/h
            c = (f2-2.*b*h-a)/2./h**2.
            z1d = z1-z0
            z2d = z1**2.-z0**2.
            z3d = z1**3.-z0**3.
            intbit = (a-b*z0+c*z0*z1)*z1d+(b/2.-c/2.*(z0+z1))*z2d+c/3.*z3d
            kz_z[i] = kz_z[i-1] + intbit
            if i == len(zpars)-2:
                # Deal with very last bin:
                z1d = z2-z1; z2d = z2**2.-z1**2.; z3d = z2**3.-z1**3.

                intbit = (a - b*z0 + c*z0*z1)*z1d + (b/2. - c/2.*(z0+z1))*z2d + c/3.*z3d
                kz_z[i+1] = kz_z[i] + intbit

    # Error checking. Sometimes when kz_z(0) << kz_z(1),
    # the interpolant can go negative. This just means that
    # we have resolved what should be kz_z(0) = 0. A simple
    # fix should suffice:
    for jj in range(0,len(kz_z)):
        if (kz_z[jj] < 0):
            kz_z[jj] = 0.

    return kz_z
## \fn kz(zpars, kzpar, gp)
# General function to calculate the Kz force law:
# Non-parametric Kz function here.
# K_z(z)= int_0^z kappa(z')dz'
#       = Sigma_z(z)*2*pi*G
#       = int_0^z 2 pi G rho(z') dz'
# ensure monotinicity in an efficient manner. Note here we
# use kz_z(0) = kzpar(0) * dz so that it can be zero, or
# otherwise just small in the case of large bins. This latter
# avoids the interpolated kz_out[0] going negative for coarse
# binning.
# @param zpars z where kzpar is defined on [pc]
# @param kzpar rho in bins [Msun/pc^3]
# @param gp global parameters


def sigz(zp, rhopar, lbaryonpar, MtoL, nupar, norm, tpar, pop, gp):
    # calculate density and Kz force:
    nutmp = rho(zp, nupar, pop, gp)
    nu_z = nutmp/np.max(nutmp)  # normalized to [1]
    gh.checkpositive(nu_z)
    rhodmtmp = rho(zp, rhopar, 0, gp) # rho_DM in linearly spaced bins
    gh.checkpositive(rhodmtmp)
    rhostartmp = rho(zp, lbaryonpar, 0, gp)
    gh.checkpositive(rhostartmp)
    rhotmp = rhodmtmp+MtoL*rhostartmp # add baryons
    kz_z = kz(zp, rhotmp, gp) # [(km/s)^2/pc]

    # add tilt correction [if required]:
    #Rsun = tparsR[0]
    #hr   = tparsR[1]
    #hsig = tparsR[2]
    #tc = sig_rz(zp, zp, tpars)
    #tc = tc * (1.0/Rsun - 1.0/hr - 1.0/hsig)
    # flag when the tilt becomes significant:
    #if abs(np.sum(tc))/abs(np.sum(kz_z)) > 0.1:
    #    print('Tilt > 10%!', abs(np.sum(tc)), abs(np.sum(kz_z)))
    #kz_z = kz_z-tc

    # do exact integral
    if not gp.quadratic:
        # linear interpolation here:
        sigint = np.zeros(len(zp))
        for i in range(1,len(zp)):
            zl = zp[i-1];  zr = zp[i];  zz = zr-zl
            b = nu_z[i-1]
            a = nu_z[i]
            q = kz_z[i-1]
            p = kz_z[i]

            intbit = (a-b)*(p-q)/(3.0*zz**2.) * (zr**3.-zl**3.) + \
                ((a-b)/(2.0*zz) * (q-(p-q)*zl/zz) + (p-q)/(2.0*zz)  * \
                 (b-(a-b)*zl/zz)) * (zr**2.-zl**2.)+\
                (b-(a-b)*zl/zz) * (q-(p-q)*zl/zz) * zz

            if i==0:
                sigint[0] = intbit
            else:
                sigint[i] = sigint[i-1] + intbit
    else:
        # quadratic interpolation
        sigint = np.zeros(len(zp))
        for i in range(1,len(zp)-1):
            z0 = zp[i-1];   z1 = zp[i];   z2 = zp[i+1]
            f0 = nu_z[i-1]; f1 = nu_z[i]; f2 = nu_z[i+1]
            fd0 = kz_z[i-1];fd1 = kz_z[i];fd2 = kz_z[i+1]
            h = z2-z1; a = f0; b = (f1-a)/h; c = (f2-2.*b*h-a)/2./h**2.
            ad = fd0; bd = (fd1-ad)/h; cd = (fd2-2.*bd*h-ad)/2./h**2.
            AA = a*bd+ad*b
            BB = c*ad+cd*a
            CC = b*cd+bd*c
            z1d = z1-z0
            z2d = z1**2.-z0**2.
            z3d = z1**3.-z0**3.
            z4d = z1**4.-z0**4.
            z5d = z1**5.-z0**5.
            intbit = z1d * (a*ad - z0*(AA-BB*z1) + z0**2.*(b*bd-CC*z1) + \
                        c*cd*z0**2.*z1**2.) + \
                        z2d * (0.5*(AA-BB*z1)-z0*(b*bd-CC*z1)-z0/2.*BB+z0**2./2.*CC - \
                        (z0*z1**2.+z0**2.*z1)*c*cd)+z3d*(1.0/3.0*(b*bd-CC*z1)+1.0/3.0*BB-2.0/3.0*z0*CC + \
                        1.0/3.0*c*cd*(z1**2.+z0**2.+4.0*z0*z1))+z4d*(1.0/4.0*CC-c*cd/2.0*(z1 + z0)) + \
                        z5d*c*cd/5.

            sigint[i] = sigint[i-1] + intbit
            if i == len(zp)-2:
                # Deal with very last bin:
                z1d = z2-z1
                z2d = z2**2.-z1**2.
                z3d = z2**3.-z1**3.
                z4d = z2**4.-z1**4.
                z5d = z2**5.-z1**5.
                intbit = z1d * (a*ad - z0*(AA-BB*z1) + z0**2.*(b*bd-CC*z1) + \
                         c*cd*z0**2.*z1**2.) + \
                  z2d * (0.5*(AA-BB*z1)-z0*(b*bd-CC*z1) - z0/2.*BB + z0**2./2.*CC - \
                         (z0*z1**2. + z0**2.*z1)*c*cd) + \
                  z3d * (1.0/3.0*(b*bd-CC*z1)+1.0/3.0*BB - 2.0/3.0*z0*CC + \
                         1.0/3.0*c*cd*(z1**2. + z0**2. + 4.0*z0*z1)) + \
                  z4d * (1.0/4.0*CC - c*cd/2.0 * (z1 + z0)) + \
                  z5d * c*cd / 5.

                sigint[i+1] = sigint[i] + intbit
    sig_z_t2 = 1.0/nu_z * (sigint + norm)
    return  np.sqrt(sig_z_t2[3:-3])
## \fn sigz(zp, rhopar, lbaryonpar, MtoL, nupar, norm, tpars, pop, gp)
# calculate z velocity dispersion
# @param zp vertical coordinate [pc]
# @param rhopar rho to be integrated to give force K_z
#SS No, rhotmp (=rhodmtmp+MtoL*rhostartmp) integrated to give K_z
# @param lbaryonpar overall baryonic light profile [Lsun]
# @param MtoL mass-to-light ratio [Msun/Lsun]
# @param nupar parameters for tracer density falloff
# @param norm value C, corresponds to integration from 0 to rmin
# @param tpars tilt parameters
# @param pop int for population (diff. halflight-radius)
# @param gp global parameters

def sig_rz(z, zpars, tpars):
    # Mirror prior
    tparsu = abs(tpars)

    # dz = zpars[2]-zpars[1]
    # sig_Rz = np.zeros(gp.nipol)
    # sig_Rz[0] = tparsu[0] * dz / 2.
    # for i in range(1,gp.nipol):
    #   sig_Rz[i] = sig_Rz[i-1] + tparsu[i] * dz
    # f = gp.ipol(zpars,sig_Rz,z)

    # Alternative here --> don't assume monotonic!
    f = gh.ipol(zpars, tparsu, z)
    return f
## \fn sig_rz(z, zpars, tpars)
# General function to describe the tilt profile
# @param z [pc]
# @param zpars [pc] z, on which sig is defined
# @param tpars tilt parameters: Rsun, hr, hsig

def nu_SUM_Sig(binmin, binmax, nudat):
    if len(binmin) != len(nudat):
        raise Exception('incompatible sizes of binmin and nudat')
    counts = nudat*(binmax-binmin)
    Sig = np.cumsum(counts)
    return Sig
## \fn nu_SUM_Sig(binmin, binmax, nudat)
# sum up nu to get surface density Sig
# @param binmin minimum of bin ranges in [pc]
# @param binmax maximum of bin ranges in [pc]
# @param nudat 3D density in bins, [Msun/pc^3]
