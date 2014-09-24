#!/usr/bin/env ipython3

##
# @file
# all integrals from gl_physics

# (c) 2013 Pascal S.P. Steger

import numpy as np
import pdb, scipy
from scipy.integrate import simps,trapz,quad
from scipy.interpolate import splrep, splev, splint
import gl_helper as gh
import gl_plot as gpl
import gl_analytic as ga
import gl_physics as phys
import gl_project as glp
import gl_plot as gpl
from pylab import *

def int_poly_inf(r0,poly):
    f = -1/poly[0]*np.exp(poly[1]+poly[0]*r0)
    gh.checknan(f, 'f in int_poly_inf')
    return f
## \fn int_poly_inf(r0,poly)
# integrate polynomial from r0[i] to infinity
# assume log(integrand) = A + B*r;    B = poly[0] < 0, A = poly[1]
# @param r0 radii
# @param poly slope


def correct_first_bin(xint, yint, k=3, s=0.01, log=True):
    if log:
        yint[1:] = np.log(yint[1:])
    tmp = splev(xint[0], splrep(xint[1:], yint[1:], k=k, s=s))
    if log:
        tmp = np.exp(tmp)
    return tmp
## \fn correct_first_bin(xint, yint, k=3, s=0.01, log=True)
# extrapolate back to first bin, with splines
# @param xint x values
# @param yint y values, of which the first entry has to be rechecked
# @param k order of the spline
# @param s smoothing of the spline. 0.01 is a *small* default value, hard to converge to
# @param log bool for working in log space


def ant_intbeta(r0, betapar, gp):
    # define function
    xint = 1.*r0
    yint = phys.beta(xint, betapar, gp)[0]/xint

    # analytic values
    # yint =  ga.beta_gaia(xint, gp)[1]/xint

    intbet = np.zeros(len(r0))
    for i in range(len(r0)):
        intbet[i] = gh.quadinf(xint, yint, r0[0]/1e5, r0[i])
    # TODO: assumption here was that integration goes down to min(r0)/1e5

    gh.checknan(intbet, 'intbet in ant_intbeta')
    return intbet                                                # [1]
## \fn ant_intbeta(r0, betapar, gp)
# integrate beta(s)/s over s
# (integrals in front of and after sigma_r^2 integral, factor 2 not in here)
# @param r0 free variable, array, in [lunit]
# @param betapar integrand, array, in [1]
# @param gp


def g(rvar, rfix, beta, dbetadr):
    tmp = 1.             # [1]
    tmp -= 2*beta*rfix**2/rvar**2       # [1]
    tmp += beta*(beta+1.)/2.*rfix**4/rvar**4 # [1]
    tmp -= rfix**4/(4.*rvar**3)*dbetadr      # last term of eq. 37 Richardson2012
    gh.checknan(tmp, 'tmp in g')
    return tmp
## \fn g(rvar, rfix, beta, dbetadr)
# function from Lokas+2005, eq. (8)
# @param rvar small r, in [pc]
# @param rfix capital R, in [pc]
# @param beta beta, in [1]
# @param dbetadr d beta/dr, in [1]


def varepsilon(r0, betapar, gp):
    su = 0
    for k in np.arange(1, gp.nbeta):
        su += betapar[k]/k*r0**k
    return np.exp(-2*su)
## \fn varepsilon(r0, betapar, gp)
# return e^{-2 int_r_min^r beta(s)/s ds}, to be used in sigma_r^2 calc
# @param r0 radius [pc]
# @param betapar [1] coefficients of beta representation
# @param gp global parameters


def ant_sigkaplos(r0, rhopar, rhostarpar, MtoL, nupar, betapar, pop, gp):
    r0nu   = gp.xfine # [pc]
    rhonu  = phys.rho(r0nu,  rhopar, 0, gp) # DM mass profile (first)
    if gp.checksig:
        loglog(r0nu, ga.rho_gaia(r0nu, gp)[0], 'b.-')
        loglog(r0nu, rhonu, 'r.-')
        #xlabel('r/pc'); ylabel('rho')
        pdb.set_trace()

    # add up tracer densities to get overall density profile
    # add rho* to take into account the baryonic addition
    # (*not* Sigma from nu_i, could miss populations, have
    # varying selection function as fct of radius
    # need a M/L parameter)
    # only if we work on real data, add up total baryonic contribution
    if gp.investigate == 'obs':
        rhostarnu = MtoL*phys.nu(r0nu, rhostarpar, pop, gp)
        rhonu += rhostarnu

    #check influence of wrong beta
    #betapar[3] -= 0.1

    betanu = phys.beta(r0nu, betapar, gp)[0]
    if gp.checksig:
        clf()
        anbeta = ga.beta_gaia(r0nu, gp)[1]
        plot(r0nu, anbeta, 'b.-')
        plot(r0nu, betanu, 'r.-')
        #xlabel('r/pc'); ylabel('beta')
        pdb.set_trace()

    #nupar[0] *= 1e3 # TODO: check that siglos**2 is not changing with scaling of nu
    nunu   = phys.nu(r0nu, nupar, pop, gp)
    if gp.checksig:
        clf()
        loglog(gp.xipol, gp.dat.nu[pop], 'b.-')
        annu = ga.rho_gaia(r0nu, gp)[pop]
        loglog(r0nu, annu, 'g.-')
        loglog(r0nu, nunu, 'r.-')
        #xlabel('r/pc'); ylabel('nu')
        pdb.set_trace()

    Signu  = glp.rho_param_INT_Sig(r0nu, nupar, pop, gp)
    if gp.checksig:
        clf()
        loglog(gp.xipol, gp.dat.Sig[pop], 'b.-')
        loglog(r0nu[:-gp.nexp], Signu, 'r.-')
        #xlabel('r/pc'); ylabel('Sigma nu')
        pdb.set_trace()


    # int beta(s)/s ds
    # -----------------------------------------------------------------------
    idnu   = ant_intbeta(r0nu, betapar, gp)
    if gp.checksig:
        clf()
        plot(r0nu, idnu, 'r.-')
        beta_star1, r_DM, gamma_star1, r_star1, r_a1, gamma_DM, rho0 = gp.files.params
        anidnu = 0.5*(np.log(r0nu**2+r_a1**2)-np.log(r_a1**2))
        plot(r0nu, anidnu, 'b.-')
        xlabel('r/pc'); ylabel('int ds beta(s)/s')
        pdb.set_trace()

    # integrate enclosed 3D mass from 3D density
    r0tmp = np.hstack([0.,r0nu])


    # M(r)
    # ----------------------------------------------------------------------
    rhoint = 4.*np.pi*r0nu**2*rhonu # set to *nunu for check of integration
    # add point to avoid 0.0 in Mrnu(r0nu[0])
    rhotmp = np.hstack([0.,rhoint])
    splpar_rho = splrep(r0tmp, rhotmp, k=1, s=0.) # not necessarily monotonic
    Mrnu = np.zeros(len(r0nu))              # work in refined model
    for i in range(len(r0nu)):              # get Mrnu
        Mrnu[i] = splint(0., r0nu[i], splpar_rho)
    gh.checkpositive(Mrnu, 'Mrnu')
    if gp.checksig:
        loglog(gp.xipol, gp.dat.Mr[pop], 'g.-')
        s = r0nu/r_DM # [1]
        anMr = 4.*np.pi*rho0*r_DM**3*(1/(1+s)+np.log(1+s)-1) # [Msun]
        clf()
        loglog(r0nu, anMr, 'b.-')
        loglog(r0nu, Mrnu, 'r.-')
        xlabel('r/pc'); ylabel('M(r)')
        pdb.set_trace()



    # sig_r^2
    # --------------------------------------------------------------------------
    # (sigr2, 3D) * nu/exp(-idnu)
    xint = r0nu                           # [pc]
    yint = gp.G1 * Mrnu / r0nu**2         # [1/pc (km/s)^2]
    yint *= nunu                          # [Munit/pc^4 (km/s)^2]
    yint *= np.exp(2*idnu)                  # [Munit/pc^4 (km/s)^2]
    gh.checkpositive(yint, 'yint sigr2')
    if gp.checksig:
        clf()
        loglog(xint, yint, 'r.-')
        loglog(xint, gp.G1 * anMr / r0nu**2 * annu * np.exp(2*anidnu), 'b--')
        xlabel('xint/pc'); ylabel('yint')
        pdb.set_trace()


    # use quadinflog or quadinfloglog here
    sigr2nu = np.zeros(len(r0nu))
    for i in range(len(r0nu)):
        # TODO: check quadinflog with walker profiles
        # TODO: use same representation for splines here
        sigr2nu[i] = np.exp(-2*idnu[i])/nunu[i]*\
                     gh.quadinflog(xint, yint, r0nu[i], gp.rinfty*r0nu[-1], True)
        if sigr2nu[i] == np.inf:
            sigr2nu[i] = 1e-100
        # last arg: warn if no convergence found
    gh.checkpositive(sigr2nu, 'sigr2nu in sigl2s')
    if gp.checksig:
        clf()
        loglog(r0nu, sigr2nu, 'r.-')
        #xlabel('r/pc');ylabel('sigr2nu')
        pdb.set_trace()

    # project back to LOS values, \sigma_{LOS}^2 * \Sigma(R)
    # --------------------------------------------------------------
    sigl2s = np.zeros(len(r0nu)-gp.nexp)
    for i in range(len(r0nu)-gp.nexp): # get sig_los^2
        xnew = np.sqrt(r0nu[i:]**2-r0nu[i]**2)             # [pc]
        ynew = 2.*(1-betanu[i]*(r0nu[i]**2)/(r0nu[i:]**2)) # TODO check
        ynew *= nunu[i:] * sigr2nu[i:]
        gh.checkpositive(ynew, 'ynew in sigl2s') # is hit several times..
        # check sigr2nu: has too many entries of inf!

        # stop integration at xnew[-1] instead of at np.inf
        # to circumvent inf when nunu has increase at fudge radii
        sigl2s[i] = gh.quadinflog(xnew, ynew, 0, xnew[-1], False)
    # for last 3 bins, we are up to a factor 2 off
    gh.checkpositive(sigl2s, 'sigl2s')
    if gp.checksig:
        clf()
        loglog(r0nu[:-gp.nexp], sigl2s, 'r.-')
        #xlabel('r/pc'); ylabel('sigl2s')
        pdb.set_trace()

    siglos2 = sigl2s/Signu
    if gp.checksig:
        clf()
        plot(r0nu[:-gp.nexp], siglos2, 'r.-')
        #xlabel('r/pc'); ylabel('siglos2')
        pdb.set_trace()

    # derefine on radii of the input vector
    splpar_sig = splrep(r0nu[:-gp.nexp], np.log(siglos2), k=3, s=0.)
    siglos2_out = np.exp(splev(r0, splpar_sig))
    gh.checkpositive(siglos2_out, 'siglos2_out')
    if gp.checksig:
        clf()
        plot(r0, np.sqrt(siglos2_out), 'r.-', label='model')
        plot(gp.xipol, gp.dat.sig[pop], 'g.-', label='data')
        fill_between(gp.xipol, gp.dat.sig[pop]-gp.dat.sigerr[pop], gp.dat.sig[pop]+gp.dat.sigerr[pop], color='g', alpha=0.6)
        xscale('log'); xlim([1,2000])
        #xlabel('r/pc');ylabel('sigma LOS')
        pdb.set_trace()

    if not gp.usekappa:
        kapl4s_out = np.ones(len(siglos2_out))
    if gp.usekappa:
        kapl4s_out = kappa(r0nu, Mrnu, nunu, sigr2nu, idnu, gp)

    zetaa = -1; zetab = -1
    if gp.usezeta:
        zetaa, zetab = zeta(r0nu[:-gp.nexp], nunu[:-gp.nexp], \
                            Signu,\
                            Mrnu[:-gp.nexp], betanu[:-gp.nexp],\
                            sigr2nu[:-gp.nexp], gp)

    gh.sanitize_vector(siglos2_out, len(r0), 0, 1e30)
    return siglos2_out, kapl4s_out, zetaa, zetab
## \fn ant_sigkaplos(r0, rhopar, rhostarpar, MtoL, nupar, betapar, pop, gp)
# calculate integral for sig_los^2 * surface density,
# correspondingly for 4th order kappa,
# and virial parameters from Richardson,Fairbairn 2014
# all on gp.nfine logarithmically spaced radii,
# and store outcome on lookup table
# @param r0 array, radial bin positions, in [pc]
# @param rhopar overall 3D density profile in [Munit/pc^3]
# @param rhostarpar rho* light profile from baryons
# @param MtoL mass-to-light ratio in [Msun/Lsun]
# @param nupar tracer density in [Munit/pc^3]
# @param betapar velocity anisotropy in [1]
# @param pop int for population to look at (really the scalrad for rhohalf)
# @param gp global parameters
# @return sig_los^2, correspondingly for 4th order kappa, zetaa, zetab


def zeta(r0nu, nunu, Signu, Mrnu, betanu, sigr2nu, gp):
    # common parameters
    N = gh.Ntot(r0nu, Signu, gp)
    vr2 = sigr2nu
    dPhidr = gp.G1*Mrnu/r0nu**2

    # zetaa scalar
    xint = r0nu
    yint = nunu*(5-2*betanu)*vr2*dPhidr*r0nu**3
    nom = gh.quadinflog(xint, yint, 0., gp.rinfty*max(gp.xepol), False)

    yint = nunu*dPhidr*r0nu**3
    denom = (gh.quadinflog(xint, yint, 0., gp.rinfty*max(gp.xepol), False))**2

    zetaa = 9*N/10. * nom/denom

    # zetab scalar
    yint = nunu*(7-6*betanu)*vr2*dPhidr*r0nu**5
    nom = gh.quadinflog(xint, yint, 0., gp.rinfty*max(gp.xepol), False)

    yint = Signu*r0nu**3
    denom *= gh.quadinflog(xint, yint, 0., gp.rinfty*max(gp.xepol), False)

    zetab = 9*N**2/35 * nom / denom

    return zetaa, zetab
## \fn zeta(r0nu, nunu, Signu, Mrnu, betanu, sigr2nu, gp)
# return zeta_A and zeta_B parameters by Richardson,Fairbairn 2014
# @param r0nu radii, [pc]
# @param nunu 3D tracer density
# @param Signu 2D surface density
# @param Mrnu enclosed mass
# @param betanu anisotropy parameter
# @param sigr2nu <v_r**2>
# @param gp global parameter


def kappa(r0nu, Mrnu, nunu, sigr2nu, idnu, gp):
    # for the following: enabled calculation of kappa

    # kappa_r^4
    kapr4nu = np.ones(len(r0nu)-gp.nexp)
    xint  = r0nu                  # [pc]
    yint  = gp.G1 * Mrnu/r0nu**2  # [1/pc (km/s)^2]
    yint *= nunu                  # [Munit/pc^4 (km/s)^2]
    yint *= sigr2nu               # [Munit/pc^4 (km/s)^4
    yint *= np.exp(idnu)          # [Munit/pc^4 (km/s)^4]
    gh.checkpositive(yint, 'yint in kappa_r^4')
    yscale = 10.**(1.-min(np.log10(yint[1:])))
    yint *= yscale
    # power-law extrapolation to infinity
    C = max(0., gh.quadinflog(xint[-3:], yint[-3:], r0nu[-1], gp.rinfty*r0nu[-1]))

    splpar_nu = splrep(xint, yint, k=3) # interpolation in real space
    for i in range(len(r0nu)-gp.nexp):
        # integrate from minimal radius to infinity
        kapr4nu[i] = 3.*(np.exp(-idnu[i])/nunu[i]) * \
            (splint(r0nu[i], r0nu[-1], splpar_nu) + C) # [(km/s)^4]

    kapr4nu /= yscale
    gh.checkpositive(kapr4nu, 'kapr4nu in kappa_r^4')

    splpar_kap = splrep(r0nu[:-gp.nexp], np.log(kapr4nu), k=3)
    kapr4ext = np.exp(splev(r0ext, splpar_kap))
    kapr4nu = np.hstack([kapr4nu, kapr4ext])
    gh.checkpositive(kapr4nu, 'kapr4nu in extended kappa_r^4')

    dbetanudr = splev(r0nu, splrep(r0nu, betanu), der=1)
    gh.checknan(dbetanudr, 'dbetanudr in kappa_r^4')

    # kappa^4_los*surfdensity
    kapl4s = np.zeros(len(r0nu)-gp.nexp)
    for i in range(len(r0nu)-gp.nexp):
        xnew = np.sqrt(r0nu[i:]**2-r0nu[i]**2)      # [pc]
        ynew = g(r0nu[i:], r0nu[i], betanu[i:], dbetanudr[i:]) # [1]
        ynew *= nunu[i:] * kapr4nu[i:]
        C = max(0., gh.quadinflog(xnew[-3:], ynew[-3:], xnew[-1], gp.rinfty*xnew[-1]))
        splpar_nu = splrep(xnew,ynew) # not s=0.1, this sometimes gives negative entries after int
        kapl4s[i] = 2. * (splint(0., xnew[-1], splpar_nu) + C)
        #kapl4s[i] /= yscale
        # LOG('ynew = ',ynew,', kapl4s =', kapl4s[i])

    # TODO: sometimes the last value of kapl4s is nan: why?
    gh.checkpositive(kapl4s, 'kapl4s in kappa_r^4')

    # project kappa4_los as well
    # only use middle values to approximate, without errors in center and far
    kapl4s_out = np.exp(splev(r0, splrep(r0nu[4:-gp.nexp], kapl4s[4:], k=3))) # s=0.
    gh.checkpositive(kapl4s_out, 'kapl4s_out in kappa_r^4')
    return sigl2s_out, kapl4s_out


def smooth_ext(x0, y0, xext, k=1, s=0.1, log=True, slope=-3., minval=0.):
    if log:
        y0 = np.log(y0)
    tmp = splev(xext, splrep(x0, y0, k=k, s=s))

    # must impose maximum value: bounded by slope=slope and min
    if log:
        tmp = np.exp(tmp)
    if min(tmp) < minval:
        # LOG(2, 'extrapolation falling below ',minval,'!')
        for i in range(len(tmp)):
            if tmp[i]<0:
                tmp[i] = 1.e-30
    return tmp
## \fn smooth_ext(x0, y0, xext, k=1, s=0.1, log=True, slope=-3., minval=0.)
# extrapolate to high radii with linear extrapolation
# NOT USED ANYMORE
# @param x0 array, free variable
# @param y0 array, dependent variable
# @param xext array, points to evaluate function on
# @param k order of the spline, 1 as a default
# @param s smoothing. 0.01 is small
# @param log bool for working in logarithmic space
# @param slope without upper limit
# @param minval minimal value allowed


def smooth_fun(x0, y0, xnew, k=3, s=0.01, log=True):
    if log:
        y0 = np.log(y0)
    tmp = splev(xnew, splrep(x0, y0, k=k, s=s))
    if log:
        tmp = np.exp(tmp)
    return tmp
## \fn smooth_fun(x0, y0, xnew, k=3, s=0.01, log=True)
# smooth function y(x)
# NOT USED ANYMORE
# @param x0 array, free variable
# @param y0 array, dependent variable
# @param xnew array, new positions (or the same as x0)
# @param k order of the spline
# @param s smoothing. 0.01 is *small*
# @param log bool for working in logarithmic space
