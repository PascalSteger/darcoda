#!/usr/bin/env ipython3

##
# @file
# all integrals from gl_physics

# (c) 2013 Pascal S.P. Steger

import numpy as np
import pdb, scipy, time
from scipy.integrate import simps,trapz,quad,romberg
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
# integrate polynomial from r0[k] to infinity
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
    # yint =  ga.beta(xint, gp)[1]/xint

    intbet = np.zeros(len(r0))
    for k in range(len(r0)):
        intbet[k] = gh.quadinf(xint, yint, r0[0]/1e5, r0[k])
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

    rmin = np.log10(min(r0))
    rmax = np.log10(max(r0)*gp.rinfty)
    r0fine = np.logspace(rmin, rmax, gp.nfine)

    # rho
    # --------------------------------------------------------------------------
    rhofine  = phys.rho(r0fine,  rhopar, 0, gp) # DM mass profile (first)
    if gp.checksig and gp.stopstep <= 1:
        clf()
        loglog(r0fine, rhofine, 'r.-', label='rederived from dn/dlogr params')
        loglog(r0fine, ga.rho(r0fine, gp)[0], 'b--', label='analytic')
        xlabel('$r/\\rm{pc}$')
        ylabel('$\\rho(r)$')
        legend(loc='lower left')
        savefig('fit_rho_hern.png')
        pdb.set_trace()

    #N = 100
    #start = time.time()
    #elapsed = (time.time()-start)/N
    #print('one iteration dsqrt takes ', elapsed, 's')
    #for k in range(N):
    #Sigfine_dsqrt  = glp.rho_param_INT_Sig(r0fine, rhopar, 0, gp)
    #elapsed = (time.time()-start)/N
    #print('one iteration dsqrt takes ', elapsed, 's')

    #start = time.time()
    #for k in range(N):
    #Sigfine_dtheta  = glp.rho_param_INT_Sig_theta(r0fine, rhopar, 0, gp)
    #elapsed = (time.time()-start)/N
    #print('one iteration dtheta takes ', elapsed, 's')

    #if gp.checksig:
    #    clf()
    #    anSig = ga.Sigma(r0fine, gp)
    #    loglog(r0fine, anSig, 'b-', label='analytic')
    #    loglog(r0fine, Sigfine_dsqrt, 'r.', label='dsqrt')
    #    loglog(r0fine, Sigfine_dtheta, 'g.', label='dtheta')
    #    xlabel('$r/\\rm{pc}$')
    #    ylabel('$\\Sigma/\\Sigma_{\\rm{analytic}}$')
    #    legend(loc='lower left')
    #    savefig('fit_Sig_hern.png')
    #    pdb.set_trace()

    # add up tracer densities to get overall density profile
    # add rho* to take into account the baryonic addition
    # (*not* Sigma from nu_i, could miss populations, have
    # varying selection function as fct of radius
    # need a M/L parameter)
    # only if we work on real data, add up total baryonic contribution
    if gp.investigate == 'obs':
        nu_baryons = MtoL*phys.nu(r0fine, rhostarpar, pop, gp)
        rhofine += nu_baryons

    # TODO: check influence of wrong beta
    # betapar[3] -= 0.1


    # beta
    # ------------------------------------------------------------------------
    betafine = phys.beta(r0fine, betapar, gp)[0]
    if gp.checksig and gp.stopstep <= 2:
        clf()
        anbeta = ga.beta(r0fine, gp)[1]
        plot(r0fine, betafine, 'r.-', label='model')
        plot(r0fine, anbeta, 'b--', label='analytic')
        xlabel('$r/\\rm{pc}$')
        ylabel('$\\beta$')
        ylim([-0.5, 1.0])
        legend(loc='lower right')
        savefig('fit_beta_hern.png')
        pdb.set_trace()

    # nu
    # ------------------------------------------------------------------------
    #nupar[0] *= 1e3 # TODO: check that siglos**2 is not changing with scaling of nu
    nufine   = phys.nu(r0fine, nupar, pop, gp)
    if gp.checksig:
        annu = ga.rho(r0fine, gp)[pop]

    if gp.checksig and gp.stopstep <= 3:
        clf()
        loglog(gp.xipol, gp.dat.nu[pop], 'g.-', label='data')
        fill_between(gp.xipol, gp.dat.nu[pop]-gp.dat.nuerr[pop], \
                     gp.dat.nu[pop]+gp.dat.nuerr[pop],\
                     color='g', alpha=0.6)
        loglog(r0fine, nufine, 'r.-', label='model')
        loglog(r0fine, annu, 'b--', label='analytic')
        legend(loc='lower left')
        xlabel('$r/\\rm{pc}$')
        ylabel('$\\nu$')
        savefig('fit_nu_hern.png')
        pdb.set_trace()

    # \Sigma
    # -----------------------------------------------------------------------
    Sigfine  = glp.rho_param_INT_Sig_theta(r0fine, nupar, pop, gp)
    if gp.checksig and gp.stopstep <= 4:
        clf()
        loglog(gp.xipol, gp.dat.Sig[pop], 'g--', label='data')

        loglog(r0fine, Sigfine, 'r.-', label='model')
        anSig = ga.Sigma(r0fine, gp)[pop]
        loglog(r0fine, anSig, 'b--', label='analytic')
        fill_between(gp.xipol, gp.dat.Sig[pop]-gp.dat.Sigerr[pop], \
                     gp.dat.Sig[pop]+gp.dat.Sigerr[pop],\
                     color='g', alpha=0.6)
        xlabel('$r/\\rm{pc}$')
        ylabel('$\\Sigma$')
        legend(loc='lower left')
        savefig('fit_Sig_hern.png')
        pdb.set_trace()

    # int beta(s)/s ds
    # -----------------------------------------------------------------------
    intbetasfine   = ant_intbeta(r0fine, betapar, gp)
    if gp.checksig:
        #beta_star1, r_DM, gamma_star1, r_star1, r_a1, gamma_DM, rho0 = gp.files.params
        # anintbetasfine = 0.5*(np.log(r0fine**2+r_a1**2)-np.log(r_a1**2)) # gaia
        anintbetasfine = 0.0*r0fine

    if gp.checksig and gp.stopstep <= 5 :
        clf()
        plot(r0fine, intbetasfine, 'r.-', label='model')
        plot(r0fine, anintbetasfine, 'b--', label='analytic')
        xlabel('$r/\\rm{pc}$')
        ylabel('$\\int ds \\beta(s)/s$')
        legend(loc='lower right')
        savefig('fit_intbeta_hern.png')
        pdb.set_trace()



    # M(r)
    # ----------------------------------------------------------------------
    rhoint = 4.*np.pi*r0fine**2*rhofine

    # add point to avoid 0.0 in Mrfine(r0fine[0])
    r0tmp = np.hstack([0.,r0fine])
    rhotmp = np.hstack([0.,rhoint])
    splpar_rho = splrep(r0tmp, rhotmp, k=1, s=0.) # not necessarily monotonic
    Mrfine = np.zeros(len(r0fine)) # work in refined model
    for i in range(len(r0fine)):
        Mrfine[i] = splint(0., r0fine[i], splpar_rho)
    gh.checkpositive(Mrfine, 'Mrfine')
    if gp.checksig:
        anMr = ga.Mr(r0fine, gp)[pop]

    if gp.checksig and gp.stopstep <= 6:
        #loglog(gp.xipol, gp.dat.Mr[pop], 'g.-', label='data')
        #s = r0fine/r_DM # [1]
        clf()
        loglog(r0fine, Mrfine, 'r.-', label='model')
        loglog(r0fine, anMr, 'b--', label='analytic')
        xlabel('$r/\\rm{pc}$')
        ylabel('$M(r)$')
        legend(loc='lower right')
        savefig('fit_M_hern_better.png')
        pdb.set_trace()


    # nu(r)\cdot\sigma_r^2(r) integrand
    # --------------------------------------------------------------------------
    # (sigr2, 3D) * nu/exp(-intbetasfine)
    xint = r0fine                           # [pc]
    yint = gp.G1 * Mrfine / r0fine**2         # [1/pc (km/s)^2]
    yint *= nufine                          # [Munit/pc^4 (km/s)^2]
    yint *= np.exp(2*intbetasfine)                  # [Munit/pc^4 (km/s)^2]
    gh.checkpositive(yint, 'yint sigr2')
    if gp.checksig and gp.stopstep <= 7:
        clf()
        loglog(xint, yint, 'r.-', label='model')
        loglog(xint, gp.G1 * anMr / r0fine**2 * annu * np.exp(2*anintbetasfine), 'b--', label='from analytic')
        xlabel('$xint/\\rm{pc}$')
        ylabel('$yint$')
        legend(loc='lower left')
        savefig('fit_nu_sigmar2_hern.png')
        pdb.set_trace()


    # actual integration
    #
    #cth = cos(theta)
    #sth = sin(theta)
    #for i in range(len(r)):
    #	rq = r[i] / cth
    #	sigmar2[i] =  intbeta(r[i]) / (r[i] * nu(r[i])) * \
    #		integrator(G*M(rq)*nu(rq)*intbeta(rq)*sth,theta)

    # --------------------------------------------------------------------------
    # use quadinflog or quadinfloglog here
    sigr2model = np.zeros(len(r0fine))
    for k in range(len(r0fine)):
        # log space
        minlog = min(np.log(yint[yint>0])) # exclude y==0 to circumvent error
        shiftlog = np.exp(2-minlog)
        # otherwise we get divergent integrals for high radii
        yshift = yint * shiftlog
        splpar_nul = splrep(xint, np.log(yshift), k=1, s=0.)
        invexp = lambda x: np.exp(splev(x, splpar_nul))
        integ = romberg(invexp, r0fine[k], gp.rinfty*r0fine[-1],\
                             rtol=1e-3, divmax=15, vec_func=True)
        integ /= shiftlog

        # loglog
        #splpar_nul = splrep(np.log(xint), np.log(yint), k=1, s=0.1) # tunable k=2; s=0, ..
        #invexp = lambda x: np.exp(splev(np.log(x), splpar_nul))
        #integ = romberg(invexp, r0fine[k], 5*r0fine[-1],\
        #                rtol=1e-3, divmax=15, vec_func=True)
        sigr2model[k] = np.exp(-2*intbetasfine[k])/nufine[k]*integ

        # old function quadinflog
        #sigr2model[k] = np.exp(-2*intbetasfine[k])/nufine[k]*\
        #                gh.quadinflog(xint, yint, r0fine[k], gp.rinfty*r0fine[-1], True)

        if sigr2model[k] == np.inf:
            sigr2model[k] = 1e-100
    gh.checkpositive(sigr2model, 'sigr2model in sigl2s')
    if gp.checksig and gp.stopstep <= 8:
        clf()
        ansigr2 = ga.sigr2(r0fine, gp)
        loglog(r0fine, sigr2model, 'r.-', label='model')
        loglog(r0fine, ansigr2, 'b--', label='analytic')
        xlabel('$r/\\rm{pc}$')
        ylabel('$\\sigma_r^2(r)\\nu(r)$')
        legend(loc='lower right')
        savefig('fit_sigr2_hern.png')
        pdb.set_trace()


    # project back to LOS values, \sigma_{LOS}^2 * \Sigma(R)
    # -------------------------------------------------------------
    start = time.time()
    sigl2s_dsqrt = np.zeros(len(r0fine)-gp.nexp)             #
    for k in range(len(r0fine)-gp.nexp): # get sig_los^2
        xnew = np.sqrt(r0fine[k:]**2-r0fine[k]**2)             # [pc]
        ynew = 2.*(1-betafine[k]*(r0fine[k]**2)/(r0fine[k:]**2)) # TODO check
        ynew *= nufine[k:] * sigr2model[k:]
        gh.checkpositive(ynew, 'ynew in sigl2s') # is hit several times..
        # check sigr2model: has too many entries of inf!

        # stop integration at xnew[-1] instead of at np.inf
        # to circumvent inf when nufine has increase at fudge radii
        sigl2s_dsqrt[k] = gh.quadinflog(xnew, ynew, 0, xnew[-1], False)
    elapsed = time.time()-start
    print('one iteration dsqrt takes ', elapsed, 's')

    start = time.time()
    sigl2s_dtheta = np.zeros(len(r0fine))
    xmin = gp.xfine[0]/1. # needed, if not: loose on first 4 bins
    r0fine = gp.xfine

    bit = 1.e-6
    theta = np.linspace(0, np.pi/2-bit, gp.nfine)
    cth = np.cos(theta)
    cth2 = cth*cth
    Rproj = 1.*r0fine

    ynew = (1-betafine*cth2)*sigr2model*nufine # is fine for Hernquist
    # for debugging purposes of the theta method
    #ynew = (1-anbeta*cth2)*ansigr2*annu
    for k in range(len(r0fine)):
        rq = Rproj[k]/cth
        ynewq = np.interp(rq, r0fine, ynew, left=0, right=0)
        sigl2s_dtheta[k] = 2.*Rproj[k]*simps(ynewq/cth2, theta)
    elapsed = time.time()-start
    print('one iteration dtheta takes ', elapsed, 's')


    #sigl2s = sigl2s_dtheta[:-gp.nexp]
    # for last 3 bins, we are up to a factor 2 off
    gh.checkpositive(sigl2s_dtheta, 'sigl2s_dtheta')
    gh.checkpositive(sigl2s_dsqrt, 'sigl2s_dsqrt')
    if gp.checksig and gp.stopstep <= 9:
        clf()
        anSigsiglos2_hern = ga.Sig_sig_los_2(r0fine, gp)
        loglog(r0fine[:-gp.nexp], sigl2s_dsqrt, 'g.-', label='dsqrt')
        loglog(r0fine, sigl2s_dtheta, 'r.-', label='dtheta')
        loglog(r0fine, anSigsiglos2_hern, 'b--', label='analytic')
        xlabel('$r/\\rm{pc}$')
        ylabel('$\\sigma_{\\rm{LOS}}^2 \Sigma$')
        legend(loc='lower left')
        savefig('fit_Sig_siglos2_hern.png')
        pdb.set_trace()


    # sigma_LOS^2
    # ----------------------------------------------------------------------
    siglos2_dtheta = sigl2s_dtheta/Sigfine
    siglos2_dsqrt = sigl2s_dsqrt/Sigfine[:-gp.nexp]
    if gp.checksig and gp.stopstep <= 10:
        clf()
        ansiglos = ga.sig_los(r0fine, gp)
        plot(r0fine[:-gp.nexp], siglos2_dsqrt, 'g.-', label='dsqrt')
        plot(r0fine, siglos2_dtheta, 'r.-', label='dtheta')
        plot(r0fine, ansiglos**2, 'b--', label='analytic')
        xlabel('$r/\\rm{pc}$')
        ylabel('$\\sigma_{\\rm{LOS}}^2$')
        legend(loc='upper right')
        savefig('fit_siglos2_hern.png')
        pdb.set_trace()

    # derefine on radii of the input vector
    splpar_sig_dtheta = splrep(r0fine, np.log(siglos2_dtheta), k=3, s=0.)
    splpar_sig_dsqrt = splrep(r0fine[:-gp.nexp], np.log(siglos2_dsqrt), k=3, s=0.)
    siglos2_out_dtheta = np.exp(splev(r0, splpar_sig_dtheta))
    siglos2_out_dsqrt = np.exp(splev(r0, splpar_sig_dsqrt))
    gh.checkpositive(siglos2_out_dtheta, 'siglos2_out_dtheta')
    gh.checkpositive(siglos2_out_dsqrt, 'siglos2_out_dsqrt')
    if gp.checksig and gp.stopstep <= 11:
        clf()
        ansiglos = ga.sig_los(r0, gp)
        plot(r0, np.sqrt(siglos2_out_dsqrt), 'g.-', label='dsqrt')
        plot(r0, np.sqrt(siglos2_out_dtheta), 'r.-', label='dtheta')
        plot(r0, ansiglos, 'b--', label='analytic')
        plot(gp.xipol, gp.dat.sig[pop], 'g.-', label='data')
        fill_between(gp.xipol, gp.dat.sig[pop]-gp.dat.sigerr[pop], gp.dat.sig[pop]+gp.dat.sigerr[pop], color='g', alpha=0.6)
        xscale('log')
        xlabel('$r/\\rm{pc}$')
        ylabel('$\\sigma_{\\rm{LOS}}$')
        legend(loc='upper right')
        savefig('fit_siglos_out_hern.png')
        pdb.set_trace()
    siglos2_out = siglos2_out_dtheta

    if not gp.usekappa:
        kapl4s_out = np.ones(len(siglos2_out))
    if gp.usekappa:
        kapl4s_out = kappa(r0fine, Mrfine, nufine, sigr2model, intbetasfine, gp)

    zetaa = -1; zetab = -1
    if gp.usezeta:
        zetaa, zetab = zeta(r0fine[:-gp.nexp], nufine[:-gp.nexp], \
                            Sigfine,\
                            Mrfine[:-gp.nexp], betafine[:-gp.nexp],\
                            sigr2model[:-gp.nexp], gp)

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


def zeta(r0fine, nufine, Sigfine, Mrfine, betafine, sigr2nu, gp):
    # common parameters
    N = gh.Ntot(r0fine, Sigfine, gp)
    vr2 = sigr2nu
    dPhidr = gp.G1*Mrfine/r0fine**2

    # zetaa scalar
    xint = r0fine
    yint = nufine*(5-2*betafine)*vr2*dPhidr*r0fine**3
    nom = gh.quadinflog(xint, yint, 0., gp.rinfty*max(gp.xepol), False)

    yint = nufine*dPhidr*r0fine**3
    denom = (gh.quadinflog(xint, yint, 0., gp.rinfty*max(gp.xepol), False))**2

    zetaa = 9*N/10. * nom/denom

    # zetab scalar
    yint = nufine*(7-6*betafine)*vr2*dPhidr*r0fine**5
    nom = gh.quadinflog(xint, yint, 0., gp.rinfty*max(gp.xepol), False)

    yint = Sigfine*r0fine**3
    denom *= gh.quadinflog(xint, yint, 0., gp.rinfty*max(gp.xepol), False)

    zetab = 9*N**2/35 * nom / denom

    return zetaa, zetab
## \fn zeta(r0fine, nufine, Sigfine, Mrfine, betafine, sigr2nu, gp)
# return zeta_A and zeta_B parameters by Richardson,Fairbairn 2014
# @param r0fine radii, [pc]
# @param nufine 3D tracer density
# @param Sigfine 2D surface density
# @param Mrfine enclosed mass
# @param betafine anisotropy parameter
# @param sigr2nu <v_r**2>
# @param gp global parameter


def kappa(r0fine, Mrfine, nufine, sigr2nu, intbetasfine, gp):
    # for the following: enabled calculation of kappa

    # kappa_r^4
    kapr4nu = np.ones(len(r0fine)-gp.nexp)
    xint  = r0fine                  # [pc]
    yint  = gp.G1 * Mrfine/r0fine**2  # [1/pc (km/s)^2]
    yint *= nufine                  # [Munit/pc^4 (km/s)^2]
    yint *= sigr2nu               # [Munit/pc^4 (km/s)^4
    yint *= np.exp(intbetasfine)          # [Munit/pc^4 (km/s)^4]
    gh.checkpositive(yint, 'yint in kappa_r^4')
    yscale = 10.**(1.-min(np.log10(yint[1:])))
    yint *= yscale
    # power-law extrapolation to infinity
    C = max(0., gh.quadinflog(xint[-3:], yint[-3:], r0fine[-1], gp.rinfty*r0fine[-1]))

    splpar_nu = splrep(xint, yint, k=3) # interpolation in real space
    for k in range(len(r0fine)-gp.nexp):
        # integrate from minimal radius to infinity
        kapr4nu[k] = 3.*(np.exp(-intbetasfine[k])/nufine[k]) * \
            (splint(r0fine[k], r0fine[-1], splpar_nu) + C) # [(km/s)^4]

    kapr4nu /= yscale
    gh.checkpositive(kapr4nu, 'kapr4nu in kappa_r^4')

    splpar_kap = splrep(r0fine[:-gp.nexp], np.log(kapr4nu), k=3)
    kapr4ext = np.exp(splev(r0ext, splpar_kap))
    kapr4nu = np.hstack([kapr4nu, kapr4ext])
    gh.checkpositive(kapr4nu, 'kapr4nu in extended kappa_r^4')

    dbetafinedr = splev(r0fine, splrep(r0fine, betafine), der=1)
    gh.checknan(dbetafinedr, 'dbetafinedr in kappa_r^4')

    # kappa^4_los*surfdensity
    kapl4s = np.zeros(len(r0fine)-gp.nexp)
    for k in range(len(r0fine)-gp.nexp):
        xnew = np.sqrt(r0fine[k:]**2-r0fine[k]**2)      # [pc]
        ynew = g(r0fine[k:], r0fine[k], betafine[k:], dbetafinedr[k:]) # [1]
        ynew *= nufine[k:] * kapr4nu[k:]
        C = max(0., gh.quadinflog(xnew[-3:], ynew[-3:], xnew[-1], gp.rinfty*xnew[-1]))
        splpar_nu = splrep(xnew,ynew) # not s=0.1, this sometimes gives negative entries after int
        kapl4s[k] = 2. * (splint(0., xnew[-1], splpar_nu) + C)
        #kapl4s[k] /= yscale
        # LOG('ynew = ',ynew,', kapl4s =', kapl4s[k])

    # TODO: sometimes the last value of kapl4s is nan: why?
    gh.checkpositive(kapl4s, 'kapl4s in kappa_r^4')

    # project kappa4_los as well
    # only use middle values to approximate, without errors in center and far
    kapl4s_out = np.exp(splev(r0, splrep(r0fine[4:-gp.nexp], kapl4s[4:], k=3))) # s=0.
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
        for k in range(len(tmp)):
            if tmp[k]<0:
                tmp[k] = 1.e-30
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
