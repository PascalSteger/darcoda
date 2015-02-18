#!/usr/bin/env ipython3

##
# @file
# all integrals from gi_physics

# (c) GPL v3 2015 Pascal S.P. Steger pascal@steger.aero

import numpy as np
import pdb
from scipy.integrate import simps,quad
from scipy.interpolate import splrep, splev, splint

import gi_units as gu
import gi_helper as gh
import gi_analytic as ga
import gi_physics as phys
import gi_project as gip

import matplotlib
matplotlib.use('pdf')
from pylab import *
ion()

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
    # assumption here is that integration goes down to min(r0)/1e5

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

def ant_sigkaplos(r0, rhodmpar, lbaryonpar, MtoL, nupar, betapar, pop, gp):
    rmin = np.log10(min(r0))
    rmax = np.log10(max(r0)*gp.rinfty)
    r0fine = np.logspace(rmin, rmax, gp.nfine)

    # rho
    # --------------------------------------------------------------------------
    rhofine  = phys.rho(r0fine,  rhodmpar, 0, gp) # DM mass profile (first)
    if gp.checksig and gp.stopstep <= 1:
        clf()
        loglog(r0fine, rhofine, 'r.-', label='rederived from dn/dlogr params')
        loglog(r0fine, ga.rho(r0fine, gp)[0], 'b--', label='analytic')
        axvline(max(gp.xipol))
        axvline(min(gp.xipol))
        axvline(gp.dat.rhalf[0], lw=2)
        xlabel('$r/\\rm{pc}$')
        ylabel('$\\rho(r)$')
        legend(loc='lower left')
        savefig('fit_rho_'+gp.investigate+'.pdf')
        pdb.set_trace()

    # add up tracer densities to get overall density profile
    # add rho* to take into account the baryonic addition
    # (*not* Sigma from nu_i, could miss populations, have
    # varying selection function as fct of radius
    # need a M/L parameter)
    # only if we work on real data, add up total baryonic contribution
    if gp.investigate == 'obs':
        nu_baryons = MtoL*phys.nu(r0fine, lbaryonpar, pop, gp)
        rhofine += nu_baryons


    # beta
    # ------------------------------------------------------------------------
    betafine = phys.beta(r0fine, betapar, gp)[0]
    if gp.checksig and gp.stopstep <= 2:
        clf()
        anbeta = ga.beta(r0fine, gp)[1]
        plot(r0fine, betafine, 'r.-', label='model')
        plot(r0fine, anbeta, 'b--', label='analytic')
        xscale('log')
        axvline(max(gp.xipol))
        axvline(min(gp.xipol))
        axvline(gp.dat.rhalf[0], lw=2)
        xlabel('$r/\\rm{pc}$')
        ylabel('$\\beta$')
        ylim([-0.5, 1.0])
        legend(loc='lower right')
        savefig('fit_beta_'+gp.investigate+'.pdf')
        pdb.set_trace()

    # nu
    # ------------------------------------------------------------------------
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
        axvline(max(gp.xipol))
        axvline(min(gp.xipol))
        axvline(gp.dat.rhalf[0], lw=2)
        xlabel('$r/\\rm{pc}$')
        ylabel('$\\nu$')
        savefig('fit_nu_'+gp.investigate+'.pdf')
        pdb.set_trace()


    # \Sigma
    # ---------------------------------------------------------------
    Sigfine  = gip.rho_param_INT_Sig_theta(r0fine, nupar, pop, gp)
    if gp.checksig and gp.stopstep <= 4:
        clf()
        anSig = ga.Sigma(r0fine, gp)[pop]
        loglog(gp.xipol, gp.dat.Sig[pop], 'g--', label='data')
        loglog(r0fine, Sigfine, 'r.-', label='model')
        loglog(r0fine, anSig, 'b--', label='analytic')
        fill_between(gp.xipol, gp.dat.Sig[pop]-gp.dat.Sigerr[pop], \
                     gp.dat.Sig[pop]+gp.dat.Sigerr[pop],\
                     color='g', alpha=0.6)
        axvline(max(gp.xipol))
        axvline(min(gp.xipol))
        axvline(gp.Xscale[0], lw=2)
        xlabel('$r/\\rm{pc}$')
        ylabel('$\\Sigma$')
        legend(loc='lower left')
        savefig('fit_Sig_'+gp.investigate+'.pdf')
        pdb.set_trace()


    # int beta(s)/s ds
    # ---------------------------------------------------------------
    intbetasfine   = ant_intbeta(r0fine, betapar, gp)
    if gp.checksig:
        if gp.investigate == 'gaia':
            beta_star1, r_DM, gamma_star1, r_star1, r_a1, gamma_DM, rho0 = gp.files.params
            anintbetasfine = 0.5*(np.log(r0fine**2+r_a1**2)-np.log(r_a1**2))
        elif gp.investigate == 'hern':
            anintbetasfine = 0.0*r0fine
    if gp.checksig and gp.stopstep <= 5 :
        clf()
        plot(r0fine, intbetasfine, 'r.-', label='model')
        plot(r0fine, anintbetasfine, 'b--', label='analytic')
        ylim([-5, 5])
        axvline(max(gp.xipol))
        axvline(min(gp.xipol))
        axvline(gp.dat.rhalf[0], lw=2)
        xscale('log')
        xlabel('$r/\\rm{pc}$')
        ylabel('$\\int ds \\beta(s)/s$')
        legend(loc='lower right')
        savefig('fit_intbeta_'+gp.investigate+'.pdf')
        pdb.set_trace()


    # M(r)
    # ---------------------------------------------------------------
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
        axvline(max(gp.xipol))
        axvline(min(gp.xipol))
        axvline(gp.dat.rhalf[0], lw=2)
        xlabel('$r/\\rm{pc}$')
        ylabel('$M(r)$')
        legend(loc='lower right')
        savefig('fit_M_hern_better.pdf')
        pdb.set_trace()


    # nu(r)\cdot\sigma_r^2(r) integrand
    # --------------------------------------------------------------------------
    # (sigr2, 3D) * nu/exp(-intbetasfine)
    xint = r0fine                           # [pc]
    yint = gu.G1__pcMsun_1km2s_2 * Mrfine / r0fine**2         # [1/pc (km/s)^2]
    yint *= nufine                          # [Munit/pc^4 (km/s)^2]
    yint *= np.exp(2*intbetasfine)                  # [Munit/pc^4 (km/s)^2]
    gh.checkpositive(yint, 'yint sigr2')
    if gp.checksig and gp.stopstep <= 7:
        clf()
        loglog(xint, yint, 'r.-', label='model')
        loglog(xint, gu.G1__pcMsun_1km2s_2 * anMr / r0fine**2 * annu * np.exp(2*anintbetasfine), 'b--', label='from analytic')
        axvline(max(gp.xipol))
        axvline(min(gp.xipol))
        axvline(gp.dat.rhalf[0], lw=2)
        xlabel('$xint/\\rm{pc}$')
        ylabel('$yint$')
        legend(loc='lower left')
        savefig('fit_nu_sigmar2_'+gp.investigate+'.pdf')
        pdb.set_trace()

    # actual integration, gives \sigma_r^2 \nu
    sigr2nu_model = np.zeros(len(r0fine))
    sigr2nu_model_new = np.zeros(len(r0fine))
    for k in range(len(r0fine)):
        #theta_old = np.linspace(0, np.arccos(r0fine[k]/(gp.rinfty*max(gp.xepol))), gp.nfine)
        theta = np.arccos(r0fine[k]/r0fine[k:])
        rq = r0fine[k]/np.cos(theta)

        Mrq = np.interp(rq, r0fine, Mrfine, left=0, right=0)
        nuq = np.interp(rq, r0fine, nufine, left=0, right=0)
        intbetaq = np.interp(rq, r0fine, intbetasfine, left=0, right=0)
        func_interp_before = Mrq*nuq*np.exp(2*intbetaq)

        func_base = Mrfine*nufine*np.exp(2*intbetasfine)
        #func_interp_after = np.interp(rq, r0fine, func_base, left=0, right=0)
        func_interp_after = func_base[k:]

        #print('median(func_interp_after / func_interp_before = ',\
        #      np.median(func_interp_after / func_interp_before))

        sigr2nu_model[k] =  np.exp(-2*intbetasfine[k])/r0fine[k] * \
                            gu.G1__pcMsun_1km2s_2*simps(func_interp_before*np.sin(theta), theta)

        sigr2nu_model_new[k] = np.exp(-2*intbetasfine[k])/r0fine[k] * \
                               gu.G1__pcMsun_1km2s_2*simps(func_interp_after*np.sin(theta), theta)

    # clean last value (which is always 0 by construction)
    sigr2nu_model[-1] = sigr2nu_model[-2]/10.
    sigr2nu_model_new[-1] = sigr2nu_model_new[-2]/10.
    gh.checkpositive(sigr2nu_model, 'sigr2nu_model in sigl2s')
    gh.checkpositive(sigr2nu_model_new, 'sigr2nu_model_new in sigl2s')
    if gp.checksig and gp.stopstep <= 8:
        clf()
        ansigr2nu = ga.sigr2(r0fine, gp)*annu
        loglog(r0fine, sigr2nu_model, 'r.-', label='model, interp each function')
        loglog(r0fine, sigr2nu_model_new, 'k--', label='model, interp product')
        loglog(r0fine, ansigr2nu, 'b--', label='analytic')
        axvline(max(gp.xipol))
        axvline(min(gp.xipol))
        axvline(gp.dat.rhalf[0], lw=2)
        xlabel('$r/\\rm{pc}$')
        ylabel('$\\sigma_r^2(r)\\nu(r)$')
        legend(loc='lower right')
        savefig('fit_sigr2_'+gp.investigate+'.pdf')

    sigr2nu_model = sigr2nu_model_new

    # project back to LOS values, \sigma_{LOS}^2 * \Sigma(R)
    # -------------------------------------------------------------
    sigl2s = np.zeros(len(r0fine))
    for k in range(len(r0fine)):
        bit = 1.e-6
        theta = np.linspace(0, np.pi/2-bit, gp.nfine)
        # work on same radii as data are given
        theta = np.arccos(r0fine[k]/r0fine[k:])
        cth = np.cos(theta)
        cth2 = cth*cth
        ynew = (1-betafine[k:]*cth2)*sigr2nu_model[k:]

        rq = r0fine[k]/cth
        ynewq = np.interp(rq, r0fine[k:], ynew, left=0, right=0)
        sigl2s[k] = 2.*r0fine[k]*simps(ynewq/cth2, theta)
    sigl2s[-1] = sigl2s[-2]/10.
    gh.checkpositive(sigl2s, 'sigl2s')
    if gp.checksig and gp.stopstep <= 9:
        clf()
        anSigsiglos2_hern = ga.Sig_sig_los_2(r0fine, gp)
        loglog(r0fine, sigl2s, 'r.-', label='model')
        loglog(r0fine, anSigsiglos2_hern, 'b--', label='analytic')
        axvline(max(gp.xipol))
        axvline(min(gp.xipol))
        axvline(gp.Xscale[0], lw=2)
        xlabel('$r/\\rm{pc}$')
        ylabel('$\\sigma_{\\rm{LOS}}^2 \Sigma$')
        legend(loc='lower left')
        savefig('fit_Sig_siglos2_'+gp.investigate+'.pdf')
        pdb.set_trace()


    # sigma_LOS^2
    # ----------------------------------------------------------------------
    siglos2 = sigl2s/Sigfine
    if gp.checksig and gp.stopstep <= 10:
        clf()
        ansiglos = ga.sig_los(r0fine, gp)
        plot(r0fine, siglos2, 'r.-', label='model')
        plot(r0fine, ansiglos**2, 'b--', label='analytic')
        axvline(max(gp.xipol))
        axvline(min(gp.xipol))
        axvline(gp.Xscale[0], lw=2)
        xscale('log')
        xlabel('$r/\\rm{pc}$')
        ylabel('$\\sigma_{\\rm{LOS}}^2$')
        legend(loc='upper right')
        savefig('fit_siglos2_'+gp.investigate+'.pdf')
        pdb.set_trace()

    # derefine on radii of the input vector
    splpar_sig = splrep(r0fine, np.log(siglos2), k=3, s=0.)
    siglos2_out = np.exp(splev(r0, splpar_sig))
    # gh.checkpositive(siglos2_out, 'siglos2_out')
    if gp.checksig and gp.stopstep <= 11:
        clf()
        ansiglos = ga.sig_los(r0, gp)
        plot(r0, np.sqrt(siglos2_out), 'r.-', label='model')
        plot(r0, ansiglos, 'b--', label='analytic')
        plot(gp.xipol, gp.dat.sig[pop], 'g.-', label='data')
        fill_between(gp.xipol, gp.dat.sig[pop]-gp.dat.sigerr[pop], gp.dat.sig[pop]+gp.dat.sigerr[pop], color='g', alpha=0.6)
        xscale('log')
        axvline(max(gp.xipol))
        axvline(min(gp.xipol))
        axvline(gp.Xscale[0], lw=2)
        xlabel('$r/\\rm{pc}$')
        ylabel('$\\sigma_{\\rm{LOS}}$')
        legend(loc='upper right')
        savefig('fit_siglos_out_'+gp.investigate+'.pdf')
        pdb.set_trace()

    if not gp.usekappa:
        kapl4s_out = np.ones(len(siglos2_out))
    if gp.usekappa:
        kapl4s_out = kappa(r0fine, Mrfine, nufine, sigr2nu_model, intbetasfine, gp)

    zetaa = -1; zetab = -1
    if gp.usezeta:
        zetaa, zetab = zeta(r0fine, nufine, \
                            Sigfine,\
                            Mrfine, betafine,\
                            sigr2nu_model, gp)

    gh.sanitize_vector(siglos2_out, len(r0), 0, 1e30, gp.debug)
    return siglos2_out, kapl4s_out, zetaa, zetab
## \fn ant_sigkaplos(r0, rhodmpar, lbaryonpar, MtoL, nupar, betapar, pop, gp)
# calculate integral for sig_los^2 * surface density,
# correspondingly for 4th order kappa,
# and virial parameters from Richardson,Fairbairn 2014
# all on gp.nfine logarithmically spaced radii,
# and store outcome on lookup table
# @param r0 array, radial bin positions, in [pc]
# @param rhodmpar overall 3D density profile in [Munit/pc^3]
# @param lbaryonpar rho* light profile from baryons
# @param MtoL mass-to-light ratio in [Msun/Lsun]
# @param nupar tracer density in [Munit/pc^3]
# @param betapar velocity anisotropy in [1]
# @param pop int for population to look at (really the scalrad for rhohalf)
# @param gp global parameters
# @return sig_los^2, correspondingly for 4th order kappa, zetaa, zetab

def zeta(r0fine, nufine, Sigfine, Mrfine, betafine, sigr2nu, gp):
    # common parameters
    N = gh.Ntot(r0fine, Sigfine, gp)
    # vr2 = sigr2nu
    # dPhidr = gu.G1__pcMsun_1km2s_2*Mrfine/r0fine**2
    # zetaa scalar
    #xint = r0fine
    #yint = nufine*(5-2*betafine)*vr2*dPhidr*r0fine**3
    #nom = gh.quadinflog(xint, yint, 0., gp.rinfty*max(gp.xepol), False)
    #yint = nufine*dPhidr*r0fine**3
    #denom = (gh.quadinflog(xint, yint, 0., gp.rinfty*max(gp.xepol), False))**2
    theta = np.arccos(r0min/r0fine)
    cth = np.cos(theta)
    sth = np.sin(theta)
    # TODO calculate nuinterp, sigr2interp, Minterp, betainterp
    yint = gu.G1__pcMsun_1km2s_2*(5-2*betainterp)*sigr2
    yint *= Minterp*rmin**2/cth**3*sth
    nom = quad(theta, yint, 0, np.pi/2)
    yint = gu.G1__pcMsun_1km2s_2**2*nuinterp*Mrinterp
    yint *= rmin**2/cth**3*sth
    denom = quad(theta, yint, 0, np.pi/2)
    zetaa = 9*N/10. * nom/denom
    # zetab scalar
    #------------------------------------------------------------
    #yint = nufine*(7-6*betafine)*vr2*dPhidr*r0fine**5
    #nom=gh.quadinflog(xint, yint, 0., gp.rinfty*max(gp.xepol), False)

    yint = gu.G1__pcMsun_1km2s_2*nuinterp*(7-6*betainterp)*sigr2interp*Mrinterp
    yint *= rmin**2/cth**3*sth
    nom = quad(theta, yint, 0, np.pi/2)

    yint = Sigfine*r0fine**3
    denom *= gh.quadinflog(xint, yint, 0., gp.rinfty*max(gp.xepol), False)

    yint = Sigmainterp*Rmin**3/cth**6*sth
    denom *= quad(theta, yint, 0, np.pi/2)
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
    yint  = gu.G1__pcMsun_1km2s_2 * Mrfine/r0fine**2  # [1/pc (km/s)^2]
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
