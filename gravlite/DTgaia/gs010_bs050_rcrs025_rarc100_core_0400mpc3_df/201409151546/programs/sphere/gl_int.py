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


def ant_intbeta(r0, r0turn, betapar, gp):
    # define function
    xint = 1.*r0
    yint = phys.beta(xint, r0turn, betapar, gp)[0]/xint

    intbet = np.zeros(len(r0))
    for i in range(len(r0)):
        intbet[i] = quad(yint, r0[0]/1.e5, r0[i])[0]
    # TODO: assumption here was that integration goes down to min(r0)/1e5

    gh.checknan(intbet, 'intbet in ant_intbeta')
    return intbet                                                # [1]
## \fn ant_intbeta(r0, r0turn, betapar, gp)
# integrate beta over r
# (integrals in front of and after sigma_r^2 integral)
# TODO: gives -9, where _old gives -1, serious speed impact
# @param r0 free variable, array, in [lunit]
# @param r0turn needed for beta calculation: turning point
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
    # add up tracer densities to get overall density profile
    # add rho* to take into account the baryonic addition 
    # (*not* Sigma from nu_i, could miss populations, have 
    # varying selection function as fct of radius, 
    # need a M/L nuisance parameter)

    # if we work on real data, add up total baryonic contribution
    if gp.investigate == 'obs':
        rhostarnu = MtoL*phys.nu(r0nu, rhostarpar, pop, gp)
        rhonu += rhostarnu
    
    #if gp.checksig:
    #    import gl_analytic as ga
    #    anrho = ga.rho_gaia(r0nu, gp)[0]
    #    rhonu = anrho
    #    #loglog(r0nu, rhonu, 'r.-')
    #    #loglog(r0nu, anrho, 'b.-')
    #    #pdb.set_trace()
    
    betanu = phys.beta(r0nu, gp.x0turn, betapar, gp)[0]
    nunu   = phys.nu(r0nu, nupar, pop, gp)
    Signu  = glp.rho_param_INT_Sig(r0nu, nupar, pop, gp)
    idnu   = ant_intbeta(r0nu, gp.x0turn, betapar, gp) # TODO check betapar are ok

    plot(r0nu, idnu, 'b.-')
    beta_star1, r_DM, gamma_star1, r_star1, r_a1, gamma_DM, rho0 = gp.files.params
    plot(gp.xfine, (np.log(gp.xfine**2+r_a1**2)-np.log(r_a1**2))/2, 'r.-')
    pdb.set_trace()

    # integrate enclosed 3D mass from 3D density
    r0tmp = np.hstack([0.,r0nu])

    rhoint = 4.*np.pi*r0nu**2*rhonu
    # add point to avoid 0.0 in Mrnu(r0nu[0])
    rhotmp = np.hstack([0.,rhoint])
    splpar_rho = splrep(r0tmp, rhotmp, k=1, s=0.) # not necessarily monotonic
    Mrnu = np.zeros(len(r0nu))              # work in refined model
    for i in range(len(r0nu)):              # get Mrnu
        Mrnu[i] = splint(0., r0nu[i], splpar_rho)
    gh.checkpositive(Mrnu, 'Mrnu')

    # (sigr2, 3D) * nu/exp(-idnu)
    xint = r0nu                           # [pc]
    yint = gp.G1 * Mrnu / r0nu**2         # [1/pc (km/s)^2]
    yint *= nunu                          # [Munit/pc^4 (km/s)^2]
    yint *= np.exp(idnu)                  # [Munit/pc^4 (km/s)^2]
    pdb.set_trace()
    gh.checkpositive(yint, 'yint sigr2')

    # use quadinflog or quadinfloglog here
    sigr2nu = np.zeros(len(r0nu))
    for i in range(len(r0nu)):
        # TODO: check quadinflog with walker profiles
        sigr2nu[i] = np.exp(-idnu[i])/nunu[i]*\
                     gh.quadinflog(xint, yint, r0nu[i], gp.rinfty*r0nu[-1], True)
        if sigr2nu[i] == np.inf:
            sigr2nu[i] = 1e-100
        # last arg: warn if no convergence found
    pdb.set_trace()
    gh.checkpositive(sigr2nu, 'sigr2nu in sigl2s')

    # project back to LOS values
    # sigl2sold = np.zeros(len(r0nu)-gp.nexp)
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
    pdb.set_trace()
    gh.checkpositive(sigl2s, 'sigl2s')

    # calculate surface density on the same r0nu as the sigl2s
    Sig = glp.rho_param_INT_Sig(r0nu, rhopar, pop, gp)
    siglos2 = sigl2s/Sig
    pdb.set_trace()

    # derefine on radii of the input vector
    splpar_sig = splrep(r0nu[:-gp.nexp], np.log(siglos2), k=3, s=0.)
    siglos2_out = np.exp(splev(r0, splpar_sig))[gp.nexp:-gp.nexp]
    gh.checkpositive(siglos2_out, 'siglos2_out')
    if not gp.usekappa:
        kapl4s_out = np.ones(len(siglos2_out))
    if gp.usekappa:
        kapl4s_out = kappa(r0nu, Mrnu, nunu, sigr2nu, idnu, gp)
        
    zetaa = np.ones(len(siglos2_out))
    zetab = np.ones(len(siglos2_out))
    if gp.usezeta:
        zetaa, zetab = zeta(r0nu[:-gp.nexp], nunu[:-gp.nexp], \
                            Signu,\
                            Mrnu[:-gp.nexp], betanu[:-gp.nexp],\
                            sigr2nu[:-gp.nexp], gp)

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


def ant_sigkaplos2surf_ipol(r0, rhopar, rhostarpar, nupar, betapar, pop, gp):
    minval = 1.e-30
    r0nu   = gh.extend_radii(r0)

    rhonu  = phys.rho(r0nu,  rhopar, 0, gp) # DM profile (first)
    # add up tracer densities to get overall density profile
    # add rho* to take into account the baryonic addition 
    # (*not* Sigma from nu_i, could miss populations, have 
    # varying selection function as fct of radius, 
    # need a M/L nuisance parameter)
    rhonu += phys.rho(r0nu, rhostarpar, 0, gp)
    
    if gp.checksig:
        import gl_analytic as ga
        rhonu = ga.rho_gaia(r0nu, gp)
        pdb.set_trace()

    nunu   = phys.nu(r0nu, nupar, pop, gp)
    betanu, dum = phys.beta(r0nu, gp.x0turn, betapar, gp)[0]

    # calculate intbeta from beta approx directly
    # TODO: check difference between ant_intbeta and ant_intbeta_old
    # usig analytical value for beta, int(beta)
    idnu   = ant_intbeta_old(r0nu, betapar, gp) # 

    # integrate enclosed 3D mass from 3D density
    r0tmp = np.hstack([0.,r0nu])
    rhoint = 4.*np.pi*r0nu**2*rhonu
    # add point to avoid 0.0 in Mrnu(r0nu[0])
    rhotmp = np.hstack([0.,rhoint])
    splpar_rho = splrep(r0tmp, rhotmp, k=3, s=0.) # not necessarily monotonic
    Mrnu = np.zeros(len(r0nu))              # work in refined model
    for i in range(len(r0nu)):              # get Mrnu
        Mrnu[i] = splint(0., r0nu[i], splpar_rho)
    gh.checkpositive(Mrnu, 'Mrnu')

    # (sigr2, 3D) * nu/exp(-idnu)
    xint = r0nu                           # [pc]
    yint = gp.G1 * Mrnu / r0nu**2         # [1/pc (km/s)^2]
    yint *= nunu                          # [Munit/pc^4 (km/s)^2]
    yint *= np.exp(idnu)                  # [Munit/pc^4 (km/s)^2]
    gh.checkpositive(yint, 'yint sigr2')

    # use quadinflog or quadinfloglog here
    sigr2nu = np.zeros(len(r0nu))
    for i in range(len(r0nu)):
        # TODO: check quadinflog with walker profiles
        sigr2nu[i] = np.exp(-idnu[i])/nunu[i]*gh.quadinflog(xint, yint, r0nu[i], np.inf, True)
        if sigr2nu[i] == np.inf:
            sigr2nu[i] = 1e-100
        # last arg: warn if no convergence found
    gh.checkpositive(sigr2nu, 'sigr2nu in sigl2s')

    # project back to LOS values
    # sigl2sold = np.zeros(len(r0nu)-gp.nexp)
    sigl2s = np.zeros(len(r0nu)-gp.nexp)
    for i in range(len(r0nu)-gp.nexp): # get sig_los^2
        xnew = np.sqrt(r0nu[i:]**2-r0nu[i]**2)                # [pc]
        ynew = 2.*(1-betanu[i]*(r0nu[i]**2)/(r0nu[i:]**2))
        ynew *= nunu[i:] * sigr2nu[i:]
        gh.checkpositive(ynew, 'ynew in sigl2s') # is hit several times..
        # check sigr2nu: has too many entries of inf!
        splpar_nu = splrep(xnew, ynew, k=1)
        # interpolation in real space for int

        sigl2s[i] = gh.quadinflog(xnew[1:], ynew[1:], xnew[0], np.inf, False)
    # for last 3 bins, we are up to a factor 2 off
    gh.checkpositive(sigl2s, 'sigl2s')
    
    # derefine on radii of the input vector
    splpar_sig = splrep(r0nu[:-gp.nexp], np.log(sigl2s), k=3, s=0.)
    sigl2s_out = np.exp(splev(r0, splpar_sig))[:-3]
    gh.checkpositive(sigl2s_out, 'sigl2s_out')
    if not gp.usekappa:
        kapl4s_out = np.ones(len(sigl2s_out))
    if gp.usekappa:
        kapl4s_out = kappa(r0nu, Mrnu, nunu, sigr2nu, idnu, gp)
        
    zetaa = np.ones(len(sigl2s_out))
    zetab = np.ones(len(sigl2s_out))
    if gp.usezeta:
        zetaa, zetab = zeta(r0nu, Mrnu, nunu, gp)

    return sigl2s_out, kapl4s_out, zetaa, zetab
## \fn ant_sigkaplos2surf_ipol(r0, rhopar, rhostarpar, nupar, betapar, pop, gp)
# calculate integral for sig_los^2 * surface density, 
# correspondingly for 4th order kappa, 
# and virial parameters from Richardson,Fairbairn 2014
# OLD VERSION, CALCULATES ALL VARIABLES ONLY ON gp.xepol
# NOT USED ANYMORE
# @param r0 array, radial bin positions, in [pc]
# @param rhopar overall 3D density profile in [Munit/pc^3]
# @param rhostarpar rho*
# @param nupar tracer density in [Munit/pc^3]
# @param betapar velocity anisotropy in [1]
# @param pop int for population to look at (really the scalrad for rhohalf)
# @param gp global parameters
# @return integral for sig_los^2 * surface density, correspondingly for 4th order kappa, zetaa, zetab


def zeta(r0nu, nunu, Signu, Mrnu, betanu, sigr2nu, gp):
    # common parameters
    N = gh.Ntot(r0nu, Signu)
    vr2 = sigr2nu
    dPhidr = gp.G1*Mrnu/r0nu**2
    
    # zetaa scalar
    xint = r0nu
    yint = nunu*(5-2*betanu)*vr2*dPhidr*r0nu**3
    nom = gh.quadinflog(xint, yint, 0., np.inf, False)
    
    yint = nunu*dPhidr*r0nu**3
    denom = (gh.quadinflog(xint, yint, 0., np.inf, False))**2
    
    zetaa = 9*N/10. * nom/denom
    
    # zetab scalar
    yint = nunu*(7-6*betanu)*vr2*dPhidr*r0nu**5
    nom = gh.quadinflog(xint, yint, 0., np.inf, False)
        
    yint = Signu*r0nu**3
    denom *= gh.quadinflog(xint, yint, 0., np.inf, False)

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
    C = max(0., gh.quadinflog(xint[-3:], yint[-3:], r0nu[-1], 1000*r0nu[-1]))
    
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
        C = max(0., gh.quadinflog(xnew[-3:], ynew[-3:], xnew[-1], 1000*xnew[-1]))
        splpar_nu = splrep(xnew,ynew) # not s=0.1, this sometimes gives negative entries after int
        kapl4s[i] = 2. * (splint(0., xnew[-1], splpar_nu) + C)
        #kapl4s[i] /= yscale
        # print('ynew = ',ynew,', kapl4s =', kapl4s[i])

    # TODO: sometimes the last value of kapl4s is nan: why?
    gh.checkpositive(kapl4s, 'kapl4s in kappa_r^4')

    # project kappa4_los as well
    # only use middle values to approximate, without errors in center and far
    kapl4s_out = np.exp(splev(r0, splrep(r0nu[4:-gp.nexp], kapl4s[4:], k=3))) # s=0.
    gh.checkpositive(kapl4s_out, 'kapl4s_out in kappa_r^4')
    return sigl2s_out, kapl4s_out


def ant_kaplos4surf(r0, beta, nu, kapr4):
    pnts = len(r0)-1
    tmp = np.zeros(pnts)

    for i in range(pnts):
        xint = r0[i:]                          # [pc]
        yint = g(r0[i:], r0[i], beta[i:])        # [1]
        yint = yint * nu[i:] * kapr4[i:] * r0[i:] # [Munit/pc^2 (km/s)^4]
        yint = yint / np.sqrt(r0[i:]**2-r0[i]**2) # [Munit/pc^3 (km/s)^4]

        # interpolate diverging bin, need 4 (sub)points to get 3rd order
        xipol = xint[0]+np.arange(4)*(xint[1]-xint[0])/4.
        dipol = beta[i]+np.arange(4)*(beta[i+1]-beta[i])/4.
        nuipol= nu[i]+np.arange(4)*(nu[i+1]-nu[i])/4.
        kapr4ipol = kapr4[i]+np.arange(4)*(kapr4[i+1]-kapr4[i])/4.

        yipol = g(xipol,r0[i],dipol)* nuipol\
                * kapr4ipol * xipol / np.sqrt(xipol**2-r0[i]**2)

        yint[0] = gh.ipol(xipol[1:4],yipol[1:4],xint[0])
        # yint[0] = gh.ipol(xint[1:3],yint[1:3],[xint[0]]) # -(xint[1]-xint[0]))
        # ^-- rather conservative, lower than above method by factor 0.6

        # then fit polynomial on second half of radii
        ml = max(len(xint)/2,2)
        x = xint[-ml:]
        y = np.log(yint[-ml:])
        polyhilo = np.polyfit(x,y, 1)

        # integrate polynomial from r0 (up to max(r_data)) to infinity
        intpoly = int_poly_inf(r0[-1],polyhilo)

        tmp[i] = 2. * (simps(yint, xint, even=gp.even) + intpoly) # [Munit/pc^2 (km/s)^4]

    # extrapolate to last 4 bins
    x = r0[-4:-1] # take values from half+some offset to last point known for tmp
    y = np.log(tmp[-3:])
    polyhilo = np.polyfit(x,y,1)

    # [Munit/pc^2 (km/s)^4]
    return np.hstack([tmp,np.exp(r0[-1:]*polyhilo[0]+polyhilo[1])])
## \fn ant_kaplos4surf(r0, beta, nu, kapr4)
# kappa_LOS^4 * surface density,
# take nu and sig_r^2, give back sig_LOS, with analytical integration to infinity
# not used anymore
# @param r0 radial bins in [pc]
# @param beta velocity dispersion in [1]
# @param nu tracer density falloff, in [Munit/pc^3]
# @param kapr4 kappa^4 in radial direction
# @return [(Munit/pc^2) (km/s)^4]


def smooth_ext(x0, y0, xext, k=1, s=0.1, log=True, slope=-3., minval=0.):
    if log:
        y0 = np.log(y0)
    tmp = splev(xext, splrep(x0, y0, k=k, s=s))

    # must impose maximum value: bounded by slope=slope and min
    if log:
        tmp = np.exp(tmp)
    if min(tmp) < minval:
        # print('extrapolation falling below ',minval,'!')
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


def int_sigr2_old(pnts, r_tot, beta_tot, M_tot, nu_tot):
    # the integrand is assumed to be zero above 2*rmax
    tmp = np.zeros(2*pnts-1)
    for i in range(2*pnts-1):
        xint  = r_tot[i:]               # [pc]
        yint  = gp.G1 * M_tot[i:]/r_tot[i:]**2 # [1/pc  Munit/msun km^2/s^2] = [1/pc (km/s)^2]
        yint *= nu_tot[i:]                     # [Munit/pc^4 (km/s)^2]
        yint *= np.exp(beta_tot[i:])          # [Munit/pc^4 (km/s)^2]
        tmp[i] = np.exp(-beta_tot[i]) * simps(yint, xint, even=gp.even)/nu[i] # [(km/s)^2]

    tmp[-1] = tmp[-2]           # [(km/s)^2]
    return tmp                  # [(km/s)^2]
## \fn int_sigr2_old(pnts, r_tot, beta_tot, M_tot, nu_tot)
# sig_r^2. no corrections at small radii
# @param pnts in [1]
# @param r_tot in [pc]
# @param beta_tot in [1]
# @param M_tot in [Munit]
# @param nu_tot in [Munit/pc^3]

def ant_intbeta_old(r0, betapar, gp):
    # extend beta by interpolating
    r0ext = [r0[0]/5., r0[0]/4., r0[0]/3., r0[0]/2.]
    r0nu = np.hstack([r0ext, r0, r0[:-1]+(r0[1:]-r0[:-1])/2.])
    r0nu.sort()
    betanu = phys.beta(r0nu, gp.x0turn, betapar, gp)[0]
    
    tmp = np.zeros(len(betanu))
    for i in range(5,len(betanu)):
        xint = r0nu[:i]                                    # [lunit]
        yint = betanu[:i]/r0nu[:i]                         # [1/lunit]
        yint[0] = correct_first_bin(xint, yint, k=3, s=0, log=False)
        tmp[i] = 2.*simps(yint, xint, even=gp.even)              # [1]
    intbet = splev(r0, splrep(r0nu[5:], tmp[5:], k=2))
    gh.checknan(intbet, 'intbet in ant_intbeta')
    return intbet                                                # [1]
## \fn ant_intbeta_old(r0, betapar, gp)
# integrate beta over r
# (integrals in front of and after sig_r^2 integral)
# @param r0 free variable, array, in [lunit]
# @param betapar integrand, array, in [1]
# @param gp


def ant_siglos2surf_old(r0, beta, nu, sigr2):
    pnts = len(r0)-1
    tmp = np.zeros(pnts)

    for i in range(pnts):
        xint = r0[i:]                          # [pc]
        yint = (1-beta[i:]*(r0[i]/r0[i:])**2) # [1]
        yint = yint * nu[i:] * sigr2[i:] * r0[i:] # [Munit/pc^2 (km/s)^2]
        yint = yint / np.sqrt(r0[i:]**2-r0[i]**2) # [Munit/pc^3 (km/s)^2]

        # interpolate diverging bin, need 4 (sub)points to get 3rd order
        xipol = xint[0]+np.arange(4)*(xint[1]-xint[0])/4.
        dipol = beta[i]+np.arange(4)*(beta[i+1]-beta[i])/4.
        nuipol= nu[i]+np.arange(4)*(nu[i+1]-nu[i])/4.
        sigr2ipol = sigr2[i]+np.arange(4)*(sigr2[i+1]-sigr2[i])/4.
        yipol = (1-dipol*(r0[i]/xipol)**2) * nuipol\
                * sigr2ipol * xipol / np.sqrt(xipol**2-r0[i]**2)
        yint[0] = gh.ipol(xipol[1:4],yipol[1:4],xint[0])

        # yint[0] = gh.ipol(xint[1:3],yint[1:3],[xint[0]]) # -(xint[1]-xint[0]))        
        # ^-- rather conservative, lower than above method by factor 0.6

        # then fit polynomial on second half of radii
        ml = max(len(xint)/2,2)
        x = xint[-ml:]
        y = np.log(yint[-ml:])
        polyhilo = np.polyfit(x,y, 1)

        # integrate polynomial from r0 (up to max(r_data)) to infinity
        intpoly = int_poly_inf(r0[-1],polyhilo)

        tmp[i] = 2. * (simps(yint, xint, even=gp.even) + intpoly) # [Munit/pc^2 (km/s)^2]

    # extrapolate to last 4 bins
    x = r0[-4:-1] # take values from half+some offset to last point known for tmp
    y = np.log(tmp[-3:])
    polyhilo = np.polyfit(x,y,1)

    # [Munit/pc^2 (km/s)^2]
    return np.hstack([tmp,np.exp(r0[-1:]*polyhilo[0]+polyhilo[1])])
## \fn ant_siglos2surf_old(r0, beta, nu, sigr2)
# take nu and sig_r^2, give back sig_LOS, with analytical integration to infinity
# not used anymore
# @param r0 radial bins in [pc]
# @param beta velocity anistropy in [1]
# @param nu 3D density falloff in [Munit/pc^3]
# @param sigr2 sig_r^2 in radial direction in [(km/s)^2]
# @return [Munit/pc^2 (km/s)^2]


def int_siglos2surf_old(pnts, r0, beta, nu, sigr2): 
    tmp = np.zeros(pnts)
    for i in range(pnts):
        xint = r0[i:]                   # [pc]
        yint = (1-beta[i:]*(r0[i]/r0[i:])**2) # [1]
        
        yint = yint * nu[i:] * sigr2[i:] * r0[i:]  # [Munit/pc^2 (km/s)^2]
        yint = yint / np.sqrt(r0[i:]**2-r0[i]**2) # [Munit/pc^3 (km/s)^2]

        yint[0] = gh.ipol(xint[1:3],yint[1:3],[xint[0]]) #-(xint[1]-xint[0]))        

        tmp[i] = 2. * simps(yint, xint, even=gp.even) # [Munit/pc^2 (km/s)^2]


    # tmp[0] = tmp[1]    # crude approximation: first bin has same sig_los2 as second bin
    # tmp[0] = gh.ipol(r_tot[1:],tmp[1:],r_tot[0])
    # tmp[0] = 2.0/surfden[0]*nu_tot[0:]*sig_tot[0:]*r_tot[0:]/r_tot[:]*rmin

    return tmp                          # [Munit/pc^2 (km/s)^2]
## \fn int_siglos2surf_old(pnts, r0, beta, nu, sigr2)
# take nu and sig_r^2, give back sig_LOS
# @param pnts integer of number of points in
# @param r0 radial bins, [pc]
# @param beta velocity anisotropy
# @param nu 3D tracer density falloff, [Munit/pc^3]
# @param sigr2 sig_r^2 in radial direction, [(km/s)^2]


