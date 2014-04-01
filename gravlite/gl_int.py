#!/usr/bin/env python3

##
# @file
# all integrals from gl_physics

# (c) 2013 Pascal S.P. Steger

import numpy as np
import pdb
import scipy
from scipy.integrate import simps,trapz,quad
from scipy.interpolate import splrep, splev, splint
import gl_chi as gc
import gl_helper as gh
import gl_plot as gpl

# from gl_analytic import *
import gl_analytic as ga
import gl_physics as phys

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
    tck = splrep(xint[1:], yint[1:], k=k, s=s)
    tmp = splev(xint[0], tck, der=0)
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


def ant_intbeta(r0, betaparam, gp):
    # define function
    yint = lambda x: phys.beta(x, betaparam, gp)/x

    intbet = np.zeros(len(r0))
    for i in range(len(r0)):
        intbet[i] = quad(yint, r0[0]/1.e5, r0[i])[0]

    gh.checknan(intbet, 'intbet in ant_intbeta')
    return intbet                                                # [1]
## \fn ant_intbeta(r0, betaparam, gp)
# integrate beta over r
# (integrals in front of and after sigma_r^2 integral)
# TODO: gives -9, where _old gives -1, serious speed impact
# @param r0 free variable, array, in [lunit]
# @param betaparam integrand, array, in [1]
# @param gp


def ant_intbeta_old(r0, betaparam, gp):
    # extend beta by interpolating
    r0ext = [r0[0]/5., r0[0]/4., r0[0]/3., r0[0]/2.]
    r0nu = np.hstack([r0ext, r0, r0[:-1]+(r0[1:]-r0[:-1])/2.])
    r0nu.sort()
    betanu = phys.beta(r0nu, betaparam, gp)
    
    tmp = np.zeros(len(betanu))
    for i in range(5,len(betanu)):
        xint = r0nu[:i]                                    # [lunit]
        yint = betanu[:i]/r0nu[:i]                         # [1/lunit]
        yint[0] = correct_first_bin(xint, yint, k=3, s=0, log=False)
        tmp[i] = 2.*simps(yint, xint, even=gp.even)              # [1]
    tckout = splrep(r0nu[5:], tmp[5:], k=2)
    intbet = splev(r0, tckout)
    gh.checknan(intbet, 'intbet in ant_intbeta')
    return intbet                                                # [1]
## \fn ant_intbeta_old(r0, betaparam, gp)
# integrate beta over r
# (integrals in front of and after sigma_r^2 integral)
# @param r0 free variable, array, in [lunit]
# @param betaparam integrand, array, in [1]
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






def introduce_points_in_between(r0):
    r0ext = np.array([r0[0]*0.25, r0[0]*0.50, r0[0]*0.75])
    dR = r0[1:] - r0[:-1]
    r0nu = np.hstack([r0ext, r0, dR/2.+r0[:-1]])
    r0nu.sort()
    return r0nu
## \fn introduce_points_in_between(r0)
# insert equally spaced points at 1/4, 1/2, and 3/4 of a bin span
# @param r0 original bin centers


def ant_sigkaplos2surf(r0, beta_param, rho_param, nu_param, gp):
    minval = 1.e-30
    r0nu   = introduce_points_in_between(r0)

    rhonu  = phys.rho(r0nu,  rho_param, gp)
    nunu   = phys.rho(r0nu,  nu_param, gp)
    betanu = phys.beta(r0nu, beta_param, gp)

    # calculate intbeta from beta approx directly
    idnu   = ant_intbeta_old(r0nu, beta_param, gp)

    # integrate enclosed 3D mass from 3D density
    r0tmp = np.hstack([0.,r0nu])
    rhoint = 4.*np.pi*r0nu**2*rhonu
    # add point to avoid 0.0 in Mrnu(r0nu[0])
    rhotmp = np.hstack([0.,rhoint])
    tck1 = splrep(r0tmp, rhotmp, k=3, s=0.) # not necessarily monotonic
    Mrnu = np.zeros(len(r0nu))              # work in refined model
    for i in range(len(r0nu)):              # get Mrnu
        Mrnu[i] = splint(0., r0nu[i], tck1)
    gh.checkpositive(Mrnu, 'Mrnu')

    # (sigr2, 3D) * nu/exp(-idnu)
    xint = r0nu                           # [pc]
    yint = gp.G1 * Mrnu / r0nu**2         # [1/pc (km/s)^2]
    yint *= nunu                          # [munit/pc^4 (km/s)^2]
    yint *= np.exp(idnu)                  # [munit/pc^4 (km/s)^2]
    gh.checkpositive(yint, 'yint sigr2')

    # use quadinflog or quadinfloglog here
    sigr2nu = np.zeros(len(r0nu))
    for i in range(len(r0nu)):
        sigr2nu[i] = np.exp(-idnu[i])/nunu[i]*gh.quadinflog(xint, yint, r0nu[i], np.inf, False)

    # project back to LOS values
    # sigl2sold = np.zeros(len(r0nu)-gp.nexp)
    sigl2s = np.zeros(len(r0nu)-gp.nexp)
    dropoffintold = 1.e30
    for i in range(len(r0nu)-gp.nexp): # get sig_los^2
        xnew = np.sqrt(r0nu[i:]**2-r0nu[i]**2)                # [pc]
        ynew = 2.*(1-betanu[i]*(r0nu[i]**2)/(r0nu[i:]**2))
        ynew *= nunu[i:] * sigr2nu[i:]
        gh.checkpositive(ynew, 'ynew in sigl2s') # is hit several times..
        # check sigr2nu: has too many entries of inf!
        tcknu = splrep(xnew, ynew, k=1)
        # interpolation in real space for int

        sigl2s[i] = gh.quadinflog(xnew[1:], ynew[1:], xnew[0], np.inf, False)
    # for last 3 bins, we are up to a factor 2 off

    gh.checkpositive(sigl2s, 'sigl2s')
    
    # derefine on radii of the input vector
    tck = splrep(r0nu[:-gp.nexp], np.log(sigl2s), k=3, s=0.)
    sigl2s_out = np.exp(splev(r0, tck))
    gh.checkpositive(sigl2s_out, 'sigl2s_out')
    if not gp.usekappa:
        return sigl2s_out, np.ones(len(sigl2s_out))

    # for the following: enabled calculation of kappa
    # TODO: include another set of anisotropy parameters beta_'
    # TODO: replace kappa^4 and beta with version using virial parameters

    # kappa_r^4
    kapr4nu = np.ones(len(r0nu)-gp.nexp)
    xint  = r0nu                  # [pc]
    yint  = gp.G1 * Mrnu/r0nu**2  # [1/pc (km/s)^2]
    yint *= nunu                  # [munit/pc^4 (km/s)^2]
    yint *= sigr2nu               # [munit/pc^4 (km/s)^4
    yint *= np.exp(idnu)          # [munit/pc^4 (km/s)^4]
    gh.checkpositive(yint, 'yint in kappa_r^4')
    yscale = 10.**(1.-min(np.log10(yint[1:])))
    yint *= yscale
    # power-law extrapolation to infinity
    C = max(0., gh.quadinflog(xint[-3:], yint[-3:], r0nu[-1], 1000*r0nu[-1]))
    
    tcknu = splrep(xint, yint, k=3) # interpolation in real space
    for i in range(len(r0nu)-gp.nexp):
        # integrate from minimal radius to infinity
        kapr4nu[i] = 3.*(np.exp(-idnu[i])/nunu[i]) * \
            (splint(r0nu[i], r0nu[-1], tcknu) + C) # [(km/s)^4]

    kapr4nu /= yscale
    gh.checkpositive(kapr4nu, 'kapr4nu in kappa_r^4')

    tcke = splrep(r0nu[:-gp.nexp], np.log(kapr4nu), k=3)
    kapr4ext = np.exp(splev(r0ext, tcke))
    kapr4nu = np.hstack([kapr4nu, kapr4ext])
    gh.checkpositive(kapr4nu, 'kapr4nu in extended kappa_r^4')
    
    tckbet = splrep(r0nu, betanu)
    dbetanudr = splev(r0nu, tckbet, der=1)
    gh.checknan(dbetanudr, 'dbetanudr in kappa_r^4')
    
    # kappa^4_los*surfdensity
    kapl4s = np.zeros(len(r0nu)-gp.nexp)
    #    gpl.start(); gpl.yscale('linear')
    for i in range(len(r0nu)-gp.nexp):
        xnew = np.sqrt(r0nu[i:]**2-r0nu[i]**2)      # [pc]
        ynew = g(r0nu[i:], r0nu[i], betanu[i:], dbetanudr[i:]) # [1]
        ynew *= nunu[i:] * kapr4nu[i:]
        C = max(0., gh.quadinflog(xnew[-3:], ynew[-3:], xnew[-1], 1000*xnew[-1]))
        tcknu = splrep(xnew,ynew) # not s=0.1, this sometimes gives negative entries after int
        kapl4s[i] = 2. * (splint(0., xnew[-1], tcknu) + C)
        #kapl4s[i] /= yscale
        # print('ynew = ',ynew,', kapl4s =', kapl4s[i])

    # TODO: sometimes the last value of kapl4s is nan: why?
    gh.checkpositive(kapl4s, 'kapl4s in kappa_r^4')

    # project kappa4_los as well
    # only use middle values to approximate, without errors in center and far
    tck = splrep(r0nu[4:-gp.nexp], kapl4s[4:], k=3) # s=0.
    kapl4s_out = np.exp(splev(r0, tck))
    gh.checkpositive(kapl4s_out, 'kapl4s_out in kappa_r^4')
    return sigl2s_out, kapl4s_out
## \fn ant_sigkaplos2surf(r0, beta_param, rho_param, nu_param)
# calculate integral for sigma_los^2 * surface density, correspondingly for 4th order kappa
# @param r0 array, radial bin positions, in [pc]
# @param beta velocity anisotropy in [1]
# @param intbeta beta/s integrated over ds from 0 to r
# @param rho overall 3D density profile in [Msun/pc^3]
# @param nu tracer density in [Msun/pc^3]
# @return integral for sigma_los^2 * surface density, correspondingly for 4th order kappa


def int_sigr2_old(pnts, r_tot, beta_tot, M_tot, nu_tot):
    # the integrand is assumed to be zero above 2*rmax
    tmp = np.zeros(2*pnts-1)
    for i in range(2*pnts-1):
        xint  = r_tot[i:]               # [pc]
        yint  = gp.G1 * M_tot[i:]/r_tot[i:]**2 # [1/pc  munit/msun km^2/s^2] = [1/pc (km/s)^2]
        yint *= nu_tot[i:]                     # [munit/pc^4 (km/s)^2]
        yint *= np.exp(beta_tot[i:])          # [munit/pc^4 (km/s)^2]
        tmp[i] = np.exp(-beta_tot[i]) * simps(yint, xint, even=gp.even)/nu[i] # [(km/s)^2]

    # could calculate that from last integral already, trying to extrapolate anyhow
    # if (sigmaprior < -1): sigmaprior=-1.
    # nusig2_outer=nusig2[-1]+(np.arange(1,pnts))*sigmaprior*nusig2[pnts-1]/(xmax-xmin)*dx
    # nusigma2_tot=np.hstack((nusig2,nusig2_outer))
    # sigr2_tot = gh.ipollog(r_tot[:-1],sigr2[:-1],r_tot) # beware of interpolations, dude!
    tmp[-1] = tmp[-2]           # [(km/s)^2]
    return tmp                  # [(km/s)^2]
## \fn int_sigr2_old(pnts, r_tot, beta_tot, M_tot, nu_tot)
# sigma_r^2. no corrections at small radii
# @param pnts in [1]
# @param r_tot in [pc]
# @param beta_tot in [1]
# @param M_tot in [Msun]
# @param nu_tot in [Msun/pc^3]

def ant_kaplos4surf(r0, beta, nu, kapr4):
    pnts = len(r0)-1
    tmp = np.zeros(pnts)

    for i in range(pnts):
        xint = r0[i:]                          # [pc]
        yint = g(r0[i:], r0[i], beta[i:])        # [1]
        yint = yint * nu[i:] * kapr4[i:] * r0[i:] # [munit/pc^2 (km/s)^4]
        yint = yint / np.sqrt(r0[i:]**2-r0[i]**2) # [munit/pc^3 (km/s)^4]

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

        tmp[i] = 2. * (simps(yint, xint, even=gp.even) + intpoly) # [munit/pc^2 (km/s)^4]

    # extrapolate to last 4 bins
    x = r0[-4:-1] # take values from half+some offset to last point known for tmp
    y = np.log(tmp[-3:])
    polyhilo = np.polyfit(x,y,1)

    # [munit/pc^2 (km/s)^4]
    return np.hstack([tmp,np.exp(r0[-1:]*polyhilo[0]+polyhilo[1])])
## \fn ant_kaplos4surf(r0, beta, nu, kapr4)
# kappa_LOS^4 * surface density,
# take nu and sig_r^2, give back sigma_LOS, with analytical integration to infinity
# not used anymore
# @param r0 radial bins in [pc]
# @param beta velocity dispersion in [1]
# @param nu tracer density falloff, in [Msun/pc^3]
# @param kapr4 kappa^4 in radial direction
# @return [(munit/pc^2) (km/s)^4]


def ant_siglos2surf_old(r0, beta, nu, sigr2):
    pnts = len(r0)-1
    tmp = np.zeros(pnts)

    for i in range(pnts):
        xint = r0[i:]                          # [pc]
        yint = (1-beta[i:]*(r0[i]/r0[i:])**2) # [1]
        yint = yint * nu[i:] * sigr2[i:] * r0[i:] # [munit/pc^2 (km/s)^2]
        yint = yint / np.sqrt(r0[i:]**2-r0[i]**2) # [munit/pc^3 (km/s)^2]

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

        tmp[i] = 2. * (simps(yint, xint, even=gp.even) + intpoly) # [munit/pc^2 (km/s)^2]

    # extrapolate to last 4 bins
    x = r0[-4:-1] # take values from half+some offset to last point known for tmp
    y = np.log(tmp[-3:])
    polyhilo = np.polyfit(x,y,1)

    # [munit/pc^2 (km/s)^2]
    return np.hstack([tmp,np.exp(r0[-1:]*polyhilo[0]+polyhilo[1])])
## \fn ant_siglos2surf_old(r0, beta, nu, sigr2)
# take nu and sig_r^2, give back sigma_LOS, with analytical integration to infinity
# not used anymore
# @param r0 radial bins in [pc]
# @param beta velocity anistropy in [1]
# @param nu 3D density falloff in [Msun/pc^3]
# @param sigr2 sigma^2 in radial direction in [(km/s)^2]
# @return [munit/pc^2 (km/s)^2]


def int_siglos2surf_old(pnts, r0, beta, nu, sigr2): 
    tmp = np.zeros(pnts)
    for i in range(pnts):
        xint = r0[i:]                   # [pc]
        yint = (1-beta[i:]*(r0[i]/r0[i:])**2) # [1]
        
        yint = yint * nu[i:] * sigr2[i:] * r0[i:]  # [munit/pc^2 (km/s)^2]
        yint = yint / np.sqrt(r0[i:]**2-r0[i]**2) # [munit/pc^3 (km/s)^2]

        yint[0] = gh.ipol(xint[1:3],yint[1:3],[xint[0]]) #-(xint[1]-xint[0]))        

        tmp[i] = 2. * simps(yint, xint, even=gp.even) # [munit/pc^2 (km/s)^2]


    # tmp[0] = tmp[1]    # crude approximation: first bin has same sigma_los2 as second bin
    # tmp[0] = gh.ipol(r_tot[1:],tmp[1:],r_tot[0])
    # tmp[0] = 2.0/surfden[0]*nu_tot[0:]*sigma_tot[0:]*r_tot[0:]/r_tot[:]*rmin

    return tmp                          # [munit/pc^2 (km/s)^2]
## \fn int_siglos2surf_old(pnts, r0, beta, nu, sigr2)
# take nu and sig_r^2, give back sigma_LOS
# @param pnts integer of number of points in
# @param r0 radial bins, [pc]
# @param beta velocity anisotropy
# @param nu 3D tracer density falloff, [Msun/pc^3]
# @param sigr2 sigma^2 in radial direction, [(km/s)^2]


def smooth_ext(x0, y0, xext, k=1, s=0.1, log=True, slope=-3., minval=0.):
    if log:
        y0 = np.log(y0)
    tck = splrep(x0, y0, k=k, s=s)
    tmp = splev(xext,tck,der=0)

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
    tck = splrep(x0, y0, k=k, s=s)
    tmp = splev(xnew,tck,der=0)
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
