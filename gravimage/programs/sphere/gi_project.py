#!/usr/bin/env ipython3

##
# @file
# Functions related to projection and deprojection of density in spherical models.
# Conventions:
# rho, r, Mr     denote 3D density, 3D radius, M(<3D radius)
# Sig, R, MR     denote 2D density, 2D radius, M(<2D radius)
# *SUM*          denotes main method = summing
# *INT*          denotes main method = integrating
# *NORM*         denotes main method = renormalization

# (c) GPL v3 2014 Pascal Steger, ETH Zurich, psteger@phys.ethz.ch

import numpy as np
import ipdb
from scipy.integrate import cumtrapz, romberg, simps, quad
from scipy.interpolate import splrep, splev, splint
import gi_helper as gh
import gi_physics as phys


def rho_INTDIRECT_Sig(r0, rho):
    # use splines, spline derivative to go beyond first radial point
    splpar_rho = splrep(r0, np.log(rho), k=3, s=0.1)
    # >= 0.1 against rising in last bin. previous: k=2, s=0.1
    # more points to the left help, but not so much
    r0ext = np.array([0.,r0[0]/6.,r0[0]/5.,r0[0]/4.,\
                      r0[0]/3.,r0[0]/2.,r0[0]/1.5])

    # introduce new points in between
    dR = r0[1:]-r0[:-1]
    r0nu = np.hstack([r0ext,r0,dR/4.+r0[:-1],dR/2.+r0[:-1],.75*dR+r0[:-1]])
    r0nu.sort()
    rhonu = np.exp(splev(r0nu,splpar_rho))

    # extend to higher radii
    splpar_lrho = splrep(r0[-3:],np.log(rho[-3:]),k=1,s=0.2) # previous: k=2
    dr0 = (r0[-1]-r0[-2])/8.
    r0ext = np.hstack([r0[-1]+dr0, r0[-1]+2*dr0, r0[-1]+3*dr0, r0[-1]+4*dr0])
    rhoext = np.exp(splev(r0ext,splpar_lrho))
    r0nu = np.hstack([r0nu, r0ext]);    rhonu = np.hstack([rhonu, rhoext])

    Sig = np.zeros(len(r0nu)-4)
    for i in range(len(r0nu)-4):
        xint = r0nu[i:]         # [lunit]
        yint = np.ones(len(xint))
        for j in range(i+1,len(r0nu)):
            yint[j-i] = r0nu[j] * rhonu[j]/np.sqrt(r0nu[j]**2-r0nu[i]**2)
            # [Munit/lunit^3]

        # Cubic-spline extrapolation to first bin
        splpar_lyint = splrep(xint[1:], np.log(yint[1:]), k=3, s=0.)
        xnew = xint
        ynew = np.hstack([np.exp(splev(xnew[0], splpar_lyint)), yint[1:]])
        # power-law extension to infinity
        splpar_llyint = splrep(xint[-4:], np.log(yint[-4:]), k=1, s=1.)
        invexp = lambda x: np.exp(splev(x, splpar_llyint))
        dropoffint = quad(invexp, r0nu[-1], np.inf)
        splpar_nu = splrep(xnew, ynew, s=0) # interpolation in real space
        Sig[i] = 2. * (splint(r0nu[i], r0nu[-1], splpar_nu) + dropoffint[0])

    Sigout = splev(r0, splrep(r0nu[:-4],Sig))     # [Munit/lunit^2]

    gh.checkpositive(Sigout, 'Sigout in rho_INTDIRECT_Sig')
    # [Munit/lunit^2]
    return Sigout
## \fn rho_INTDIRECT_Sig(r0, rho)
# take 3D density, calculate projected surface density
# @param r0 bin radii, [pc]
# @param rho 3D density, [Munit/pc^3]


def rho_INT_Sig(r0, rho, gp):
    # use splines on variable transformed integral
    # \Sigma(R) = \int_{r=R}^{R=\infty} \rho(r) d \sqrt(r^2-R^2)
    gh.checknan(rho, 'rho_INT_Sig')

    # >= 0.1 against rising in last bin. previous: k=2, s=0.1
    splpar_rho = splrep(r0, np.log(rho), k=3, s=0.01)
    r0ext = np.array([0., r0[0]*0.25, r0[0]*0.50, r0[0]*0.75])
    dR = r0[1:]-r0[:-1]
    r0nu = np.hstack([r0ext,r0])
    # points in between possible, but not helping much:
    # ,dR*0.25+r0[:-1],dR*0.50+r0[:-1],dR*0.75+r0[:-1]])
    r0nu.sort()
    rhonu = np.exp(splev(r0nu, splpar_rho))

    # extend to higher radii
    splpar_lrhor   = splrep(r0[-3:],np.log(rho[-3:]),k=1,s=1.) # k=2 gives NaN!
    dr0    = (r0[-1]-r0[-2])/8.
    r0ext  = np.hstack([r0[-1]+dr0, r0[-1]+2*dr0, r0[-1]+3*dr0, r0[-1]+4*dr0])
    rhoext = np.exp(splev(r0ext, splpar_lrhor))
    r0nu   = np.hstack([r0nu, r0ext])
    rhonu  = np.hstack([rhonu, rhoext])
    gh.checkpositive(rhonu, 'rhonu in rho_INT_Sig')

    Sig = np.zeros(len(r0nu)-4)
    for i in range(len(r0nu)-4):
        xnew = np.sqrt(r0nu[i:]**2-r0nu[i]**2)         # [lunit]
        ynew = 2.*rhonu[i:]
        yscale = 10.**(1.-min(np.log10(ynew)))
        ynew *= yscale

        # power-law extension to infinity
        C = gh.quadinflog(xnew[-4:],ynew[-4:],xnew[-1], gp.rinfty*max(gp.xepol))

        splpar_nu  = splrep(xnew,ynew,k=3) # interpolation in real space. previous: k=2, s=0.1
        Sig[i] = splint(0., xnew[-1], splpar_nu) + C

    Sig /= yscale
    gh.checkpositive(Sig, 'Sig in rho_INT_Sig')
    Sigout = splev(r0, splrep(r0nu[:-4], Sig))     # [Munit/lunit^2]

    gh.checkpositive(Sigout, 'Sigout in rho_INT_Sig')
    return Sigout
## \fn rho_INT_Sig(r0, rho, gp)
# take 3D density, calculate projected surface density
# @param r0 radii of bins, [pc]
# @param rho 3D density, [Munit/pc^3]
# @param gp global parameters


def rho_param_INT_Sig(r0, rhodmpar, pop, gp):
    # use splines on variable transformed integral
    # \Sigma(R) = \int_{r=R}^{R=\infty} \rho(r) d \sqrt(r^2-R^2)
    #gh.sanitize_vector(rhodmpar, gp.nrho, -gp.nrtol, \
    #                   max(gp.maxrhoslope, 10**(np.log10(gp.rhohalf)+gp.log10rhospread)), gp.debug)
    xmin = gp.xfine[0]/15. # needed, if not: loose on first 4 bins
    r0nu = gp.xfine

    rhonu = phys.rho(r0nu, rhodmpar, pop, gp)
    Sig = np.zeros(len(r0nu)-gp.nexp)
    for i in range(len(r0nu)-gp.nexp):
        xnew = np.sqrt(r0nu[i:]**2-r0nu[i]**2)         # [lunit]
        ynew = 2.*rhonu[i:]

        # power-law extension to infinity
        C = gh.quadinflog(xnew[-gp.nexp:], ynew[-gp.nexp:], xnew[-1], gp.rinfty*xnew[-1])
        Sig[i] = gh.quadinflog(xnew[1:], ynew[1:], xmin, xnew[-1]) + C # np.inf)

    gh.checkpositive(Sig, 'Sig in rho_param_INT_Sig')

    # interpolation onto r0
    tck1 = splrep(np.log(gp.xfine[:-gp.nexp]), np.log(Sig))
    Sigout = np.exp(splev(np.log(r0), tck1))
    return Sigout
## \fn rho_param_INT_Sig(r0, rhodmpar, pop, gp)
# take 3D density parameters, calculate projected surface density
# @param r0 radii of bins, [pc]
# @param rhodmpar 3D density, (nrho entries) [Munit/pc^3]
# @param pop int population to take halflight radius from (0 both, 1, 2)
# @param gp global parameters


def rho_param_INT_Sig_theta(Rproj, rhodmpar, pop, gp):
    # use splines on variable transformed integral
    # \Sigma(R) = \int_{r=R}^{R=\infty} \rho(r) d \sqrt(r^2-R^2)
    gh.sanitize_vector(rhodmpar, gp.nrho, -gp.nrtol, \
                       max(gp.maxrhoslope, 10**(np.log10(gp.rhohalf)+gp.log10rhospread)), gp.debug)
    bit = 1.e-6
    theta = np.linspace(0, np.pi/2-bit, gp.nfine)
    cth = np.cos(theta)
    cth2 = cth*cth
    rhonu = phys.rho(Rproj, rhodmpar, pop, gp)
    Sig = np.zeros(len(Rproj))
    for i in range(len(Rproj)):
        rq = Rproj[i]/cth
        rhoq = np.interp(rq, Rproj, rhonu, left=0, right=0)
        #right=rhonu[-1]/1e10) # best for hern
        #rhoq = phys.rho(rq, rhodmpar, pop, gp)
        Sig[i] = 2.*Rproj[i]*simps(rhoq/cth2, theta)
    gh.checkpositive(Sig, 'Sig in rho_param_INT_Sig')

    # interpolation onto r0
    #tck1 = splrep(np.log(gp.xfine), np.log(Sig))
    #Sigout = np.exp(splev(np.log(r0), tck1))
    return Sig
## \fn rho_param_INT_Sig_theta(Rproj, rhodmpar, pop, gp)
# take 3D density parameters, calculate projected surface density with theta substitution
# @param Rproj radii of 2D bins, [pc]
# @param rhodmpar 3D density parameters, (nrho entries) [Munit/pc^3]
# @param pop int population to take halflight radius from (0 both, 1, 2)
# @param gp global parameters


def Sig_SUM_MR(r0, Sig):
    MR = np.zeros(len(r0))
    MR[0] = Sig[0]*np.pi*r0[0]**2
    for i in range(1,len(r0)):
        MR[i] = MR[i-1] + Sig[i]*np.pi*(r0[i]**2 - r0[i-1]**2)
    return MR
## \fn Sig_SUM_MR(r0, Sig)
# calculate enclosed mass from summing up (density in rings)*(ring surface)
# @param r0 radii, [pc]
# @param Sig 2D density, [Munit/pc^2]


def rho_INT_Sum_MR(r0, rho, gp):
    surf_tot = rho_INT_Sig(r0, rho, gp)                # gives [rho0, 2D]
    surfmass = Sig_SUM_MR(r0, surf_tot)            # [Munit,2D]
    # [Munit, 2D]
    return surfmass
## \fn rho_INT_Sum_MR(r0, rho, gp)
# take 3D (mass) density, convert it to surface density,
# and give back enclosed mass in rings
# @param r0 radii, [pc]
# @param rho 3D mass density, [Munit/pc^3]
# @param gp global parameters


def Sig_NORM_rho(R0, Sig, Sigerr, gp):
    rho =  Sig_INT_rho(R0, Sig, gp)         # [Munit/lunit^3]

    if min(rho)<0.:
        gh.LOG(1, '*** Sig_NORM_rho: got bin with negative 3D density! ***')
        for i in range(len(rho)):
            if rho[i] < 0.: rho[i] = rho[i-1]
            gh.LOG(2, 'corrected iteratively to last valid value')

    # normalization: calc tot mass from 2D and 3D, set equal,
    # get center 3D density rho0 right
    # and convert it into [Munit/lunit^3]
    Mr = rho_SUM_Mr(R0, rho)                   # [Munit]

    # R0 != r0, too bad
    # TODO: correction for 3D radii
    r0  = R0

    MR = Sig_SUM_MR(R0, Sig)                   # [Munit]

    corr = MR[-1]/Mr[-1]
    gh.LOG(2, ' * Sig_NORM_rho:  no correction by ', corr)
    # rho *= corr                                      # [Munit/lunit^3]

    # fractional error propagation
    # TODO right error determination needs to be included in Sig_INT_rho()
    rhoerr = rho * Sigerr/Sig                        # [Munit/lunit^3]

    # [pc], [Munit/lunit^3], [Munit/lunit^3]
    return r0, rho, rhoerr, Mr
## \fn Sig_NORM_rho(R0, Sig, Sigerr, gp)
# take 2D density and density error, get back 3D density and error'
# @param R0 radius, [pc]
# @param Sig 2D density, [Munit/pc^2]
# @param Sigerr 2D density error, [Munit/pc^2]
# @param gp global parameters


def Sig_INT_rho(R0, Sig, gp):
    splpar_Sig = splrep(R0, np.log(Sig)) # get spline in log space to circumvent flattening
    J = np.zeros(len(Sig)-gp.nexp)
    for i in range(len(Sig)-gp.nexp):
        xint = np.sqrt(R0[i:]**2-R0[i]**2)
        yint = np.exp(splev(np.sqrt(xint**2+R0[i]**2), splpar_Sig))
        J[i] = gh.quadinflog(xint, yint, 0., gp.rinfty*max(gp.xepol))

        #yf = lambda r: np.exp(splev(np.sqrt(r**2+R0[i]**2), splpar_Sig))
        #J[i] = quad(yf, 0, np.inf)[0]

        # OLD: direct integration with integrand as fct of x
        #xint = R0[i:]
        #yint = Sig[i:]*R0[i:]/np.sqrt(R0[i:]**2-R0[i]**2)
        #J[i] = gh.quadinflog(xint[1:], yint[1:], R0[i], np.inf)
    gh.checkpositive(J)
    sm0 = 0.02
    splpar_J = splrep(R0[:-gp.nexp], np.log(J), s=sm0) # smoothing determined for Hernquist profile. can go below to 0.01 for very smooth profiles
    rho = -1./np.pi/R0[:-gp.nexp]*J*splev(R0[:-gp.nexp], splpar_J, der=1)
    sm = sm0*1.
    while(min(rho)<0):
        print('min(rho)<0, rho = ', rho)
        sm *= 2
        splpar_J = splrep(R0[:-gp.nexp], np.log(J), s=sm)
        rho = -1./np.pi/R0[:-gp.nexp]*J*splev(R0[:-gp.nexp], splpar_J, der=1)
        #if(sm>1):
        #    raise Exception('Very irregular profile')

    gh.checkpositive(rho)

    rhoright = gh.expolpowerlaw(R0[:-gp.nexp], rho, R0[-3:], -3.001)
    rho = np.hstack([rho, rhoright])
    return rho
## \fn Sig_INT_rho(R0, Sig, gp)
# Abel transformation to get 3D mass density from 2D surface density
# see appendix in Mamon Boue 2009
# @param R0 radii in [pc]
# @param Sig surface density
# @param gp global parameters

def rho_SUM_Mr(r0max, rho):
    r0max = np.hstack([0., r0max])
    deltavol = 4./3.*np.pi*(r0max[1:]**3 - r0max[:-1]**3)    # [lunit^3]
    deltaM   = rho * deltavol                       # [Munit]
    Mr = np.cumsum(deltaM)                          # [Munit]
    return Mr                                       # [Munit]
## \fn rho_SUM_Mr(r0max, rho)
# take 3D density, calculate summed 3D mass
# @param r0max radii in [pc]
# @param rho 3D density in [Munit/pc^3]
