#!/usr/bin/env ipython3

##
# @file
# Functions related to projection and deprojection of density in spherical models.
# Conventions:
# rho, r, Mr     denote 3D density, 3D radius, M(<3D radius)
# Rho, R, MR     denote 2D density, 2D radius, M(<2D radius)
# *SUM*          denotes main method = summing
# *INT*          denotes main method = integrating
# *NORM*         denotes main method = renormalization

# (c) 2013 Pascal Steger, ETH Zurich, psteger@phys.ethz.ch

import numpy as np
import pdb
from scipy.integrate import simps
from scipy.integrate import quad, fixed_quad, quadrature, romberg, cumtrapz
from scipy.interpolate import splrep, splev, splint
import gl_helper as gh
import gl_plot as gpl
import gl_physics as phys


def rho_INTDIRECT_Rho(r0, rho):
    # use splines, spline derivative to go beyond first radial point
    tck0 = splrep(r0,np.log(rho),k=3,s=0.1) # >= 0.1 against rising in last bin. previous: k=2, s=0.1
    # more points to the left help, but not so much
    r0ext = np.array([0.,r0[0]/6.,r0[0]/5.,r0[0]/4.,\
                      r0[0]/3.,r0[0]/2.,r0[0]/1.5])
    
    # introduce new points in between
    dR = r0[1:]-r0[:-1]
    r0nu = np.hstack([r0ext,r0,dR/4.+r0[:-1],dR/2.+r0[:-1],.75*dR+r0[:-1]])
    r0nu.sort()
    rhonu = np.exp(splev(r0nu,tck0))

    # extend to higher radii
    tckr = splrep(r0[-3:],np.log(rho[-3:]),k=1,s=0.2) # previous: k=2
    dr0 = (r0[-1]-r0[-2])/8.
    r0ext = np.hstack([r0[-1]+dr0, r0[-1]+2*dr0, r0[-1]+3*dr0, r0[-1]+4*dr0])
    rhoext = np.exp(splev(r0ext,tckr))
    r0nu = np.hstack([r0nu, r0ext]);    rhonu = np.hstack([rhonu, rhoext])

    Rho = np.zeros(len(r0nu)-4)
    for i in range(len(r0nu)-4):
        xint = r0nu[i:]         # [lunit]
        yint = np.ones(len(xint))
        for j in range(i+1,len(r0nu)):
            yint[j-i] = r0nu[j] * rhonu[j]/np.sqrt(r0nu[j]**2-r0nu[i]**2) 
            # [Munit/lunit^3]

        # Cubic-spline extrapolation to first bin
        tck = splrep(xint[1:],np.log(yint[1:]),k=3,s=0.)
        xnew = xint
        ynew = np.hstack([np.exp(splev(xnew[0],tck,der=0)),yint[1:]])
        # power-law extension to infinity
        tck = splrep(xint[-4:], np.log(yint[-4:]), k=1, s=1.)
        invexp = lambda x: np.exp(splev(x, tck, der=0))
        dropoffint = quad(invexp, r0nu[-1], np.inf)
        tcknu = splrep(xnew, ynew, s=0) # interpolation in real space
        Rho[i] = 2. * (splint(r0nu[i], r0nu[-1], tcknu) + dropoffint[0])

        # gpl.plot(r0nu[:-4], Rho, '.', color='green')
    tcke = splrep(r0nu[:-4],Rho)
    Rhoout = splev(r0,tcke)     # [Munit/lunit^2]

    gh.checkpositive(Rhoout, 'Rhoout in rho_INTDIRECT_Rho')
    # [Munit/lunit^2]
    return Rhoout
## \fn rho_INTDIRECT_Rho(r0, rho)
# take 3D density, calculate projected surface density
# @param r0 bin radii, [pc]
# @param rho 3D density, [Munit/pc^3]


def rho_INT_Rho(r0, rho):
    # use splines on variable transformed integral
    # \Sigma(R) = \int_{r=R}^{R=\infty} \rho(r) d \sqrt(r^2-R^2)
    gh.checknan(rho, 'rho_INT_Rho')

    # >= 0.1 against rising in last bin. previous: k=2, s=0.1
    tck0 = splrep(r0,np.log(rho),k=3,s=0.01)
    r0ext = np.array([0., r0[0]*0.25, r0[0]*0.50, r0[0]*0.75])
    dR = r0[1:]-r0[:-1]
    r0nu = np.hstack([r0ext,r0])
    # points in between possible, but not helping much:
    # ,dR*0.25+r0[:-1],dR*0.50+r0[:-1],dR*0.75+r0[:-1]]) 
    r0nu.sort()
    rhonu = np.exp(splev(r0nu,tck0))

    # extend to higher radii
    tckr   = splrep(r0[-3:],np.log(rho[-3:]),k=1,s=1.) # k=2 gives NaN!
    dr0    = (r0[-1]-r0[-2])/8.
    r0ext  = np.hstack([r0[-1]+dr0, r0[-1]+2*dr0, r0[-1]+3*dr0, r0[-1]+4*dr0])
    rhoext = np.exp(splev(r0ext,tckr))
    r0nu   = np.hstack([r0nu, r0ext])
    rhonu  = np.hstack([rhonu, rhoext])
    gh.checkpositive(rhonu, 'rhonu in rho_INT_Rho')
    
    Rho = np.zeros(len(r0nu)-4)
    for i in range(len(r0nu)-4):
        xnew = np.sqrt(r0nu[i:]**2-r0nu[i]**2)         # [lunit]
        ynew = 2.*rhonu[i:]
        yscale = 10.**(1.-min(np.log10(ynew)))
        ynew *= yscale

        # power-law extension to infinity
        C = gh.quadinflog(xnew[-4:],ynew[-4:],xnew[-1],np.inf)
        #print('C[',i,'] = ',C)
        tcknu  = splrep(xnew,ynew,k=3) # interpolation in real space. previous: k=2, s=0.1
        Rho[i] = splint(0., xnew[-1], tcknu) + C

    Rho /= yscale
    gh.checkpositive(Rho, 'Rho in rho_INT_Rho')
    # gpl.plot(r0nu[:-4],Rho,'.')
    tcke = splrep(r0nu[:-4], Rho)
    Rhoout = splev(r0, tcke)     # [Munit/lunit^2]

    gh.checkpositive(Rhoout, 'Rhoout in rho_INT_Rho')
    return Rhoout
## \fn rho_INT_Rho(r0, rho)
# take 3D density, calculate projected surface density
# @param r0 radii of bins, [pc]
# @param rho 3D density, [Munit/pc^3]


def rho_param_INT_Rho(r0, rhopar, pop, gp):
    # use splines on variable transformed integral
    # \Sigma(R) = \int_{r=R}^{R=\infty} \rho(r) d \sqrt(r^2-R^2)
    # gh.checknan(rhopar, 'rho_param_INT_Rho')
    xmin = r0[0]/30. # tweaked. r0[0]/1e4 gives error in quad()
    r0left = np.array([xmin, r0[0]*0.25, r0[0]*0.50, r0[0]*0.75])
    r0nu = np.hstack([r0left, r0])

    rhonu = phys.rho(r0nu, rhopar, pop, gp)
    Rho = np.zeros(len(r0nu)-gp.nexp)
    for i in range(len(r0nu)-gp.nexp):
        xnew = np.sqrt(r0nu[i:]**2-r0nu[i]**2)         # [lunit]
        ynew = 2.*rhonu[i:]

        # power-law extension to infinity
        C = gh.quadinflog(xnew[-gp.nexp:], ynew[-gp.nexp:], xnew[-1], 1e6*max(xnew)) #np.inf)
        # tcknu  = splrep(xnew, ynew, k=3)
        # interpolation in real space, not log space
        # problem: splint below could give negative values
        # reason:  for high radii (high i), have a spline that goes negative
        # workaround: multiply by const/add const to keep spline positive @ all times
        #             or set to log (but then integral is not straightforward
        # Rho[i] = splint(0., xnew[-1], tcknu) + C
        Rho[i] = gh.quadinflog(xnew[1:], ynew[1:], xmin, xnew[-1]) + C # np.inf)

    gh.checkpositive(Rho, 'Rho in rho_param_INT_Rho')
    return Rho[len(r0left):] # @r0 (r0nu without r0left, and without 3 extension bins)
## \fn rho_param_INT_Rho(r0, rhopar, pop, gp)
# take 3D density parameters, calculate projected surface density
# @param r0 radii of bins, (nepol-nexp entries) [pc]
# @param rhopar 3D density, (nepol entries) [Munit/pc^3]
# @param pop int population to take halflight radius from (0 both, 1, 2)
# @param gp global parameters


def rho_param_INT_Rho_disc(r0, rhopar, pop, gp):
    # TODO: check integration for z direction only

    # use splines on variable transformed integral
    # \Sigma(R) = \int_{r=0}^{R} \rho(r) dr
    xmin = r0[0]/30. # tweaked. r0[0]/1e4 gives error in quad()
    r0left = np.array([xmin, r0[0]*0.25, r0[0]*0.50, r0[0]*0.75])
    r0nu = np.hstack([r0left, r0])

    rhonu = phys.rho(r0nu, rhopar, pop, gp) # rho takes rho(rhalf) and n(r) parameters
    Rho = np.zeros(len(r0nu)-gp.nexp)
    for i in range(len(r0nu)-gp.nexp):
        Rho[i] = gh.quadinflog(r0nu, rhonu, xmin, r0nu[i])

    gh.checkpositive(Rho, 'Rho in rho_param_INT_Rho')
    return Rho[len(r0left):] # @r0 (r0nu without r0left, and without 3 extension bins)
## \fn rho_param_INT_Rho_disc(r0, rhopar, pop, gp)
# take 3D density parameters, calculate projected surface density
# @param r0 radii of bins, (nepol-nexp entries) [pc]
# @param rhopar 3D density, (nepol entries) [Munit/pc^3]
# @param pop int for population (0 both, 1, 2, ..)
# @param gp global parameters


def nu_param_INT_Sig_disc(r0, nupar, pop, gp):
    # TODO: check integration for z direction only
    # use splines on variable transformed integral
    # \Sigma(R) = \int_{r=0}^{R} \rho(r) dr
    r0nu = r0

    nunu = phys.nu_decrease(r0nu, nupar, gp)
    Sig = np.zeros(len(r0nu))
    for i in range(len(r0nu)):
        Sig[i] = gh.quadinflog(r0nu, nunu, r0nu[0], r0nu[i])

    gh.checkpositive(Sig, 'Sig in nu_param_INT_Sig_disc')
    return Sig
## \fn nu_param_INT_Sig_disc(r0, nupar, pop, gp)
# take 3D density parameters, calculate projected surface density
# @param r0 radii of bins, (nepol-nexp entries) [pc]
# @param nupar 3D density, (nepol entries) [Munit/pc^3]
# @param pop int for population (0 both, 1, 2, ..)
# @param gp global parameters

    
def rho_INTIPOL_Rho(r0, rho, gp):
    # use splines on variable transformed integral
    # \Sigma(R) = \int_{r=R}^{r=\infty} \rho(r) d \sqrt(r^2-R^2)
    gh.checkpositive(rho, 'rho in rho_INTIPOL_Rho')
    rmax = r0[-1]
    rext = [2*rmax, 4*rmax, 8*rmax]
    r0nu = np.hstack([r0, rext])
    nexp = 3
    rhoe = rho[-1]
    rhonu = np.hstack([rho, rhoe/1e2, rhoe/1e4, rhoe/1e8])
    Rho = np.zeros(len(r0))
    for i in range(len(r0)):
        xnew = np.sqrt(r0nu[i:]**2-r0nu[i]**2)         # [lunit]
        ynew = 2.*rhonu[i:]

        # power-law extension to infinity
        C = gh.quadinflog(xnew[-nexp:], ynew[-nexp:], xnew[-1], gp.rinfty*xnew[-1])
        Rho[i] = gh.quadinflog(xnew[1:], ynew[1:], xnew[0], xnew[-1]) + C

    gh.checkpositive(Rho, 'Rho in rho_INTIPOL_Rho')
    return Rho
## \fn rho_INTIPOL_Rho(r0, rho, gp)
# take 3D density, compute 2D density in bins, with interpolation
# \Sigma(R) = 2\int_R^\infty dr \frac{\nu(r)}{\sqrt{r^2-R^2}}
# \Sigma(R) = \int_{r=R}^{r=\infty} \rho(r) d \sqrt(r^2-R^2)
# @param r0 radii in [pc]
# @param rho 3D density, [Munit/pc^3]


def Rho_SUM_MR(r0, Rho):
    MR = np.zeros(len(r0))
    MR[0] = Rho[0]*np.pi*r0[0]**2
    for i in range(1,len(r0)):
        MR[i] = MR[i-1] + Rho[i]*np.pi*(r0[i]**2 - r0[i-1]**2)
    return MR
## \fn Rho_SUM_MR(r0, Rho)
# calculate enclosed mass from summing up (density in rings)*(ring surface)
# @param r0 radii, [pc]
# @param Rho 2D density, [Munit/pc^2]


def rho_INT_Sum_MR(r0, rho):
    surf_tot = rho_INT_Rho(r0, rho)                # gives [rho0, 2D]
    surfmass = Rho_SUM_MR(r0, surf_tot)            # [Munit,2D]
    # [Munit, 2D]
    return surfmass
## \fn rho_INT_Sum_MR(r0, rho)
# take 3D (mass) density, convert it to surface density, 
# and give back enclosed mass in rings
# @param r0 radii, [pc]
# @param rho 3D mass density, [Munit/pc^3]


def Rho_NORM_rho(R0, Rho, Rhoerr):
    try:
        rho =  Rho_INT_rho(R0, Rho)                      # [Munit/lunit^3]
    except:
        raise Exception('error in Rho_INT_rho')
        
    if min(rho)<0.:
        print('*** got negative 3D density! ***')
        for i in range(len(rho)):
            if rho[i] < 0.: rho[i] = rho[i-1]
        print('corrected iteratively to last valid value')

    # normalization: calc tot mass from 2D and 3D, set equal,
    # get center 3D density rho0 right
    # and convert it into [Munit/lunit^3]
    Mr = rho_SUM_Mr(R0, rho)                   # [Munit]

    # R0 != r0, too bad
    # TODO: correction for 3D radii
    r0  = R0

    MR = Rho_SUM_MR(R0, Rho)                   # [Munit]

    corr = MR[-1]/Mr[-1]
    rho *= corr                                      # [Munit/lunit^3]
    
    # fractional error propagation
    rhoerr = rho * Rhoerr/Rho                        # [Munit/lunit^3]

    # [pc], [Munit/lunit^3], [Munit/lunit^3]
    return r0, rho, rhoerr, Mr
## \fn Rho_NORM_rho(R0, Rho, Rhoerr)
# take 2D density and density error, get back 3D density and error'
# @param R0 radius, [pc]
# @param Rho 2D density, [Munit/pc^2]
# @param Rhoerr 2D density error, [Munit/pc^2]



def Rho_INT_rho(R0, Rho):
    # TODO: deproject with variable transformed, x = sqrt(r^2-R^2)
    pnts = len(Rho)                    # [1]

    # use splines, spline derivative to go beyond first radial point
    tck0 = splrep(R0, np.log(Rho), s=0.01)
    R0ext = np.array([R0[0]/6.,R0[0]/5.,R0[0]/4.,R0[0]/3.,R0[0]/2.,R0[0]/1.5])
    
    # introduce new points in between
    R0nu = np.hstack([R0ext,R0,(R0[1:]-R0[:-1])/2.+R0[:-1]])
    R0nu.sort()
    Rhonu = np.exp(splev(R0nu,tck0))
    
    # extend to higher radii
    tckr = splrep(R0[-3:],np.log(Rho[-3:]),k=1,s=0.3) # previous: k=2, s=0.1
    dR0 = (R0[-1]-R0[-2])/8.
    R0ext = np.hstack([R0[-1]+dR0, R0[-1]+2*dR0, R0[-1]+3*dR0, R0[-1]+4*dR0])
    Rhoext = np.exp(splev(R0ext,tckr))

    R0nu = np.hstack([R0nu, R0ext])
    Rhonu = np.hstack([Rhonu, Rhoext])
    gh.checkpositive(Rhonu, 'Rhonu in Rho_INT_rho')

    tck1 = splrep(R0nu, Rhonu, k=3, s=0.)
    dnubydR = splev(R0nu, tck1, der=1)  # Attention, numerical derivative

    rho = np.zeros(len(R0nu)-4)
    for i in range(len(R0nu)-4):
        xint = R0nu[i:]         # [lunit]
        yint = np.ones(len(xint))
        for j in range(i+1,len(R0nu)):
            yint[j-i] = dnubydR[j]/np.sqrt(R0nu[j]**2-R0nu[i]**2)
            # [Munit/lunit^3/lunit]

        # Cubic-spline extrapolation to first bin
        tck = splrep(xint[1:],yint[1:],k=3,s=0.01)
        xnew = xint
        ynew = splev(xnew,tck,der=0)
        tcknu = splrep(xnew,ynew,s=0) # interpolation in real space
        rho[i] = -1./np.pi * splint(0.,max(R0)*2.,tcknu)
    tcke = splrep(R0nu[:-4],rho)
    rhoout = splev(R0,tcke)     # [Munit/lunit^3]
    return rhoout               # [Munit/lunit^3]
## \fn Rho_INT_rho(R0, Rho):
# take surface density, deproject, 2D => 3D with radius r'
# @param R0 radius in [pc]
# @param Rho 2D density in [Munit/pc^2]


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

