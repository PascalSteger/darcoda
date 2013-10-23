#!/usr/bin/env python3


##
# @file
# functions related to projection and deprojection of density in spherical models
# (c) 2013 Pascal Steger, ETH Zurich, psteger@phys.ethz.ch
# conventions:
# rho, r, Mr     denote 3D density, 3D radius, M(<3D radius)
# Rho, R, MR     denote 2D density, 2D radius, M(<2D radius)
# *SUM*          denotes main method = summing
# *INT*          denotes main method = integrating
# *NORM*         denotes main method = renormalization

import numpy as np
import pdb
from scipy.integrate import simps
from scipy.integrate import quad, fixed_quad, quadrature, romberg, cumtrapz
from scipy.interpolate import splrep, splev, splint

import gl_params as gp
import gl_helper as gh
from gl_int import int_poly_inf
import gl_plot as gpl

## stop execution if using wrong geometry
def enforce_sphere():
    if gp.geom == 'disc':
        print('attention: using spherical part of code for disc!')
        pdb.set_trace()
    return

## take 3D density, calculate projected surface density
# @param r0 bin radii, [pc]
# @param rho 3D density, [Msun/pc^3]
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
            # [munit/lunit^3]

        # Cubic-spline extrapolation to first bin
        tck = splrep(xint[1:],np.log(yint[1:]),k=3,s=0.)
        xnew = xint
        ynew = np.hstack([np.exp(splev(xnew[0],tck,der=0)),yint[1:]])
        # gpl.start(); pdb.set_trace()

        # power-law extension to infinity
        tck = splrep(xint[-4:],np.log(yint[-4:]),k=1,s=1.)
        invexp = lambda x: np.exp(splev(x,tck,der=0))
        dropoffint = quad(invexp,r0nu[-1],np.inf)
        tcknu = splrep(xnew,ynew,s=0) # interpolation in real space
        Rho[i] = 2. * (splint(r0nu[i], r0nu[-1], tcknu) + dropoffint[0])

        # gpl.plot(r0nu[:-4], Rho, '.', color='green')
    tcke = splrep(r0nu[:-4],Rho)
    Rhoout = splev(r0,tcke)     # [munit/lunit^2]

    # TODO: debug: raise Rho overall with factor 1.05?!
    gh.checknan(Rhoout)
    # [munit/lunit^2]

    return Rhoout

## take 3D density, calculate projected surface density
# @param r0 radii of bins, [pc]
# @param rho 3D density, [Msun/pc^3]
def rho_INT_Rho(r0, rho):
    # use splines on variable transformed integral
    # \Sigma(R) = \int_{r=R}^{R=\infty} \rho(r) d \sqrt(r^2-R^2)
    gh.checknan(rho)

    # pdb.set_trace()
    tck0 = splrep(r0,np.log(rho),k=3,s=0.01) # >= 0.1 against rising in last bin. previous: k=2, s=0.1
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
    gh.checkpositive(rhonu)
    
    Rho = np.zeros(len(r0nu)-4)
    for i in range(len(r0nu)-4):
        xnew = np.sqrt(r0nu[i:]**2-r0nu[i]**2)         # [lunit]
        ynew = 2.*rhonu[i:]
        yscale = 10.**(1.-min(np.log10(ynew)))
        ynew *= yscale

        # power-law extension to infinity
        C = gh.quadinf(xnew[-4:],ynew[-4:],xnew[-1],np.inf)
        #print('C[',i,'] = ',C)
        tcknu  = splrep(xnew,ynew,k=3) # interpolation in real space. previous: k=2, s=0.1
        Rho[i] = splint(0., xnew[-1], tcknu) + C

    Rho /= yscale
    gh.checkpositive(Rho)
    # gpl.plot(r0nu[:-4],Rho,'.')
    tcke = splrep(r0nu[:-4],Rho)
    Rhoout = splev(r0,tcke)     # [munit/lunit^2]

    gh.checkpositive(Rhoout)
    return Rhoout



## take 3D density, compute 2D density in bins, with interpolation
# @param r0 radii in [pc]
# @param rho 3D density, [Msun/pc^3]
def rho_INTIPOL_Rho(r0, rho):
    # assume len(r) = 2*pnts, or at least pnts+4
    pnts = len(r0)-1 # start with one missing bin, s.t. interpolation on sub-bin possible
    Rho = np.zeros(pnts)
    for i in range(pnts):
        xint = r0[i:]                                      # [lunit]
        yint = np.zeros(len(xint))
        for j in range(i+1,len(xint)):
            yint[j] = r0[j] * rho[j]/np.sqrt(r0[j]**2-r0[i]**2) # [munit/lunit^3]
        
        # extrapolating with 1:3 is linear only! two points (1,2) to work with!
        yint[0] = gh.ipol(xint[1:3], yint[1:3], xint[0]) # linear extrapolation
        # yint[0] = gh.ipol(xipol[1:3], yipol[1:3], xint[0], smooth=0.1) # higher smooth means lower sigmalos
        # yint[0] = gh.ipollog(xipol[1:3], yipol[1:3], xint[0]) # too jumpy
        # yint[0] = gh.ipol(xipol[1:4], yipol[1:4], xint[0]) # increase, still missing at high radii
        # yint[0] = gh.ipol(xipol[1:3],yipol[1:3],xint[0]*1.1) # overall increase, kappa off @ 0 and max
        # yint[0] = gh.ipol(xipol[1:3],yipol[1:3],xint[0])*0.5 # overall increase, but too high @ center

        # use other interpolation scheme
        # igi = InterpolatedUnivariateSpline(xint[1:5],yint[1:5])
        # yint[0] = igi(xint[0])

        # then fit polynomial on second half of radii
        ml = max(len(xint)/4,2)
        x = xint[-ml:]
        y = np.log(yint[-ml:])
        polyhilo = np.polyfit(x,y,1)

        # integrate polynomial from r0 (up to max(r_data)) to infinity
        intpoly = int_poly_inf(r0[-1],polyhilo)

        # attention: can be that intpoly is negative, if rho[i+1] > rho[i]
        # and even that integration + intpoly < 0 in the wash for Rho[i] for the last bin
        # circumvent that by requesting intpoly in [0,infty[
        if intpoly<0.: intpoly = 0.

        Rho[i]  = 2.*(simps(yint, xint, even=gp.even) + intpoly) # [munit/lunit^2]

    # better: extrapolation using polynomes
    # surfden = gh.ipollog(r_tot[1:-1],surfden[1:-1],r_tot)
    # Rho[0]  = gh.ipol(r0[1:3],Rho[1:3],[r0[0]]) # TODO: check that it's fine to disable


    # extrapolate to last 4 bins, which were neglected for sake of extrapolating yint[0]
    x = r0[-4:-1]                 # take values from half+some offset to last point known for Rho
    y = np.log(Rho[-3:])          # [log(munit/lunit^2)]
    polyhilo = np.polyfit(x,y,1)

    return np.hstack([Rho,np.exp(r0[-1:]*polyhilo[0]+polyhilo[1])])


## calculate enclosed mass from summing up (density in rings)*(ring surface)
# @param r0 radii, [pc]
# @param Rho, 2D density, [Msun/pc^2]
def Rho_SUM_MR(r0, Rho):
    MR = np.zeros(len(r0))
    MR[0] = Rho[0]*np.pi*r0[0]**2
    for i in range(1,len(r0)):
        MR[i] = MR[i-1] + Rho[i]*np.pi*(r0[i]**2 - r0[i-1]**2)
    return MR

## take 3D (mass) density, convert it to surface density, 
# and give back enclosed mass in rings
# @param r0 radii, [pc]
# @param rho 3D mass density, [Msun/pc^3]
def rho_INT_Sum_MR(r0, rho):
    enforce_sphere()
    surf_tot = rho_INT_Rho(r0, rho)                # gives [dens0, 2D]
    surfmass = Rho_SUM_MR(r0, surf_tot)            # [munit,2D]
    # [munit, 2D]
    return surfmass


## take 2D density and density error, get back 3D density and error'
# @param R0 radius, [pc]
# @param Rho 2D density, [Msun/pc^2]
# @param Rhoerr 2D density error, [Msun/pc^2]
def Rho_NORM_rho(R0, Rho, Rhoerr):
    rho =  Rho_INT_rho(R0, Rho)                      # [munit/lunit^3]
    if min(rho)<0.:
        print('*** got negative 3D density! ***')
        for i in range(len(rho)):
            if rho[i] < 0.: rho[i] = rho[i-1]
        print('corrected iteratively to last valid value')

    # normalization: calc tot mass from 2D and 3D, set equal,
    # get center 3D density rho0 right
    # and convert it into [munit/lunit^3]
    Mr = rho_SUM_Mr(R0, rho)[-1]                   # [munit]

    r0  = R0 # TODO: correction for different radii!     R0 != r0, too bad
    MR = Rho_SUM_MR(R0, Rho)[-1]                   # [munit]

    corr = MR/Mr
    rho *= corr                                      # [munit/lunit^3]
    
    # fractional error propagation
    rhoerr = rho * Rhoerr/Rho                        # [munit/lunit^3]

    # [pc], [munit/lunit^3], [munit/lunit^3]
    return r0, rho, rhoerr

## take surface density, deproject, 2D => 3D with radius r'
# @param R0 radius in [pc]
# @param Rho 2D density in [Msun/pc^2]
def Rho_INT_rho(R0, Rho):
    # TODO: deproject with variable transformed, x=sqrt(r^2-R^2)
    # TODO: any good? have 1/x in transformation, for x=0 to infty
    pnts = len(Rho)                    # [1]

    # use splines, spline derivative to go beyond first radial point
    tck0 = splrep(R0,np.log(Rho),s=0.01)
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

    R0nu = np.hstack([R0nu, R0ext]);    Rhonu = np.hstack([Rhonu, Rhoext])
    gh.checkpositive(R0nu);             gh.checkpositive(Rhonu)
    # gpl.start(); pdb.set_trace()
    tck1 = splrep(R0nu,Rhonu,k=3,        s=0.)
    # TODO: adapt s for problem at hand ==^
    dnubydR = splev(R0nu,tck1,der=1)
    # Rhofit = np.exp(splev(R0nu,tck1))   # debug for outer part

    rho = np.zeros(len(R0nu)-4)
    for i in range(len(R0nu)-4):
        xint = R0nu[i:]         # [lunit]
        yint = np.ones(len(xint))
        for j in range(i+1,len(R0nu)):
            yint[j-i] = dnubydR[j]/np.sqrt(R0nu[j]**2-R0nu[i]**2)
            # [munit/lunit^3/lunit]

        # Cubic-spline extrapolation to first bin
        tck = splrep(xint[1:],yint[1:],k=3,s=0.01)
        xnew = xint
        ynew = splev(xnew,tck,der=0)
        tcknu = splrep(xnew,ynew,s=0) # interpolation in real space
        rho[i] = -1./np.pi * splint(0.,max(R0)*2.,tcknu)
    tcke = splrep(R0nu[:-4],rho)
    rhoout = splev(R0,tcke)     # [munit/lunit^3]
    return rhoout               # [munit/lunit^3]






## take 3D density, calculate summed 3D mass
# @param r0 radii in [pc]
# @param rho 3D density in [Msun/pc^3]
def rho_SUM_Mr(r0, rho):

    # TODO: take radii from binmin, binmax for this purpose, not gp.xipol!
    
    r0 = np.hstack([0,r0])                        # [lunit]
    deltavol = 4./3.*np.pi*(r0[1:]**3 - r0[:-1]**3) # [lunit^3] 
    # TODO: check R, r, binmax, binmin
    deltaM   = rho * deltavol                  # [munit]
    Mr = np.cumsum(deltaM)                          # [munit]
    
    # old, explicit way. somehow, this is not the same as above, overestimates M
    # Mr = np.zeros(len(rdat))
    # Mr[0] = rho[0]*4./3.*np.pi*r0[1]**3
    # for i in range(1,len(M_r)):
    #     Mr[i] = Mr[i-1] + rho[i]*4./3.*np.pi*(r0[i]**3 - r0[i-1]**3)

    # assuming rho(r) = M(<r)/V(r)
    # r0 = rdat
    # vol = 4./3.*np.pi*r0**3
    # M_r = denspars * vol
        
    # if gp.cprior < 1.e29: Mr[0] = gp.cprior # [munit]
    if (Mr[0] < 0): Mr[0] = 0.              # [munit]

    # no mass can be negative (at least on these scales :o))
    if min(Mr) < 0.:
        print('M sign error! setting M_out_min = 0 for all negative values')
        for i in range(gp.nipol):
            if(Mr[i]<0.):
                Mr[i] == 0.0            # [munit]
    return Mr                           # [munit]
