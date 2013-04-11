#!/usr/bin/python2.7
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
'''all functions working on spherical coordinates'''
import pdb
import numpy as np
import gl_params as gp
import gl_funs as gfun
import gl_plot as gpl
from gl_int import *
from gl_analytic import *




def rhodm_hernquist(r,rho0,r_DM,alpha_DM,beta_DM,gamma_DM):
    'return hernquist DM density from generalized Hernquist profile'
    return rho0*(r/r_DM)**(-gamma_DM)*(1+(r/r_DM)**alpha_DM)**((gamma_DM-beta_DM)/alpha_DM)




def deproject(x, nu2d, err2d): #[pc], [munit/lunit^2], [munit/lunit^2]
    'take 2D density and density error, get back 3D density and error'
    nu3d =  int_2D3D(x, nu2d) #[munit/lunit^3]
    if min(nu3d)<0.:
        print 'got negative 3D density!'
        exit(1)

    # normalization: calc tot mass from 2D and 3D, set equal, get center 3D density rho0 right
    # and convert it into [munit/lunit^3]
    totmass_3D = Mr3D(x, nu3d)[-1] # [munit]
    totmass_2D = Mr2D(x, nu2d)[-1] # [munit]
    corr = totmass_2D/totmass_3D
    nu3d *= corr
    
    # fractional error propagation
    err3d = nu3d * err2d/nu2d #[munit/lunit^3]

    return x, nu3d, err3d #[pc], [munit/lunit^3], [munit/lunit^3]




def densdefault(denspars):
    'return density from denspars, assuming default gp.xipol as radius'
    return dens(gp.xipol,denspars)        # [TODO]




# def dens(xipol, denspars):
#     tmp = np.zeros(len(xipol))
#     # tmp += denspars[0] # (*xipol**0) # TODO: check unnecessary, included below
#     for i in range(0,len(denspars)):
#         tmp += denspars[i]*(xipol)**(-i)
#     # gpl.plot(gp.ipol.densr,gp.ipol.densdat); gpl.plot(gp.xipol,10.**tmp); gpl.yscale('log')
#     return 10.0**tmp



def dens(xipol, denspars):                # [pc], [1] or [msun/pc^3]
    'take denspars for polynomial coefficients, calculate polynomial, give back density'
    if gp.checksigma:
        return gp.ipol.densdat
    if not gp.poly:
        return denspars

    scale = gp.scaledens*max(xipol)     # [pc]

    tmp = np.zeros(len(xipol))
    for i in range(0,len(denspars)):
        tmp += denspars[i]*((scale-xipol)/scale)**i # [log10(msun/pc^3)]
    # gpl.plot(gp.ipol.densr,gp.ipol.densdat); gpl.plot(gp.xipol,10.**tmp); gpl.yscale('log')
    return 10.0**tmp                    # [msun/pc^3]






def calculate_dens(r, M):                           # [lunit], [munit]
    'take enclosed mass, radii, and compute 3D density in shells'
    r0 = np.hstack([0,r])                         # [lunit]
    M0 = np.hstack([0,M])

    deltaM   = M0[1:]-M0[:-1]                     # [munit]
    if min(deltaM)<0.:
        print 'unphysical negative mass increment encountered'
        exit(1)

    deltavol = 4./3.*np.pi*(r0[1:]**3-r0[:-1]**3) # [lunit^3]
    dens     = deltaM / deltavol                  # [munit/lunit^3]

    return dens                                   # [munit/lunit^3]






def calculate_surfdens(r, M):                       # [lunit], [munit]
    'take mass(<r) in bins, calc mass in annuli, get surface density'
    r0 = np.hstack([0,r])                         # [lunit]
    M0 = np.hstack([0,M])                         # [munit]

    deltaM = M0[1:]-M0[:-1]                       # [munit]
    if min(deltaM)<0.:
        print 'unphysical negative mass increment encountered'
        exit(1)
    
    deltavol = np.pi*(r0[1:]**2 - r0[:-1]**2)        # [lunit^2]
    dens = deltaM/deltavol                           # [munit/lunit^2]
    return dens                                      # [munit/lunit^2]





# TODO: check not used anymore anywhere, use Mr3D() below instead
# def Mr3D(rdat, denspars): #[munit/lunit^3], [lunit]
#     'calculate enclosed mass from integrating up density in spherical shells'
#     M_r = np.zeros(len(rdat)) # [munit]
#     M_r[0] = denspars[0]*4./3.*np.pi*rdat[0]**3 # [munit]
#     for i in range(1,len(rdat)):
#         M_r[i] = M_r[i-1] + denspars[i]*4./3.*np.pi*(rdat[i]**3 - rdat[i-1]**3) # [munit]
#     return M_r  # [munit]





def Mr2D(rdat, dens2D):
    'calculate enclosed mass from integrating up density in rings'
    M_r = np.zeros(len(rdat))
    M_r[0] = dens2D[0]*np.pi*rdat[0]**2
    for i in range(1,len(rdat)):
        M_r[i] = M_r[i-1] + dens2D[i]*np.pi*(rdat[i]**2 - rdat[i-1]**2)
    return M_r






def Mzdefault(denspars):                # [msun/pc^3]
    'get 3D enclosed mass from 3D mass density'
    return Mr3D(gp.xipol, dens(gp.xipol,denspars[:]))






def Mr3D(rdat, denspars):              # [lunit], [munit/lunit^3]
    'take density denspars(rdat), calculate integrated 3D mass'
    # denspars = dens(rdat,denspars)  # [munit/lunit^3] # do not do that twice!
    # done already in gl_funs.calc_M_nu_sig, and Mzdefault() up above

    if gp.densint:
        # TODO: does not exist anymore, either rewrite or delete
        M_r = int_dens(rdat,denspars)
        if gp.cprior<1e29: M_r[0] = gp.cprior

    else : # assuming rho(r) is given as a local measure, in shells
        r0 = np.hstack([0,rdat])                        # [lunit]
        deltavol = 4./3.*np.pi*(r0[1:]**3 - r0[:-1]**3) # [lunit^3]
        deltaM   = denspars * deltavol                  # [munit]
        M_r = np.cumsum(deltaM)                         # [munit]

        # old, explicit way. somehow, this is not the same as above, overestimates M
        # M_r = np.zeros(len(rdat))
        # M_r[0] = denspars[0]*4./3.*np.pi*rdat[0]**3
        # for i in range(1,len(M_r)):
        #     M_r[i] = M_r[i-1] + denspars[i]*4./3.*np.pi*(rdat[i]**3 - rdat[i-1]**3)


    
    # else:                               # assuming rho(r) = M(<r)/V(r)
    #     r0 = rdat
    #     vol = 4./3.*np.pi*r0**3
    #     M_r = denspars * vol
        
    if gp.cprior < 1.e29: M_r[0] = gp.cprior # [munit]
    if (M_r[0] < 0): M_r[0] = 0.        # [munit]

    # no mass can be negative (at least on these scales :o))
    if (min(M_r) < 0):
        gp.LOG.error( 'M sign error! setting M_out_min = 0 for all')
        # gpl.save_plot()
        for i in range(len(M_r)):
            if(M_r[i]<0):
                M_r[i] == 0.0           # [munit]
    return M_r                          # [munit]







def nu(nupars):
    '''General function to describe the density profile.\
    before: Fully non-parametric decreasing function [c.f. Kz function].\
    now:    just return interpolated nupars directly'''
    # drdat = rdat[1:]-rdat[:-1]
    nuout = np.zeros(len(nupars))

    # nuout[0] = nupars[0]
    # for i in range(1,len(rdat)):
    #     nuout[i] = nuout[i-1] * np.exp(-nupars[i]) # * drdat[i-1]

    for i in range(len(nupars)):
        nuout[i] = nupars[i]
    nuout[-1] = 0.9*nuout[-2] # or interpolate via gp.ipol.nudat1[-1]
    # f = gfun.ipol(rdat,nuout,xipol)
    return nuout




def get_zarrays(r_in):
    'set up regularly spaced r0 array for integration, with extension to 2rmax'
    pnts = gp.nipol; rmin = min(r_in); rmax = max(r_in);
    dr = (rmax-rmin)/(pnts-1.); r0 = np.arange(pnts)*dr+rmin
    # arrays of length pnts.   A[-1] means last element of A
    # delta is assumed to be constant from rmax to 2*rmax
    r_outer     = rmax+r0-rmin
    r_tot       = np.hstack((r0[:-1], r_outer))  # [1:] all but first, [:-1] all but last
    return r0,r_tot,rmin,rmax,dr,r_outer





def get_nu(nu_r, r0, r_outer): # [densunit, lunit, lunit]
    'return extrapolated nu to 2rmax'
    pnts = gp.nipol #[1]
    # nu goes down linearly in log space from nu[-1] at rmax
    # towards zero at 2*rmax, not reaching it
    pntsam = gp.nipol/2.05
    nuipolbase = nu_r[::pntsam] # [densunit]
    xipolbase  = r0[::pntsam] #[lunit]
    nu_outer   = gh.ipollog(xipolbase,nuipolbase,r_outer)
    # nu_outer   = gh.ipollog(xipolbase,nuipolbase,r_outer)
    # nu_outer   = gh.ipollog(r0,nu_r,r_outer,smooth=1.e-8)
    # nu_outer   = gh.ipollog(r0,nu_r,r_outer,smooth=10.)
    # nu_outer   = nu_r[-1]*np.ones(pnts-1)
    # nu_outer   = nu_r[-1]/((nu_r[-2]/nu_r[-1])**(np.arange(pnts))) #[densunit], extrapolation from last two points

    if gp.analytic: nu_outer   = rho_anf(r_outer)

    nu_tot = np.hstack((nu_r[:-1], nu_outer)) #[densunit]
    return nu_tot #[densunit]





def get_mprior(M): # [munit]
    if len(M)>gp.nipol:
        print 'too long M in get_mprior!'
        exit(1)

    # choose one of the following lines

    # mprioru  = (M[-1]-M[3*gp.nipol/4])/(gp.xipol[-1]-gp.xipol[3*gp.nipol/4])
    mprioru  = (M[-1]-M[3*gp.nipol/4])/(gp.xipol[-1]-gp.xipol[3*gp.nipol/4]) # tuned variant
    # mprioru = (M[-1]-M[gp.nipol-2])/(gp.xipol[-1]-gp.xipol[gp.nipol-2])
    # mprioru  = 0.

    if not gp.checksigma:
        mprioru = abs(gp.parst.Msl)

    return mprioru #[1]





def get_M(M_r,mprioru,r_outer,dr):
    'get M extrapolated to 2rmax'
    pnts = gp.nipol
    # use slope in last quarter
    # M rises from M[-1] at rmax up to mprioru*M[-1] at 2*rmax)
    M_outer = M_r[-1] + (np.arange(pnts))*mprioru*dr
    if gp.analytic: M_outer = M_anf(r_outer)
    # M_outer = M_r[-1] + (np.arange(1,pnts))*mprioru*M_r[-1]/(rmax-rmin)*dr # for use with mprior
    #M_outer = gfun.ipol(r0,M_r,r_outer[1:])
    return np.hstack((M_r[:-1],     M_outer))




def get_delta(delta_r):                 # [1]
    'get extrapolated delta out to 2rmax'
    delta_outer = delta_r[-1] * np.ones(gp.nipol-1)
    return np.hstack((delta_r,  delta_outer))





def sig_los(pop, x, M_x, nu_r, delta_r): #[1], [pc], [munit, 3D], [dens0,3D], [1]
    'General function to calculate \sigma_{los}:'

    pnts = gp.nipol
    r0,r_tot,rmin,rmax,dr,r_outer = get_zarrays(x) #6*[pc]
    
    # Calculate density and M force:
    mprioru = get_mprior(M_x) #[1], slope

    nu_tot    = get_nu(nu_r,r0,r_outer) #[dens0,3D]
    M_tot     = get_M(M_x,mprioru,r_outer,dr) #[munit,3D]
    delta_tot = get_delta(delta_r) #[1]

    if(gp.analytic):
        sigr2_an  = sig2_anf(r_tot)              #[TODO]
        siglos2surf_an = surfden_sig2_anf(r_tot) #[TODO]
        surfden_an  = surfden_anf(r_tot)         #[TODO]
        siglos_an = sig_los_anf(r_tot)           #[TODO]


    delta_r_t   = int_delta_r_t(pnts, r_tot, delta_tot) # [1]
    
    sigr2       = int_sigr2( pnts, r_tot, delta_r_t, M_tot, nu_tot)
    # takes [1], [pc], [1], [munit], [munit/pc^3], gives back [(km/s)^2]

    siglos2surf = int_siglos2surf(pnts, r_tot, delta_tot, nu_tot, sigr2)
    # takes [1], [pc], [1], [munit/pc^3], [(km/s)^2], gives back # [munit/pc^2 (km/s)^2]

    surfden     = int_surfden( pnts, r_tot, nu_tot)
    # takes [1], [pc], [munit/pc^3], gives back [munit/pc^2]
    

    posi        = siglos2surf/surfden #/1.05    # TODO: check sigr2, remove 1.05 factor
    # takes [munit/pc^2 (km/s)^2], gives back [(km/s)^2]
    
    for i in range(len(posi)):
        if posi[i]<0:
            posi[i] =1e-20              # correct for too low and possibly negative bogus values

    siglos      = np.sqrt(posi)         # [km/s]

    # Interpolate back to input array:
    #siglos_ipol = gfun.smooth(gp.xipol,siglos)

    # if #bins < 30, use following fix to get first bin approx. right
    # TODO: more investigations needed in int_sigr2,surfden to get rid of the following line
    # siglos[0] = siglos[1]*0.95
    return siglos                       # [km/s]
