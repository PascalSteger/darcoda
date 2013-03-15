#!/usr/bin/python
import pdb
import numpy as np
import gl_params as gp
import gl_funs as gfun
import gl_plot as gplot
from gl_int import *
from gl_analytic import *


def rhodm_hernquist(r,rho0,r_DM,alpha_DM,beta_DM,gamma_DM):
    return rho0*(r/r_DM)**(-gamma_DM)*(1+(r/r_DM)**alpha_DM)**((gamma_DM-beta_DM)/alpha_DM)

def Mzdefault(denspars):
    return Mr(denspars, gp.xipol)

def densdefault(denspars):
    return dens(denspars,gp.xipol)

# def dens(denspars, xipol):
#     tmp = np.zeros(len(xipol))
#     # tmp += denspars[0] # (*xipol**0) # TODO: check unnecessary, included below
#     for i in range(0,len(denspars)):
#         tmp += denspars[i]*(xipol)**(-i)
#     # gplot.plot(gp.ipol.densr,gp.ipol.densdat); gplot.plot(gp.xipol,10.**tmp); gplot.yscale('log')
#     return 10.0**tmp

def dens(denspars, xipol):
    if not gp.poly:
        return denspars
    if gp.checksigma:
        return gp.ipol.densdat

    tmp = np.zeros(len(xipol))
    for i in range(0,len(denspars)):
        tmp += denspars[i]*((gp.scaledens-xipol)/gp.scaledens)**i
    # gplot.plot(gp.ipol.densr,gp.ipol.densdat); gplot.plot(gp.xipol,10.**tmp); gplot.yscale('log')
    return 10.0**tmp

def calculate_dens(M,r): #[totmass], [rcore]
    r0 = np.hstack([0,r]) #[rcore]
    deltaM = [];
    deltaM.append(M[0]) #[totmass]
    for i in range(1,len(r)):
        deltaM.append(M[i]-M[i-1]) #[totmass]
    
    dens = [];  deltavol = []
    for i in range(len(r)):
        deltavol.append( 4./3.*np.pi*(r0[i+1]**3-r0[i]**3) ) #[rcore**3]
        dens.append( deltaM[i]/deltavol[i] ) #[totmass/rcore**3]
    return np.array(dens) #[totmass/rcore**3]

def Mr(denspars, rdat):
    'General function to calculate enclosed mass M from density denspars(rdat)'
    if gp.poly:
        denspars = dens(denspars,rdat)

    if gp.densint:
        M_r = int_dens(denspars,rdat)
        if gp.cprior<1e29: M_r[0] = gp.cprior
    else:
        M_r = np.zeros(len(rdat))
        M_r[0] = denspars[0]*4./3.*np.pi*rdat[0]**3
        for i in range(1,len(M_r)):
            M_r[i] = M_r[i-1] + denspars[i]*4./3.*np.pi*(rdat[i]**3 - rdat[i-1]**3)
        if gp.cprior < 1.e29: M_r[0] = gp.cprior

#     if(set(rdat)==set(xipol)):
#         M_out = gfun.smoothlog(xipol,M_r)
#     else:
#         M_out = gfun.ipollog(rdat,M_r,xipol)
    M_out = M_r
    
    # Sometimes when kz_z(0) << kz_z(1), the interpolant goes negative. This just means that 
    # we have resolved what should be kz_z(0) = 0.
    if (M_out[0] < 0): M_out[0] = 0.

    # This should never happen: 
    if (min(M_out) < 0):
        gp.LOG.error( 'Kz sign error! setting M_out_min = 0 for all')
        # gplot.save_plot()
        for i in range(len(M_out)):
            if(M_out[i]<0):
                M_out[i]==0.0
    return M_out

def nu(nupars):
    'General function to describe the density profile.\
    before: Fully non-parametric decreasing function [c.f. Kz function].\
    now:    just return interpolated nupars directly'
    # drdat = rdat[1:]-rdat[:-1]
    # nuout = np.zeros(len(nupars))
    # nuout[0] = nupars[0]
    # for i in range(1,len(rdat)):
    #     nuout[i] = nuout[i-1] * np.exp(-nupars[i]) # * drdat[i-1]

    nuout = nupars
    nuout[-1] = 0.9*nuout[-2] # or interpolate via gp.ipol.nudat1[-1]
    # f = gfun.ipol(rdat,nuout,xipol)
    return nuout

def get_zarrays(r_in):
    # Set up regularly spaced r0 array for integration:
    pnts = gp.nipol; rmin = min(r_in); rmax = max(r_in);
    dr = (rmax-rmin)/(pnts-1.); r0 = np.arange(pnts)*dr+rmin
    # arrays of length pnts.   A[-1] means last element of A
    # delta is assumed to be constant from rmax to 2*rmax
    r_outer     = rmax+r0-rmin
    r_tot       = np.hstack((r0[:-1], r_outer))  # [1:] all but first, [:-1] all but last
    return r0,r_tot,rmin,rmax,dr,r_outer

def get_nu(nu_r,r0,r_outer):
    pnts = gp.nipol
    # nu goes down linearly in log space from nu[-1] at rmax
    # towards zero at 2*rmax, not reaching it
    pntsam = gp.nipol/2.05
    nuipolbase = nu_r[::pntsam]
    xipolbase  = r0[::pntsam]
    # nu_outer   = gfun.ipol(xipolbase,nuipolbase,r_outer)
    # nu_outer   = gfun.ipollog(xipolbase,nuipolbase,r_outer)
    nu_outer   = nu_r[-1]/((nu_r[-2]/nu_r[-1])**(np.arange(pnts)))
    # nu_outer = gfun.ipollog(r0,nu_r,r_outer,smooth=10.)
    if gp.analytic: nu_outer   = rho_anf(r_outer)
    # nu_outer   = nu_r[-1]*np.ones(pnts-1)
    # nu_outer   = gfun.ipollog(r0,nu_r,r_outer,smooth=1)
    nu_tot = np.hstack((nu_r[:-1], nu_outer))
    return nu_tot

def get_M(M_r,mprioru,r_outer,dr):
    pnts = gp.nipol
    # use slope in last quarter
    # M rises from M[-1] at rmax up to mprioru*M[-1] at 2*rmax)
    M_outer = M_r[-1] + (np.arange(pnts))*mprioru*dr
    if gp.analytic: M_outer = M_anf(r_outer)
    # M_outer = M_r[-1] + (np.arange(1,pnts))*mprioru*M_r[-1]/(rmax-rmin)*dr # for use with mprior
    #M_outer = gfun.ipol(r0,M_r,r_outer[1:])
    return np.hstack((M_r[:-1],     M_outer))

def get_delta(delta_r):
    delta_outer = delta_r[-1]*np.ones(gp.nipol-1)
    return np.hstack((delta_r,  delta_outer))

def sig_los(pop):
    'General function to calculate \sigma_{los}:'
    pnts = gp.nipol
    r0,r_tot,rmin,rmax,dr,r_outer = get_zarrays(gp.xipol)
    
    # Calculate density and M force:
    if gp.checksigma:
        M_r    = gp.ipol.Mdat # or Mzdefault(gp.ipol.densdat) to check integration
        if pop==1:
            nu_r   = gp.ipol.nudat1
        elif pop==2:
            nu_r   = gp.ipol.nudat2
        if gp.analytic:
            M_r    = M_anf(r0)
            nu_r   = rho_anf(r0)
        if gp.investigate == 'walker':
            delta_r   = betawalker(gp.xipol)[pop-1] #gp.delta0 # gfun.ipol(gp.dat.nur1,   gp.delta0, r0)
        else:
            delta_r = gp.delta0
        mprioru  = (M_r[pnts-1]-M_r[3*pnts/4])/(r0[pnts-1]-r0[3*pnts/4]) # tuned variant
        # mprioru = (M_r[pnts-1]-M_r[pnts-2])/(r0[pnts-1]-r0[pnts-2])
    else:
        M_r     = Mzdefault(gp.parst.dens)
        mprioru = abs(gp.parst.Msl)
        if pop==1:
            nu_r   = nu(gp.parst.nu1)
            delta_r = gp.parst.delta1
        elif pop==2:
            # change name to include both, for two runs
            nu_r   = nu(gp.parst.nu2)
            delta_r = gp.parst.delta2

    nu_tot   = get_nu(nu_r,r0,r_outer)
    M_tot    = get_M(M_r,mprioru,r_outer,dr)
    delta_tot = get_delta(delta_r)

    if(gp.analytic):
        sigr2_an  = sig2_anf(r_tot)
        siglos2surf_an = surfden_sig2_anf(r_tot)
        surfden_an  = surfden_anf(r_tot)
        siglos_an = sig_los_anf(r_tot)

    delta_r_t   = int_beta_r_t(pnts, r_tot, delta_tot)
    sigr2       = int_sigr2( pnts, r_tot, delta_r_t, M_tot, nu_tot)
    siglos2surf = int_siglos2surf(pnts, r_tot, delta_tot, nu_tot, sigr2)
    surfden     = int_surfden( pnts, r_tot, nu_tot)
    # TODO: check sigr2, remove 1.05 factor
    posi        = siglos2surf/surfden #/1.05
    for i in range(len(posi)):
        if posi[i]<0:
            posi[i] =1e-20
    siglos      = np.sqrt(posi)

    # Interpolate back to input array:
    #siglos_ipol = gfun.smooth(gp.xipol,siglos)

    # if #bins < 30, use following fix to get first bin approx. right
    # TODO: more investigations needed in int_sigr2 to get rid of this
    # siglos[0] = siglos[1]*0.95 # TODO: remove, use better interpolation above for sigr2, surfden
    return siglos
