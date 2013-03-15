#!/usr/bin/python
'''all integrals from physics_sphere'''
import numpy as np
import gl_params as gp
import gl_funs as gfun
import gl_helper as gh
import scipy
import scipy.integrate
from scipy.integrate import simps,trapz


def int_rho(rho,r):
    'Solve mass integral'
    tmp = np.zeros(len(r))
    ri  = np.hstack([0,r,r[-1]+(r[-1]-r[-2])])
    rhoi= gh.ipollog(r,rho,ri)
    # tmp[0] = 4.*np.pi/3.*rho[0]*r[0]**3
    for i in range(2,len(ri)):
        xint = ri[:i]
        yint = rhoi[:i]*ri[:i]**2
        tmp[i-2] = simps(yint,xint,even=gp.even)
    return 4.*np.pi*tmp # *0.92, factor needed to account for overestimation of simpson vs rectangles


# TODO: correct first bin
def int_beta_r_t(pnts, r, b):
    'Solve Beta integral'
    tmp = np.zeros(2*pnts-1)
    for i in range(1,2*pnts-1):
        xint = r[0:i]
        yint = b[0:i]/r[0:i]
        tmp[i] = simps(yint,xint,even=gp.even)
    return tmp

def int_sigr2(pnts, r, b, M, nu):
    'Solve cumulative integral.'
    # the integrand is assumed to be zero above 2*rmax
    tmp = np.zeros(2*pnts-1)
    for i in range(2*pnts-1):
        prefac  = np.exp(-2*b[i])
        xint = r[i:]
        yint = M[i:]*nu[i:] / (r[i:]**2) * np.exp(2*b[i:])
        tmp[i] = prefac * simps(yint,xint,even=gp.even)/nu[i]

    # gp.LOG.info( 'compute sigma outer')
    # could calculate that from last integral already, trying to extrapolate anyhow
    # if (sigmaprior < -1): sigmaprior=-1.
    # if gp.checksigma: sigmaprior=-1.
    # nusig2_outer=nusig2[-1]+(np.arange(1,pnts))*sigmaprior*nusig2[pnts-1]/(xmax-xmin)*dx
    # nusigma2_tot=np.hstack((nusig2,nusig2_outer))
    # sigr2_tot = gh.ipollog(r_tot[:-1],sigr2[:-1],r_tot) # beware of interpolations, dude!
    tmp[-1] = tmp[-2]
    return tmp
    
def int_surfden(pnts, r, nu):
    'compute surface density, assume len(r) = 2*pnts'
    tmp = np.zeros(pnts)
    for i in range(pnts):
        xint = r[i:]
        yint = r[i:]*nu[i:]/np.sqrt(r[i:]**2-r[i]**2)
        # interpolate diverging bin, need 4 points to get 3rd order
        yint[0] = gh.ipol(xint[1:4],yint[1:4],[xint[0]-(xint[1]-xint[0])])
        tmp[i]  = 2.0*simps(yint,xint,even=gp.even)

    # approximate innermost surface density as (total mass inside sphere)/(projected area)
    # based on the assumption that the density inside the innermost bin is much higher than further out
    # this neglects all contributions from before and behind 0 +/- rmin!
    # surfden[0]=4.0/3.0*r_tot[0]**3*nu_tot[0]/\
    #            2.0*r_tot[0]**2
    # better: extrapolation using polynomes
    # surfden = gh.ipollog(r_tot[1:-1],surfden[1:-1],r_tot)
    tmp[0]  = gh.ipol(r[1:3],tmp[1:3],[r[0]])
    return tmp
    
def int_siglos2surf(pnts,r,b,nu,sig2):
    'Compute line of sight integral'
    tmp = np.zeros(pnts)
    for i in range(pnts):
        xint = r[i:]
        yint = (1-b[i:]*(r[i]/r[i:])**2) * nu[i:] * sig2[i:]*r[i:] / np.sqrt(r[i:]**2-r[i]**2)
        yint[0] = gh.ipol(xint[1:3],yint[1:3],[xint[0]])#-(xint[1]-xint[0]))        
        tmp[i] = 2.0*simps(yint,xint,even=gp.even)  #/surfden

    # crude approximation: first bin has same sigma_los2 as second bin
    # tmp[0] = tmp[1]
    # tmp[0] = gh.ipol(r_tot[1:],tmp[1:],r_tot[0])
    # tmp[0] = 2.0/surfden[0]*nu_tot[0:]*sigma_tot[0:]*r_tot[0:]/r_tot[:]*rmin
    return tmp
