#!/usr/bin/python
'''all analytic profiles from physics_sphere'''
import numpy as np
from gl_int import int_surfden
#from physics_sphere import get_nu,get_zarrays
import physics_sphere as phys
import gl_params as gp
import pdb

asech = lambda x: np.arccosh(1./x)
asec  = lambda x: np.arccos(1./x)

def X(R):
    Xout = np.zeros(len(R))
    for i in range(len(R)):
        r = R[i]
        if r<=1.0:
            Xout[i] = (1.0-r**2)**(-0.5)*asech(r)
        else:
            Xout[i] = (r**2-1.0)**(-0.5)*asec(r)
    return Xout

def rho_anf(r):
    return 1./(2.*np.pi*r*(1.+r)**3)

def surfden_anf(r):
    return ((2.+r**2)*X(r)-3.)/(2.*np.pi*(1.-r**2)**2)

def sig2_anf(r):
    return r*(1+r)**3*np.log((1.+r)/r)-\
           (r*(25.+52.*r+42.*r**2+12.*r**3))/(12.*(1.+r))

def M_anf(r):
    return r**2/(r+1.)**2

def nusigr2_anf(r):
    return 1./(12.*np.pi)*(0.5/(1.-r**2)**3*\
                           (-3.*r**2*X(r)*(8.*r**6-28.*r**4+35.*r**2-20)\
                            -24.*r**6+68.*r**4-65.*r**2+.6)-6.*np.pi*r)
def surfden_sig2_anf(r):
    return 1./(24.*np.pi*(1.-r**2)**3)*\
           (3.*r**2*(20.-35.*r**2+28.*r**4-8.*r**6)*X(r)+\
            (6.-65.*r**2+68.*r**4-24.*r**6))-r/2.

def sig_los_anf(r):
    return np.sqrt(surfden_sig2_anf(r)/surfden_anf(r))

def rhohern(rho0, alpha, beta, gamma, rDM, r): # [munit/pc^3], [pc], [pc]
    # determine M(r) for a generalized Hernquist model
    # rho_DM(r) = rho_0 (r/r_DM)^(-gamma_DM) *
    #             (1+(r/r_DM)^alpha_DM)^((gamma_DM-beta_DM)/alpha_DM)
    tmp = 1.*rho0
    tmp *= (r/rDM)**(-gamma)
    tmp *= (1.+(r/rDM)**alpha)**((gamma-beta)/alpha)
    return tmp
    
def rhowalker(rad): #[munit/pc^2], [pc]
    # Walker model: read values from theoretical params file, calculate from eq. 2
    # generalized Hernquist profiles

    # rho_DM(r) = rho_0 (r/r_DM)^(-gamma_DM) *
    #             (1+(r/r_DM)^alpha_DM)^((gamma_DM-beta_DM)/alpha_DM)
    # need values for rho0, r_DM, alpha_DM, beta_DM, gamma_DM
    # and the corresponding variables for the stellar component

    A = np.loadtxt(gp.files.analytic,unpack=False)
    gamma_star1 = A[7]
    beta_star1  = A[8]
    alpha_star1 = 2.0
    r_star1     = A[9] #[pc]
    
    gamma_star2 = A[11]
    beta_star2  = A[12]
    alpha_star2 = 2.0
    r_star2     = A[13] #[pc]

    gamma_DM    = A[15]
    beta_DM     = A[16]
    r_DM        = A[17] # [pc]
    alpha_DM    = A[18]
    rho0        = A[19] # [Msun/pc^3]
    
    rhostar1 = rhohern(rho0, alpha_star1, beta_star1, gamma_star1, r_star1, rad) #[munit/pc^3]
    rhostar2 = rhohern(rho0, alpha_star2, beta_star2, gamma_star2, r_star2, rad) #[munit/pc^3]
    rhodm    = rhohern(rho0, alpha_DM, beta_DM, gamma_DM, r_DM, rad) #[munit/pc^3]

    # calculate 2D surface density from this    
    surfdm, surfstar1, surfstar2 = surf_walker(gp.nipol, rad, rhodm, rhostar1, rhostar2)

    return surfdm, surfstar1, surfstar2 #[munit/pc^2]


#surf_dm, surf_star1, surf_star2 = surf_walker(pnts, r, M):
def surf_walker(pnts, r, rhodm, rhostar1, rhostar2):
    r0,r_tot,dummy,dummy,dr,r_outer = phys.get_zarrays(r)
    rhodm_tot    = phys.get_nu(rhodm,    r0, r_outer)
    rhostar1_tot = phys.get_nu(rhostar1, r0, r_outer)
    rhostar2_tot = phys.get_nu(rhostar2, r0, r_outer)
    surf_dm    = int_surfden(pnts, r_tot, rhodm_tot)
    surf_star1 = int_surfden(pnts, r_tot, rhostar1_tot)
    surf_star2 = int_surfden(pnts, r_tot, rhostar2_tot)
    return surf_dm, surf_star1, surf_star2

# def betawalker():
#     A = np.loadtxt(gp.files.analytic,unpack=False)
#     beta_star1  = A[8]
#     beta_star2  = A[12]
#     print 'beta1, beta2 =',beta_star1, beta_star2
#     return beta_star1, beta_star2

# this gives velocity anisotropy from model

def betawalker(rad): # rad in pc
    A = np.loadtxt(gp.files.analytic,unpack=False)
    # Osipkov-Merritt anisotropy profile with r_a/r_* =
    rars1  = A[10]
    rs1    = A[9]  # pc
    rbyrs1 = rad/rs1
    beta1  = rbyrs1**2/(rars1**2+rbyrs1**2)

    rars2  = A[14]
    rs2    = A[13] # pc
    rbyrs2 = rad/rs2
    beta2  = rbyrs2**2/(rars2**2+rbyrs2**2)
    
    return beta1, beta2
