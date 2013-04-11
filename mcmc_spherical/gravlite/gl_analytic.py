#!/usr/bin/python2.7
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
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
    'equation XX from Hernquist 1990'
    Xout = np.zeros(len(R))
    for i in range(len(R)):
        r = R[i]
        if r<=1.0:
            Xout[i] = (1.0-r**2)**(-0.5)*asech(r)
        else:
            Xout[i] = (r**2-1.0)**(-0.5)*asec(r)
    return Xout



def rho_anf(r):
    'equation XX from Hernquist 1990'
    return 1./(2.*np.pi*r*(1.+r)**3)



def surfden_anf(r):
    'equation XX from Hernquist 1990'
    return ((2.+r**2)*X(r)-3.)/(2.*np.pi*(1.-r**2)**2)




def sig2_anf(r):
    'equation XX from Hernquist 1990'
    return r*(1+r)**3*np.log((1.+r)/r)-\
           (r*(25.+52.*r+42.*r**2+12.*r**3))/(12.*(1.+r))



def M_anf(r):
    'equation XX from Hernquist 1990'
    return r**2/(r+1.)**2




def nusigr2_anf(r):
    'equation XX from Hernquist 1990'
    return 1./(12.*np.pi)*(0.5/(1.-r**2)**3*\
                           (-3.*r**2*X(r)*(8.*r**6-28.*r**4+35.*r**2-20)\
                            -24.*r**6+68.*r**4-65.*r**2+.6)-6.*np.pi*r)



def surfden_sig2_anf(r):
    'equation XX from Hernquist 1990'
    return 1./(24.*np.pi*(1.-r**2)**3)*\
           (3.*r**2*(20.-35.*r**2+28.*r**4-8.*r**6)*X(r)+\
            (6.-65.*r**2+68.*r**4-24.*r**6))-r/2.




def sig_los_anf(r):
    'equation XX from Hernquist 1990'
    return np.sqrt(surfden_sig2_anf(r)/surfden_anf(r))




def walker_delta(pop):
    # TODO: check call of walker_delta, if not used in case investigate != walker, delete if
    if gp.investigate == 'walker':
        from gl_analytic import betawalker
        x = gp.xipol
        delta_r   = betawalker(x)[pop-1]
        # delta_r = gp.delta0
        # delta_r = gfun.ipol(gp.dat.nur1,   gp.delta0, r0) # extrapolate to r_tot?
    else:
        delta_r = gp.delta0

    return delta_r






def rhohern(r0, rscale, rho0, alpha, beta, gamma): # 2*[pc], [munit/pc^3], 3*[1]
    'determine rho(r) for a generalized Hernquist model'
    # rho_DM(r) = rho_0 (r/r_DM)^(-gamma_DM) *
    #             (1+(r/r_DM)^alpha_DM)^((gamma_DM-beta_DM)/alpha_DM)
    tmp = 1.*rho0                       # [munit/pc^3]
    tmp = tmp * (r0/rscale)**(-gamma)   # [munit/pc^3]
    tmp = tmp * (1.+(r0/rscale)**alpha)**((gamma-beta)/alpha) # [munit/pc^3]
    return tmp                          # [munit/pc^3]







def rhowalker_3D(rad):                                          # [pc]
    'Walker model: read values from theoretical params file, calculate from eq. 2 generalized Hernquist profiles'

    # rho_DM(r) = rho_0 (r/r_DM)^(-gamma_DM) *
    #             (1+(r/r_DM)^alpha_DM)^((gamma_DM-beta_DM)/alpha_DM)
    # need values for rho0, r_DM, alpha_DM, beta_DM, gamma_DM
    # and the corresponding variables for the stellar component

    A = np.loadtxt(gp.files.analytic, unpack=False)
    
    gamma_star1 = A[7]
    beta_star1  = A[8]
    alpha_star1 = 2.0
    r_star1     = 1000.*A[9] #[pc] # both description of filename and file content is wrong
    # name is rstar/10 in three digits, file content is given in kpc
    
    gamma_star2 = A[11]
    beta_star2  = A[12]
    alpha_star2 = 2.0
    r_star2     = 1000.*A[13] #[pc]

    gamma_DM    = A[15]
    beta_DM     = A[16]
    r_DM        = A[17] # [pc]
    alpha_DM    = A[18]
    rho0        = A[19] # [Msun/pc^3]

    # TODO: set to right values nu0 (given somewhere?)
    rho0star1   = rho0/1.e6
    rho0star2   = rho0/1.e6

    # rhohern(r, rscale, rho0, alpha, beta, gamma): #  (2*[pc], or 2*[rcore]), [munit/pc^3], 3*[1]
    rhodm    = rhohern(rad, r_DM, rho0, alpha_DM, beta_DM, gamma_DM) # [msun/pc^3]
    rhostar1 = rhohern(rad, r_star1, rho0star1, alpha_star1, beta_star1, gamma_star1) # [msun/pc^3]
    rhostar2 = rhohern(rad, r_star2, rho0star2, alpha_star2, beta_star2, gamma_star2) # [msun/pc^3]

    return rhodm, rhostar1, rhostar2                 # 3* [munit/pc^3]




def rhowalkertot_3D(rad):               # [pc]
    'return total mass density, stars+DM'
    rhodm, rhostar1, rhostar2 = rhowalker_3D(rad)     # 3* [msun/pc^3]
    return rhodm + rhostar1 + rhostar2                # [msun/pc^3]




def Mwalkertot(rad):
    rhotot = rhowalkertot_3D(rad)
    Mtot = phys.Mr3D(rad, rhotot)
    return Mtot





def rhowalker_2D(rad):                  # [pc]
    'calculate 2D surface density from radius'
    rhodm, rhostar1, rhostar2 = rhowalker_3D(rad)     # 3* [msun/pc^3]

    pnts = len(rad)

    surfdm    = int_surfden(pnts, rad, rhodm)    # [msun/pc^2]
    surfstar1 = int_surfden(pnts, rad, rhostar1) # [msun/pc^2]
    surfstar2 = int_surfden(pnts, rad, rhostar2) # [msun/pc^2]

    return surfdm, surfstar1, surfstar2               # 3* [msun/pc^2]



    

def rhowalkertot_2D(rad):                      # [pc]
    'return total surface density, stars+DM'
    surfdm, surfstar1, surfstar2 = rhowalker_2D(rad)   # 3*[msun/pc^2]
    return surfdm + surfstar1 + surfstar2              # [msun/pc^2]





# def betawalker():
#     A = np.loadtxt(gp.files.analytic,unpack=False)
#     beta_star1  = A[8]
#     beta_star2  = A[12]
#     print 'beta1, beta2 =',beta_star1, beta_star2
#     return beta_star1, beta_star2

# this gives velocity anisotropy from model



def betawalker(rad):                                            # [pc]
    'calculate Osipkov-Merritt velocity anisotropy profile'
    A = np.loadtxt(gp.files.analytic,unpack=False)
    # Osipkov-Merritt anisotropy profile with r_a/r_* = 10^4 for isotropic models
    rars1  = A[10]
    rs1    = A[9] * 1000. # [pc]
    rbyrs1 = rad/rs1
    beta1  = rbyrs1**2 / (rars1**2+rbyrs1**2)

    rars2  = A[14]
    rs2    = A[13] * 1000. # [pc]
    rbyrs2 = rad/rs2
    beta2  = rbyrs2**2 / (rars2**2+rbyrs2**2)

    return beta1, beta2
