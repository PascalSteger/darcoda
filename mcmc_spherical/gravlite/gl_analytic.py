#!/usr/bin/python
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


def X(s0):                              # [1]
    '''equation 33, 34 from Hernquist 1990'''
    Xout = np.zeros(len(s0))
    for i in range(len(s0)):
        s = s0[i]                       # [1]
        
        if s<=1.:
            Xout[i] = (1.-s**2)**(-0.5)*asech(s) # [1]
        else:
            Xout[i] = (s**2-1.)**(-0.5)*asec(s) # [1]
        
    return Xout                         # [1]




def rho_anf(r0, a=gp.ascale, M=gp.Mscale): # 2*[pc], [Msun]
    '''equation 2b from Baes&Dejonghe 2002'''
    '''equation 2 from Hernquist 1990, 3D mass density'''
    s = r0/a                              # [1]
    return M/(2.*np.pi*a**3)/(s*(1+s)**3) # [Msun/pc^3]



def surfden_anf(r0, a=gp.ascale, M=gp.Mscale): # 2*[pc], [Msun]
    '''equation 3 from Baes&Dejonghe 2002'''
    s = r0/a
    return M/(2.*np.pi*a**2)*((2.+s**2)*X(s)-3.)/((1.-s**2)**2) # [Msun/pc^2]




def sig2_anf(r0, a=gp.ascale, M=gp.Mscale): # 2*[pc], [Msun]
    '''equation 10 of Hernquist 1990, sigma_r^2'''
    s = r0/a                            # [pc/1000pc]
    return gp.G1*M/a*s*(1+s)**3*np.log(1.+1./s)-\
           gp.G1*M/(12.*a)*s/(s+1)*(25.+52.*s+42.*s**2+12.*s**3) # [(km/s)^2]




def M_anf(r0, a=gp.ascale, M=gp.Mscale): # 2*[pc], [Msun]
    '''equation 3 from Hernquist 1990'''
    s = r0/a                            # [1]
    return M*s**2/(1.+s)**2             # [Msun]






def surfden_sig2_anf(r0, a=gp.ascale, M=gp.Mscale): # 2*[pc], [Msun]
    '''equation 21 from Baes&Dejonghe 2002'''
    s = r0/a                            # [1]
    return gp.G1*M**2/(12.*np.pi*a**3)*(1./(2.*(1.-s**2)**3)*(-3.*s**2*X(s)*(8.*s**6-28.*s**4+35.*s**2-20.)-24.*s**6+68.*s**4-65.*s**2+6.)-6.*np.pi*s) # [(km/s)^2 * Msun/pc^2]




def sig_los_anf(r0, a=gp.ascale, M=gp.Mscale):   # 2*[pc], [Msun]
    '''determined from analytic surfden*sig2 and surfden'''
    return np.sqrt(surfden_sig2_anf(r0,a,M)/surfden_anf(r0,a,M))




def walker_delta(pop):
    # TODO: check call of walker_delta, if not used in case investigate != walker: delete 'if'
    if gp.investigate == 'walker':
        x = gp.xipol
        delta_r   = betawalker(x)[pop-1]
        # delta_r = gp.delta0
        # delta_r = gfun.ipol(gp.dat.nur1,   gp.delta0, r0) # extrapolate to r_tot?
    else:
        delta_r = gp.delta0

    return delta_r






def rhohern(r0, rscale, rho0, alpha, beta, gamma): # 2*[pc], [munit/pc^3], 3*[1]
    '''determine rho(r) for a generalized Hernquist model'''
    # rho_DM(r) = rho_0 (r/r_DM)^(-gamma_DM) *
    #             (1+(r/r_DM)^alpha_DM)^((gamma_DM-beta_DM)/alpha_DM)
    tmp = 1.*rho0                       # [munit/pc^3]
    tmp = tmp * (r0/rscale)**(-gamma)   # [munit/pc^3]
    tmp = tmp * (1.+(r0/rscale)**alpha)**((gamma-beta)/alpha) # [munit/pc^3]
    return tmp                          # [munit/pc^3]





def rhotriax(rad):
    alpha = 1.
    beta = 4.
    rs = 1500.                          # [pc]
    if gp.triaxcase == 0:               # core
        gamma = 0.23
        rhos  = 1.177E-1                # [Msun/pc^3]
    elif gp.triaxcase == 1:             # cusp
        gamma = 1.
        rhos  = 5.522E-2                # [Msun/pc^3]

    rho = rhos
    rho /= (rad/rs)**gamma
    rho /= (1+(rad/rs)**(1/alpha))**(alpha*(beta-gamma))
    return rho





def rhowalker_3D(rad):                                          # [pc]
    '''Walker model: read values from theoretical params file, calculate from eq. 2 generalized Hernquist profiles'''

    # rho_DM(r) = rho_0 (r/r_DM)^(-gamma_DM) *
    #             (1+(r/r_DM)^alpha_DM)^((gamma_DM-beta_DM)/alpha_DM)
    # need values for rho0, r_DM, alpha_DM, beta_DM, gamma_DM
    # and the corresponding variables for the stellar component

    A = np.loadtxt(gp.files.analytic, unpack=False)
    
    gamma_star1 = A[7]
    beta_star1  = A[8]
    alpha_star1 = 2.0
    r_star1     = 1000*A[9] #[pc] # both description of filename and file content is wrong
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

    # TODO: set to right values nu0 (given somewhere?) tracer numbers?
    ntracer1 = np.loadtxt(gp.files.get_ntracer_file(1),unpack=True)
    ntracer2 = np.loadtxt(gp.files.get_ntracer_file(2),unpack=True)
    rho0star1   = rho0*ntracer1/1.e6
    rho0star2   = rho0*ntracer2/1.e6

    # rhohern(r, rscale, rho0, alpha, beta, gamma): #  (2*[pc], or 2*[rcore]), [munit/pc^3], 3*[1]
    rhodm    = rhohern(rad, r_DM, rho0, alpha_DM, beta_DM, gamma_DM) # [msun/pc^3]
    rhostar1 = rhohern(rad, r_star1, rho0star1, alpha_star1, beta_star1, gamma_star1) # [msun/pc^3]
    rhostar2 = rhohern(rad, r_star2, rho0star2, alpha_star2, beta_star2, gamma_star2) # [msun/pc^3]

    return rhodm, rhostar1, rhostar2                 # 3* [munit/pc^3]




def rhowalkertot_3D(rad):               # [pc]
    '''return total mass density, stars+DM'''
    rhodm, rhostar1, rhostar2 = rhowalker_3D(rad)     # 3* [msun/pc^3]
    return rhodm + rhostar1 + rhostar2                # [msun/pc^3]




def Mwalkertot(rad):
    '''return total mass in Walker model'''
    rhotot = rhowalkertot_3D(rad)
    Mtot = phys.Mr3D(rad, rhotot)
    return Mtot





def rhowalker_2D(rad):                  # [pc]
    '''calculate 2D surface density from radius'''
    rhodm, rhostar1, rhostar2 = rhowalker_3D(rad)     # 3* [msun/pc^3]

    surfdm    = int_surfden(rad, rhodm)    # [msun/pc^2]
    surfstar1 = int_surfden(rad, rhostar1) # [msun/pc^2]
    surfstar2 = int_surfden(rad, rhostar2) # [msun/pc^2]

    return surfdm, surfstar1, surfstar2               # 3* [msun/pc^2]



    

def rhowalkertot_2D(rad):                      # [pc]
    '''return total surface density, stars+DM'''
    surfdm, surfstar1, surfstar2 = rhowalker_2D(rad)   # 3*[msun/pc^2]
    return surfdm + surfstar1 + surfstar2              # [msun/pc^2]





# def betawalker():
#     A = np.loadtxt(gp.files.analytic,unpack=False)
#     beta_star1  = A[8]
#     beta_star2  = A[12]
#     print 'beta1, beta2 =',beta_star1, beta_star2
#     return beta_star1, beta_star2

# this gives velocity anisotropy from model


def betatriax(rad):
    '''calculate velocity anisotropy for Dehnen-Wilkinson triaxial model'''
    eta = 0.5
    rsbeta = 810.                       # [pc]
    beta0 = 0.
    betainf = 0.5
    beta = (rsbeta**eta*beta0+rad**eta*betainf)/(rad**eta+rsbeta**eta)
    return beta



def betawalker(rad):                                            # [pc]
    '''calculate Osipkov-Merritt velocity anisotropy profile'''
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
