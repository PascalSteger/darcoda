#!/usr/bin/env python3

##
# @file
# @ingroup gravlite
# all analytic profiles from gl_physics

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import gl_params as gp
from gl_project import rho_INT_Rho, rho_SUM_Mr
import pdb

asech = lambda x: np.arccosh(1./x)
asec  = lambda x: np.arccos(1./x)

def X(s0):
    Xout = np.zeros(len(s0))
    for i in range(len(s0)):
        s = s0[i]                       # [1]
        
        if s<=1.:
            Xout[i] = (1.-s**2)**(-0.5)*asech(s) # [1]
        else:
            Xout[i] = (s**2-1.)**(-0.5)*asec(s) # [1]
        
    return Xout                         # [1]
## \fn X(s0)
# equation 33, 34 from Hernquist 1990
# @param s0: radius, [1]


def rho_anf(r0, a=gp.ascale, M=gp.Mscale):
    s = r0/a                              # [1]
    return M/(2.*np.pi*a**3)/(s*(1+s)**3) # [Msun/pc^3]
## \fn rho_anf(r0, a=gp.ascale, M=gp.Mscale)
# equation 2 from Hernquist 1990, 3D mass density
# equation 2b from Baes&Dejonghe 2002
# @param r0 radius in [pc]
# @param a scale radius in [pc]
# @param M scale mass in [Msun]
# @return density in [Msun/pc^3]

def rho_2_anf(r0, a, M):
    gamma = 1. # TODO: change if other than gamma_DM = 1 used!
    rho0 = M*(3.-gamma)/(4.*np.pi*a**3)
    rho0 *= (r0/a)**(-gamma)
    rho0 *= (1.+r0/a)**(gamma-4.)
    return rho0
## \fn rho_2_anf(r0, a=gp.ascale, M=gp.Mscale)
# find theoretical mass profile of double-Dehnen spheres with stars,
# according to http://www.astrosim.net/code/doku.php?id=home:codetest:massmodel:sphericaltest
# @param r0 radius in pc
# @param a  scale radius of the overall mass distribution in pc
# @param M  M_\infty of all stars+DM

def Sigma_anf(r0, a=gp.ascale, M=gp.Mscale):
    s = np.array(r0)/a                                          # [1]
    return M/(2.*np.pi*a**2)*((2.+s**2)*X(s)-3.)/((1.-s**2)**2) # [Msun/pc^2]
## \fn Sigma_anf(r0, a=gp.ascale, M=gp.Mscale)
# surface density of Hernquist profile, equation 3 from Baes&Dejonghe 2002
# @param r0 radius in [pc]
# @param a scale radius in [pc]
# @param M scale mass in [Msun]
# @return 2D surface density in [Msun/pc^2]


def sigr2_anf(r0, a=gp.ascale, M=gp.Mscale): 
    s = np.array(r0)/a          # [pc/1000pc]
    return gp.G1*M/a*s*(1+s)**3*np.log(1.+1./s)-\
           gp.G1*M/(12.*a)*s/(s+1)*(25.+52.*s+42.*s**2+12.*s**3) # [(km/s)^2]
## \fn sig_r_2_anf(r0, a=gp.ascale, M=gp.Mscale)
# sigma_r^2, equation 10 of Hernquist 1990
# @param r0 radius in [pc]
# @param a scale radius in [pc]
# @param M scale mass in [Msun]
# @return sigma_r^2 in [(km/s)^2]


def kappa_anf(r0):
    return 3.*np.ones(len(r0))          # [1]
## \fn kappa_anf(r0)
# analyitc value for the fourth order momentum of the LOS velocity
# @param r0 radius in [pc]
# @return 3 for Gaussian velocity distribution

def M_anf(r0, a=gp.ascale, M=gp.Mscale): 
    s = r0/a                            # [1]
    return M*s**2/(1.+s)**2             # [Msun]
## \fn M_anf(r0, a=gp.ascale, M=gp.Mscale)
# equation 3 from Hernquist 1990
# @param r0 radius in [pc]
# @param a scale radius in [pc]
# @param M scale mass in [Msun]
# @return 3D mass in [Msun]


    
def Sigma_sig_los_2_anf(r0, a=gp.ascale, M=gp.Mscale):
    # \sigma_p = \sigma_projected = \sigma_{LOS}
    s = r0/a                            # [1]
    return gp.G1*M**2/(12.*np.pi*a**3)*(1./(2.*(1.-s**2)**3)\
                                          *(-3.*s**2*X(s)\
                                          *(8.*s**6-28.*s**4+35.*s**2-20.)\
                                          -24.*s**6+68.*s**4-65.*s**2+6.)\
                                          -6.*np.pi*s) # [(km/s)^2 * Msun/pc^2]
## \fn Sigma_sig_los_2_anf(r0, a=gp.ascale, M=gp.Mscale)
# equation 21 from Baes&Dejonghe 2002
# @param r0 radius in [pc]
# @param a scale radius in [pc]
# @param M scale mass in [Msun]
# @return surface density * sigma_LOS^2 in [(km/s)^2 * Msun/pc^2]


def sig_los_anf(r0, a=gp.ascale, M=gp.Mscale):
    return np.sqrt(Sigma_sig_los_2_anf(r0,a,M)/Sigma_anf(r0,a,M))
## \fn sig_los_anf(r0, a=gp.ascale, M=gp.Mscale)
# sig_los determined from analytic surfden*sig2 and surfden
# @param r0 radius in [pc]
# @param a scale radius in [pc]
# @param M scale mass in [Msun]
# @return sigma_LOS in [km/s]


def walker_delta(pop, r0):
    r0
    delta_r   = betawalker(r0)[pop-1]
    # delta_r = gp.delta0
    # delta_r = gh.ipol(gp.dat.nur1,   gp.delta0, r0) # extrapolate to r_tot?
    return delta_r
## \fn walker_delta(pop)
# analytic value for the Osipkov-Merrit velocity anisotropy as defined in the Walker dataset
# @param pop which population?
# @param r0 radii for evaluation
# @return beta in [1], from 0 up to 1


def rhohern(r0, rscale, rho0, alpha, beta, gamma):
    # rho_DM(r) = rho_0 (r/r_DM)^(-gamma_DM) *
    #             (1+(r/r_DM)^alpha_DM)^((gamma_DM-beta_DM)/alpha_DM)
    tmp = 1.*rho0                       # [munit/pc^3]
    tmp = tmp * (r0/rscale)**(-gamma)   # [munit/pc^3]
    tmp = tmp * (1.+(r0/rscale)**alpha)**((gamma-beta)/alpha) # [munit/pc^3]
    return tmp                          # [munit/pc^3]
## \fn rhohern(r0, rscale, rho0, alpha, beta, gamma):
# determine rho(r) for a generalized Hernquist model
# @param r0 radius in [pc]
# @param rscale scale radius in [pc]
# @param rho0 central 3D mass density, in [Msun/pc^3]
# @param gamma inner slope TODO: check
# @param alpha outer slope TODO: check
# @param beta steepness of turnover TODO: check
# @return 3D Hernquist density in [Msun/pc^3]

    
def rhotriax(rad):
    alpha = 1.
    beta = 4.
    rs = 1500.                          # [pc]
    if gp.case == 0:               # core
        gamma = 0.23
        rhos  = 1.177E-1                # [Msun/pc^3]
    elif gp.case == 1:             # cusp
        gamma = 1.
        rhos  = 5.522E-2                # [Msun/pc^3]

    rho = rhos
    rho /= (rad/rs)**gamma
    rho /= (1+(rad/rs)**(1/alpha))**(alpha*(beta-gamma))
    return rho
## \fn rhotriax(rad)
# density for a triaxial halo in Gaia challenge
# @param rad radius in [pc]
# @return density in [Msun/pc^3]


def rhowalk_3D(rad):
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

    ntracer1 = np.loadtxt(gp.files.get_ntracer_file(1)+'_3D',unpack=True)
    # [msun/pc^3]
        
    # 1)
    # rho0star1   = rho0/1.e6*ntracer1 #*10 * ntracer2 # too high, too low
    # 2)
    # rho0star1   = rho0/1.e6 #*gp.G1*np.pi*1.5 # *ntracer1/(ntracer1+ntracer2)
    # rho0star2   = rho0/1.e6 #*gp.G1*np.pi/2. # *ntracer2/(ntracer1+ntracer2)
    # 3)
    rho0star1    = rho0/1.e6*ntracer1

    # rhohern(r, rscale, rho0, alpha, beta, gamma): 
    #  (2*[pc], or 2*[rcore]), [munit/pc^3], 3*[1]
    rhodm    = rhohern(rad, r_DM, rho0, alpha_DM, beta_DM, gamma_DM) # [msun/pc^3]
    rhostar1 = rhohern(rad, r_star1, rho0star1, alpha_star1, beta_star1, gamma_star1)

    # [msun/pc^3]

    # TODO: adjust for pops=1 case: only using rho0star1
    ntracer2  = np.loadtxt(gp.files.get_ntracer_file(2)+'_3D', unpack=True)
    rho0star2 = rho0/1.e6*ntracer2
    rhostar2  = rhohern(rad, r_star2, rho0star2, alpha_star2, beta_star2, gamma_star2)

    return rhodm, rhostar1, rhostar2                 # 3* [munit/pc^3]
## \fn rhowalk_3D(rad)
# Walker model: read values from theoretical params file,
# calculate from eq. 2 generalized Hernquist profiles
# @param rad radius in [pc]
# @return 3D density in [Msun/pc^3], for each DM, stellar pop 1, stellar pop 2


def rhogaiatot_3D(rad):
    alpha_star1 = 2.
    alpha_DM = 1.;    beta_DM = 3.
    beta_star1,r_DM,gamma_star1,r_star1,r_a1,gamma_DM,rho0 = gp.files.params
    rhodm = rhohern(rad, r_DM, rho0, alpha_DM, beta_DM, gamma_DM)
    rhostar1 = rhohern(rad, r_star1, rho0/1.e6, alpha_star1, beta_star1, gamma_star1)
    return rhodm+rhostar1
## \fn rhogaiatot_3D(rad)
# give total density for the spherical Gaia challenge dataset
# @param rad radius in [pc]
# @return 3D overall density (DM+stellar population) in [Msun/pc^3]


def nugaiatot_3D(rad):
    alpha_star1 = 2.
    alpha_DM = 1.;    beta_DM = 3.
    beta_star1, r_DM, gamma_star1, r_star1, r_a1, gamma_DM, rho0 = gp.files.params
    rhostar1 = rhohern(rad, r_star1, rho0, alpha_star1, beta_star1, gamma_star1)
    return rhostar1
## \fn rhogaiatot_3D(rad)
# give total density for the spherical Gaia challenge dataset
# @param rad radius in [pc]
# @return 3D overall density (DM+stellar population) in [Msun/pc^3]


def nrwalktot_3D_deriv(rad):
    lrho = np.log(rhowalktot_3D(rad))
    lr   = np.log(rad)
    import gl_helper as gh
    return -gh.derivcoarse(lrho, lr)
## \fn nrwalktot_3D(rad)
# plot d log rho/d log r
# @param rad radius in pc, not in log
    

def nrgaiatot_3D_deriv(rad):
    lrho = np.log(rhogaiatot_3D(rad))
    lr   = np.log(rad)
    import gl_helper as gh
    return -gh.derivcoarse(lrho, lr)
## \fn nrgaiatot_3D(rad)
# plot d log rho/d log r
# @param rad radius in pc, not in log


def nrwalktot_3D(rad):
    # TODO: find a way to get at these values for walk investigation
    beta_star1, r_DM, gamma_star1, r_star1, r_a1, gamma_DM, rho0 = gp.files.params
    if gamma_DM == 0:
        nr = 3*rad/(r_DM+rad)           # TODO: update formula
    elif gamma_DM == 1:
        nr = 2*rad/(r_DM+rad)+1         # TODO: update formula
    else:
        print('unknown gamma_DM = ', gamma_DM)
        nr = 0.*rad
    return nr
## \fn nrwalktot_3D(rad)
# plot - d log rho/ d log r for Walker models
# @param rad radius in pc, not log.


def nrgaiatot_3D(rad):
    beta_star1, r_DM, gamma_star1, r_star1, r_a1, gamma_DM, rho0 = gp.files.params
    if gamma_DM == 0:
        nr = 3*rad/(r_DM+rad)
    elif gamma_DM == 1:
        nr = 2*rad/(r_DM+rad)+1
    else:
        print('unknown gamma_DM = ', gamma_DM)
        nr = 0.*rad
    return nr
## \fn nrgaiatot_3D(rad)
# plot - d log rho/ d log r for Gaia models
# @param rad radius in pc, not log.


def rhowalktot_3D(rad):
    rhodm, rhostar1, rhostar2 = rhowalk_3D(rad)     # 3* [msun/pc^3]
    return rhodm + rhostar1 + rhostar2                # [msun/pc^3]
## \fn rhowalktot_3D(rad)
# return total mass density, stars+DM
# @param rad radius in [pc]
# @return total 3D density in [Msun/pc^3]

def Mwalkertot(rbin):
    # based on underlying binned data
    # rhotot = rhowalktot_3D(rad)
    # Mtot = rho_SUM_Mr(rad, rhotot)

    # better method using quad numeric integration scheme with continuous rhowalktot
    from scipy.integrate import quad

    def igra(r):
        return 4.*np.pi*r**2*rhowalktot_3D(r)

    Mtot = np.zeros(gp.nipol); meps = np.zeros(gp.nipol)
    for i in range(gp.nipol):
        out=quad(igra,0,rbin[i],limit=100,full_output=0) # or np.inf
        Mtot[i] = out[0]
        meps[i] = out[1]
    # print('errors: ',meps)
    return Mtot
## \fn Mwalkertot(rbin)
# return total mass in Walker model
# @param rbin radii of bins
# @return total 3D mass for the Walker dataset in [Msun]


def q(i,j):
    from scipy.misc import factorial
    if i >= j and j >= 0:
        return (-1.)**j*factorial(i)/(factorial(j)*factorial(i-j))
    else:
        return 0.
## \fn q(i,j)
# start analytic generalized Hernquist profiles following equation (17) from Zhao 1996
# @param i first integer parameter
# @param j second integer parameter
# @return q quantity (see paper)


def apar(i,alpha,beta,gamma,rho0):
    C = rho0                    # TODO
    return 4.*np.pi*alpha*C/(alpha*(3-gamma)+i)*q(alpha*(beta-3.)-1.,i)
## \fn apar(i,alpha,beta,gamma,rho0)
# apar
# @param i
# @param alpha
# @param beta
# @gamma
# rho0



def bpar(i,alpha,beta,gamma,rho0):
    out = 0.
    if i>0:
        print('adjust range stepsize to +1 in analytic b()!')
        pdb.set_trace()
    for j in range(0,i-1,-1):
        out += q(alpha-1,j)*apar(i-j,alpha,beta,gamma,rho0)
    return out
## \fn bpar(i,alpha,beta,gamma,rho0)
# bpar
# @param i
# @param alpha
# @param beta
# @param gamma
# @param rho0


def S(i,chi):
    if i==0:
        return -np.log(chi)
    else:
        return (1.-chi**i)/i
## \fn S(i,chi)
# S
# @param i
# @param chi


    
def chi(r,alpha):
    ra = r**(1./alpha)
    return ra/(1.+ra)
## \fn chi(r,alpha)
# chi
# @param r radius in [pc]
# @param alpha


def Phi(r,alpha,beta,gamma,rho0):
    out=0.
    if alpha*(beta-2)-2>=0:
        print('adjust range stepsize to +1 in analytic Phi!')
        pdb.set_trace()
    for i in range(0,alpha*(beta-2)-2-1,-1):
        out += bpar(i,alpha,beta,gamma,rho0)*S(alpha*(2-gamma)+i,chi(r,alpha))
    return -out
## \fn Phi(r,alpha,beta,gamma,rho0)
# Phi
# @param r in [pc]
# @param alpha [1]
# @param beta [1]
# @param gamma [1]
# @param rho0 [Msun/pc^3]



def read_abc():
    A = np.loadtxt(gp.files.analytic, unpack=False)
    
    alpha_star1 = 2
    beta_star1  = int(A[8])
    gamma_star1 = int(A[7])
    r_star1     = 1000*A[9] #[pc] # both description of filename and file contents are wrong
    # name is rstar/10 in three digits, file content is given in kpc
    
    alpha_star2 = 2
    beta_star2  = int(A[12])
    gamma_star2 = int(A[11])
    r_star2     = 1000.*A[13] #[pc]

    alpha_DM    = int(A[18])
    beta_DM     = int(A[16])
    gamma_DM    = int(A[15])
    r_DM        = A[17] # [pc]

    rho0        = A[19] # [Msun/pc^3]
    return alpha_DM, beta_DM, gamma_DM, rho0, r_DM
## \fn read_abc()
# read analytic values


def Mabc_analytic(r):
    # read in rho0, alpha, beta, gamma
    alpha,beta,gamma,rho0,r_DM = read_abc()

    return -Phi(r,alpha,beta,gamma,rho0)*r/gp.G1    # [Msun]
## \fn Mabc_analytic(r)
# calculate enclosed mass
    
def PhiBeta(r,beta,C):
    return -C/r*(1.-(1.+r)**(3-beta))
## \fn PhiBeta(r,beta,C):
# equation in Zhang table p.2
# @param r in [pc]
# @param beta [1]
# @param C


def MBetaAnalytic(r):
    alpha,beta,gamma,rho0,r_DM = read_abc()
    r /= r_DM
    return -PhiBeta(r,beta,rho0)*r/gp.G1 # [Msun]
## \fn MBetaAnalytic(r)
# MBetaAnalytic
# @param r in [pc]

    
def Manalytic(r):
    alpha,beta,gamma,rho0,r_DM = read_abc()
    return 4*np.pi*rho0*(1./(1.+r/r_DM)+np.log(1.+r/r_DM)-1.)
## \fn Manalytic(r)
# analytic mass profile
# @param r radius in [pc]


def Rhowalk(rad):
    rhodm, rhostar1, rhostar2 = rhowalk(rad)     # 3* [msun/pc^3]

    surfdm    = rho_INT_Rho(rad, rhodm)    # [msun/pc^2]
    surfstar1 = rho_INT_Rho(rad, rhostar1) # [msun/pc^2]
    surfstar2 = rho_INT_Rho(rad, rhostar2) # [msun/pc^2]

    return surfdm, surfstar1, surfstar2               # 3* [msun/pc^2]
## \fn Rhowalk(rad)
# calculate 2D surface density from radius
# @param rad radius in [pc]
    
def Rhowalktot(rad):
    surfdm, surfstar1, surfstar2 = Rhowalk(rad)   # 3*[msun/pc^2]
    return surfdm + surfstar1 + surfstar2              # [msun/pc^2]
## \fn Rhowalktot(rad)
# return total surface density, stars+DM
# @param rad radius in [pc]

    
def betatriax(rad):
    eta = 0.5
    rsbeta = 810.                       # [pc]
    beta0 = 0.
    betainf = 0.5
    beta = (rsbeta**eta*beta0+rad**eta*betainf)/(rad**eta+rsbeta**eta)
    return beta
## \fn betatriax(rad)
# calculate velocity anisotropy for Dehnen-Wilkinson triaxial model
# @param rad in [pc]


def betawalker(rad):
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
## \fn betawalker(rad)
# calculate Osipkov-Merritt velocity anisotropy profile
# @param rad radius in [pc]


def betagaia(rad):
    beta_star1,r_DM,gamma_star1,r_star1,r_a1,gamma_DM,rho0 = gp.files.params
    return rad**2/(rad**2+r_a1**2)
## \fn betagaia(rad)
# Osipkov-Merritt velocity anisotropy profile
# @param rad radius in [pc]
