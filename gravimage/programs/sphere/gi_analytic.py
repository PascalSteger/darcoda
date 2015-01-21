#!/usr/bin/env ipython3

##
# @file
# @ingroup gravimage
# all analytic profiles from gi_physics

# (c) GPL v3 2014 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import pdb

import gi_units as gu
import gi_project as gip
import gi_helper as gh

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

def rho_general(r0, rscale, rho0, alpha, beta, gamma):
    # rho_DM(r) = rho_0 (r/r_DM)^(-gamma_DM) *
    #             (1+(r/r_DM)^alpha_DM)^((gamma_DM-beta_DM)/alpha_DM)
    tmp = 1.*rho0                       # [Munit/pc^3]
    tmp = tmp * (r0/rscale)**(-gamma)   # [Munit/pc^3]
    tmp = tmp * (1.+(r0/rscale)**alpha)**((gamma-beta)/alpha) # [Munit/pc^3]
    return np.array(tmp)                          # [Munit/pc^3]
## \fn rho_general(r0, rscale, rho0, alpha, beta, gamma):
# determine rho(r) for a generalized Hernquist model
# @param r0 radius in [pc]
# @param rscale scale radius in [pc]
# @param rho0 central 3D mass density, in [Munit/pc^3]
# @param gamma inner slope
# @param alpha outer slope
# @param beta steepness of turnover
# @return 3D Hernquist density in [Munit/pc^3]

def rho_triax(rad, gp):
    alpha = 1.
    beta = 4.
    rs = 1500.                          # [pc]
    if gp.case <= 4:               # core
        gamma = 0.23
        rhos  = 1.177E-1                # [Munit/pc^3]
    elif gp.case > 4 and gp.case <=8:             # cusp
        gamma = 1.
        rhos  = 5.522E-2                # [Munit/pc^3]
    else:
        raise Exception('wrong case for triax system')

    rho = rhos
    rho /= (rad/rs)**gamma
    rho /= (1+(rad/rs)**(1/alpha))**(alpha*(beta-gamma))
    return rho
## \fn rho_triax(rad, gp)
# density for a triaxial halo in Gaia challenge
# @param rad radius in [pc]
# @param gp global parameters
# @return density in [Munit/pc^3]

def rho_walk(rad, gp, mf1=1, mf2=1):
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
    rho0        = A[19] # [Munit/pc^3]

    # attention: the stellar tracers in Walker's datasets do not have mass density profile!
    rho0star1    = rho0*mf1

    # rho(r, rscale, rho0, alpha, beta, gamma):
    #  (2*[pc], or 2*[rcore]), [Munit/pc^3], 3*[1]
    rhodm    = rho_general(rad, r_DM, rho0, alpha_DM, beta_DM, gamma_DM) # [msun/pc^3]
    rhostar1 = rho_general(rad, r_star1, rho0star1, alpha_star1, beta_star1, gamma_star1)

    if gp.pops==1:
        return rhodm, rhostar1
    # [msun/pc^3]

    #ntracer2  = gp.ntracer[2]
    rho0star2 = rho0*mf2#/1.e6*ntracer2
    rhostar2  = rho_general(rad, r_star2, rho0star2, alpha_star2, beta_star2, gamma_star2)

    return rhodm, rhostar1, rhostar2                 # 3* [Munit/pc^3]
## \fn rho_walk(rad, gp, mf1, mf2)
# Walker model: read values from theoretical params file,
# calculate from eq. 2 generalized Hernquist profiles
# @param rad radius in [pc]
# @param gp global parameters
# @param mf1 factor to go from rho to rho_star1
# @param mf2 factor to go from rho to rho_star2
# @return 3D density in [Munit/pc^3], for each DM, stellar pop 1, stellar pop 2 (if available)


def rho_gaia(rad, gp):
    if gp.investigate != 'gaia':
        raise Exception('wrong investigation!')
    alpha_star1 = 2.
    alpha_DM = 1.
    beta_DM = 3.
    if gp.case == 9 or gp.case == 10:
        alpha_star1 = 0.5
        beta_DM = 4.

    beta_star1, r_DM, gamma_star1, r_star1, r_a1,\
        gamma_DM,rho0 = gp.files.params

    if gamma_star1 == 0.1:
        nu0 = 2.2e7/r_star1**3
    elif gamma_star1 == 1.0:
        nu0 = 1.5e7/r_star1**3
    gh.LOG(2, '  analytic rho_gaia:')
    gh.LOG(2, '   rho0 = ',rho0)
    gh.LOG(2, '   r_DM = ', r_DM)
    gh.LOG(2, '   r_star1 = ', r_star1)

    rhodm = rho_general(rad, r_DM, rho0, alpha_DM, beta_DM, gamma_DM)
    rhostar1 = rho_general(rad, r_star1, nu0, \
                   alpha_star1, beta_star1, gamma_star1)
    return rhodm, rhostar1
## \fn rho_gaia(rad, gp)
# give densities for the spherical Gaia challenge dataset
# @param rad radius in [pc]
# @param gp global parameters
# @return 3D densities (DM+stellar population) in [Munit/pc^3]


def M_gaia(rad, gp):
    rhodm, rhostar = rho_gaia(rad, gp)
    # 3D radius here

    beta_star1, r_DM, gamma_star1, r_star1, r_a1, gamma_DM, rho0 = gp.files.params
    if gamma_DM == 1:
        s = rad/r_DM
        Mstar = 4.*np.pi*rho0*r_DM**3*(1/(1+s)+np.log(1+s)-1) # [Msun]
    else:
        #Mdm = gip.rho_INT_Sum_MR(rad, rhodm, gp)
        Mstar = gip.rho_INT_Sum_MR(rad, rhostar, gp)
    # from here on, we assume to work on 2D radii
    return Mstar, Mstar
## \fn M_gaia(rad, gp)
# calculate mass in 2D radial bins
# @param rad radius in pc, 3D
# @param gp global parameters


def rhotot_gaia(rad, gp):
    rhodm, rhostar1 = rho_gaia(rad, gp)
    return rhodm
## \fn rhotot_gaia(rad, gp)
# return total mass density for Gaia challenge models
# @param rad radii in [pc]
# @param gp global parameters


def Sig_gaia(rad, gp):
    rhodm, rhostar1 = rho_gaia(rad, gp)
    Sigdm = gip.rho_INT_Sig(rad, rhodm, gp)
    Sigstar = gip.rho_INT_Sig(rad, rhostar1, gp)
    return Sigdm, Sigstar
## \fn Sig_gaia(rad, gp)
# get projected surface density for Gaia tracer population
# not based on analytic integral
# @param rad radii in [pc]
# @param gp global parameters


def nu3Dtot_gaia(rad, gp):
    if gp.investigate != 'gaia':
        raise Exception('wrong investigation!')
    alpha_star1 = 2.
    #alpha_DM = 1.
    beta_DM = 3.
    if gp.case == 9 or gp.case == 10:
        alpha_star1 = 0.5
        #beta_DM = 4.
    beta_star1, r_DM, gamma_star1, r_star1, r_a1, gamma_DM, rho0 = gp.files.params
    nustar1 = rho_general(rad, r_star1, rho0, alpha_star1, beta_star1, gamma_star1)
    return nustar1
## \fn nu3Dtot_gaia(rad, gp)
# give total density for the spherical Gaia challenge dataset
# @param rad radius in [pc]
# @param gp global parameters
# @return 3D overall density (DM+stellar population) in [Munit/pc^3]


def nr3Dtot_deriv_walk(rad, gp):
    # TODO too high for walk1, core
    lrho = np.log(rhotot_walk(rad, gp))
    lr   = np.log(rad)
    import gi_helper as gh
    return -gh.derivcoarse(lrho, lr)
## \fn nr3Dtot_deriv_walk(rad, gp)
# plot d log rho/d log r
# @param rad radius in pc, not in log
# @param gp global parameters


def nr3Dtot_deriv_triax(rad, gp):
    lrho = np.log(rho_triax(rad, gp))
    lr   = np.log(rad)
    import gi_helper as gh
    return -gh.derivcoarse(lrho, lr)
## \fn nr3Dtot_deriv_triax(rad, gp)
# plot d log rho/d log r
# @param rad radius in pc, not in log
# @param gp global parameters

def nr3Dtot_deriv_gaia(rad):
    lrho = np.log(rhotot_gaia(rad))
    lr   = np.log(rad)
    import gi_helper as gh
    return -gh.derivcoarse(lrho, lr)
## \fn nr3Dtot_deriv_gaia(rad)
# plot d log rho/d log r
# @param rad radius in pc, not in log


def nr3Dtot_walk(rad):
    beta_star1, r_DM, gamma_star1, r_star1, r_a1, \
      gamma_DM, rho0 = gp.files.params
    if gamma_DM == 0:
        nr = 3*rad/(r_DM+rad)
    elif gamma_DM == 1:
        nr = 2*rad/(r_DM+rad)+1
    else:
        gh.LOG(1, 'unknown gamma_DM = ', gamma_DM)
        nr = 0.*rad
    return nr
## \fn nr3Dtot_walk(rad)
# plot - d log rho/ d log r for Walker models
# @param rad radius in pc, not log.


def nr3Dtot_gaia(rad, gp):
    beta_star1, r_DM, gamma_star1, r_star1, r_a1, gamma_DM, rho0 = gp.files.params
    if gamma_DM == 0:
        nr = 3*rad/(r_DM+rad)
    elif gamma_DM == 1:
        nr = 2*rad/(r_DM+rad)+1
    else:
        gh.LOG(1, 'unknown gamma_DM = ', gamma_DM)
        nr = 0.*rad
    return nr
## \fn nr3Dtot_gaia(rad, gp)
# plot - d log rho/ d log r for Gaia models
# @param rad radius in pc, not log.
# @param gp global parameters

def rhotot_walk(rad, gp, mf1=1, mf2=1):
    if gp.pops==1:
        rhodm, rhostar1 = rho_walk(rad, gp, mf1, mf2)         # 3* [msun/pc^3]
    elif gp.pops==2:
        rhodm, rhostar1, rhostar2 = rho_walk(rad, gp, mf1, mf2)
    return rhodm               # [msun/pc^3]
## \fn rhotot_walk(rad, gp, mf1, mf2)
# return total mass density, stars+DM
# @param rad radius in [pc]
# @param mf1 factor to go from DM rho0 to rho0_star1
# @param mf2 factor to go from DM rho0 to rho0_star2
# @param gp global parameters
# @return total 3D density in [Munit/pc^3]. in case of Walker dataset, we have that the tracer particles do not contribute to the overall density, thus we use the DM density only

def Mtot_walk(rbin, gp, mf1 = 1, mf2 = 1):
    # based on underlying binned data
    # rhotot = rhotot_walk(rad, gp, mf1, mf2)
    # Mtot = gip.rho_SUM_Mr(rad, rhotot)
    # better method using quad numeric integration scheme with continuous rhotot_walk
    from scipy.integrate import quad
    def igra(r):
        return 4.*np.pi*r**2*rhotot_walk(r, gp, mf1, mf2)
    Mtot = np.zeros(gp.nipol)
    meps = np.zeros(gp.nipol)
    for i in range(gp.nipol):
        out = quad(igra, 0, rbin[i], limit=100, full_output=0) #or np.inf
        Mtot[i] = out[0]
        meps[i] = out[1]
    # gh.LOG(1, 'errors: ',meps)
    return Mtot
## \fn Mtot_walk(rbin, gp, mf1, mf2)
# return total mass in Walker model
# NOT USED ANYMORE
# @param rbin radii of bins
# @param gp global parameters
# @param mf1 factor rho0 pop 1
# @param mf2 factor rho0 pop 2
# @return total 3D mass for the Walker dataset in [Munit]

def Sig_walk(rad, gp, mf1=1, mf2=1):
    rhodm, rhostar1, rhostar2 = rho_walk(rad, gp, mf1, mf2)     # 3* [msun/pc^3]
    Sig_star1 = gip.rho_INT_Sig(rad, rhostar1, gp) # [msun/pc^2]
    Sig_star2 = gip.rho_INT_Sig(rad, rhostar2, gp) # [msun/pc^2]
    return Sig_star1+Sig_star2, Sig_star1, Sig_star2               # 3* [msun/pc^2]
## \fn Sig_walk(rad, gp, mf1, mf2)
# calculate 2D surface density from radius for pop 0, 1, 2
# @param rad radius in [pc]
# @param gp global parameters
# @param mf1 factor rho0 pop 1
# @param mf2 factor rho0 pop 2

def beta_walk(rad, gp):
    # Osipkov-Merritt anisotropy profile with r_a/r_* = 10^4 for isotropic models
    A = np.loadtxt(gp.files.analytic, unpack=False)
    import pdb
    rars1  = A[10]
    rs1    = A[9] * 1000. # [pc]
    rbyrs1 = rad/rs1
    beta1  = rbyrs1**2 / (rars1**2+rbyrs1**2)
    rars2  = A[14]
    rs2    = A[13] * 1000. # [pc]
    rbyrs2 = rad/rs2
    beta2  = rbyrs2**2 / (rars2**2+rbyrs2**2)
    return -1, beta1, beta2
## \fn beta_walk(rad, gp)
# calculate Osipkov-Merritt velocity anisotropy profile, for pop 0 (dummy), 1, 2
# @param rad radius in [pc]
# @param gp

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

def apar(i, alpha, beta, gamma, rho0):
    C = rho0
    return 4.*np.pi*alpha*C/(alpha*(3-gamma)+i)*q(alpha*(beta-3.)-1.,i)
## \fn apar(i,alpha,beta,gamma,rho0)
# a_parameter
# @param i int
# @param alpha see Zhao 1996
# @param beta see Zhao 1996
# @param gamma see Zhao 1996
# @param rho0 central density

def bpar(i, alpha, beta, gamma, rho0):
    out = 0.
    if i > 0:
        raise Exception('adjust range stepsize to +1 in analytic b()!')
    for j in range(0,i-1,-1):
        out += q(alpha-1,j)*apar(i-j,alpha,beta,gamma,rho0)
    return out
## \fn bpar(i,alpha,beta,gamma,rho0)
# bpar
# @param i see Zhao 1996
# @param alpha see Zhao 1996
# @param beta see Zhao 1996
# @param gamma see Zhao 1996
# @param rho0 central density

def S(i,chi):
    if i==0:
        return -np.log(chi)
    else:
        return (1.-chi**i)/i
## \fn S(i,chi)
# S function
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
        raise Exception('adjust range stepsize to +1 in analytic Phi!')
    for i in range(0,alpha*(beta-2)-2-1,-1):
        out += bpar(i,alpha,beta,gamma,rho0)*S(alpha*(2-gamma)+i,chi(r,alpha))
    return -out
## \fn Phi(r,alpha,beta,gamma,rho0)
# Phi
# @param r in [pc]
# @param alpha [1]
# @param beta [1]
# @param gamma [1]
# @param rho0 [Munit/pc^3]

def read_abc(gp):
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
    rho0        = A[19] # [Munit/pc^3]
    return alpha_DM, beta_DM, gamma_DM, rho0, r_DM
## \fn read_abc(gp)
# read analytic values for Walker mock data
# @param gp global parameters

def Mabc_analytic(r, gp):
    # read in rho0, alpha, beta, gamma
    alpha,beta,gamma,rho0,r_DM = read_abc(gp)
    return -Phi(r,alpha,beta,gamma,rho0)*r/gu.G1__pcMsun_1km2s_2    # [Munit]
## \fn Mabc_analytic(r)
# calculate enclosed mass
# @param gp global parameters

def PhiBeta(r,beta,C):
    return -C/r*(1.-(1.+r)**(3-beta))
## \fn PhiBeta(r,beta,C):
# equation in Zhang table p.2
# @param r in [pc]
# @param beta [1]
# @param C

def MBetaAnalytic(r, gp):
    alpha,beta,gamma,rho0,r_DM = read_abc(gp)
    r /= r_DM
    return -PhiBeta(r,beta,rho0)*r/gu.G1__pcMsun_1km2s_2 # [Munit]
## \fn MBetaAnalytic(r)
# MBetaAnalytic
# @param r in [pc]
# @param gp global parameters

def Manalytic(r, gp):
    alpha,beta,gamma,rho0,r_DM = read_abc(gp)
    return 4*np.pi*rho0*(1./(1.+r/r_DM)+np.log(1.+r/r_DM)-1.)
## \fn Manalytic(r)
# analytic mass profile
# @param r radius in [pc]
# @param gp global parameters

def beta_gaia(rad, gp):
    if gp.investigate != 'gaia':
        raise Exception('wrong investigation!')
    if gp.case == 9 or gp.case == 10:
        # constant tangential velocity anisotropy of beta=-0.5
        return rad*0.0-0.5, rad*0.0-0.5
    beta_star1, r_DM, gamma_star1, r_star1, r_a1, gamma_DM, rho0 = gp.files.params
    beta = rad**2/(rad**2+r_a1**2)
    return beta, beta
## \fn beta_gaia(rad, gp)
# Osipkov-Merritt velocity anisotropy profile, for pop 0, pop 1
# @param rad radius in [pc]
# @param gp global parameters

def beta_hern(rad):
    return 0.*rad, 0.*rad
## \fn beta_hern(rad)
# analytic value for beta in Hernquist case
# @param rad radius [pc]

def beta_triax(rad):
    eta = 0.5
    rsbeta = 810.                       # [pc]
    beta0 = 0.
    betainf = 0.5
    beta = (rsbeta**eta*beta0+rad**eta*betainf)/(rad**eta+rsbeta**eta)
    return beta
## \fn beta_triax(rad)
# calculate velocity anisotropy for Dehnen-Wilkinson triaxial model
# @param rad in [pc]

def beta(rad, gp):
    if gp.investigate =='hern':
        return beta_hern(rad)
    elif gp.investigate =='gaia':
        return beta_gaia(rad, gp)
    elif gp.investigate == 'walk':
        return beta_walk(rad, gp)
    elif gp.investigate == 'triax':
        return beta_triax(rad)
    else:
        gh.LOG(1, 'ga.beta not defined')
        pdb.set_trace()
## \fn beta(rad, gp)
# analytic profile for all investigations
# @param rad radii in pc
# @param gp global parameters

def rho_hern(r0, gp):
    s = r0/gp.ana                              # [1]
    rho = gp.anM/(2.*np.pi*gp.ana**3)/(s*(1+s)**3) # [Munit/pc^3]
    return rho, rho
## \fn rho_hern(r0, gp)
# equation 2 from Hernquist 1990, 3D mass density
# equation 2b from Baes&Dejonghe 2002
# @param r0 radius in [pc]
# @param gp global parameters
# @return density in [Munit/pc^3]

def rho_hern2(r0, gp, gamma = 1):
    rho0 = M*(3.-gamma)/(4.*np.pi*gp.ana**3)
    rho0 *= (r0/gp.ana)**(-gamma)
    rho0 *= (1.+r0/gp.ana)**(gamma-4.)
    return rho0, rho0
## \fn rho_hern2(r0, gp, gamma = 1)
# find theoretical mass profile of double-Dehnen spheres with stars,
# according to http://www.astrosim.net/code/doku.php?id=home:codetest:massmodel:sphericaltest
# @param r0 radius in pc
# @param gp global parameters
# @param gamma inner density slope. default: cusp

def rho(r0, gp):
    if gp.investigate == 'hern':
        return rho_hern(r0, gp)
    elif gp.investigate == 'gaia':
        return rho_gaia(r0, gp)
    elif gp.investigate == 'walk':
        return rho_walk(r0, gp)
    else:
        gh.LOG(1, 'ga.rho not defined')
        pdb.set_trace()
## \fn rho(rad, gp)
# analytic total mass density profile for all investigations
# @param rad radii in pc
# @param gp global parameters

def M_hern(r0, gp):
    s = r0/gp.ana                            # [1]
    Mr =  gp.anM*s**2/(1.+s)**2             # [Munit]
    return Mr, Mr
## \fn M_hern(r0, gp)
# equation 3 from Hernquist 1990
# @param r0 radius in [pc]
# @param gp global parameters
# @return 3D mass in [Munit]

def Mr(r0, gp):
    if gp.investigate == 'hern':
        return M_hern(r0, gp)
    elif gp.investigate == 'gaia':
        return M_gaia(r0, gp)
    else:
        gh.LOG(1, 'ga.Mr not defined')
        pdb.set_trace()
## \fn Mr(rad, gp)
# analytic total mass profile for all investigations
# @param rad radii in pc
# @param gp global parameters

def Sig_hern(r0, gp):
    s = np.array(r0)/gp.ana                                          # [1]
    Sigma =  gp.anM/(2.*np.pi*gp.ana**2)*((2.+s**2)*X(s)-3.)/((1.-s**2)**2) # [Munit/pc^2]
    return Sigma, Sigma
## \fn Sig_hern(r0, gp)
# surface density of Hernquist profile, equation 3 from Baes&Dejonghe 2002
# @param r0 radius in [pc]
# @param gp global parameters
# @return 2D surface density in [Munit/pc^2]

def Sigma(r0, gp):
    if gp.investigate == 'hern':
        return Sig_hern(r0, gp)
    elif gp.investigate == 'gaia':
        return Sig_gaia(r0, gp)
    elif gp.investigate == 'walk':
        return Sig_walk(r0, gp)
    else:
        gh.LOG(1, 'ga.Sigma not defined')
        pdb.set_trace()
## \fn Sigma(rad, gp)
# analytic total mass surface density profile for all investigations
# @param rad radii in pc
# @param gp global parameters
# @return Sigma_0, Sigma_1

def sigr2_hern(r0, gp):
    s = np.array(r0)/gp.ana          # [pc/1000pc]
    G1 = 1
    return G1*gp.anM/gp.ana*s*(1+s)**3*np.log(1.+1./s)-\
           G1*gp.anM/(12.*gp.ana)*s/(s+1)*(25.+52.*s+42.*s**2+12.*s**3) # [(km/s)^2]
## \fn sigr2_hern(r0, gp)
# sig_r^2, equation 10 of Hernquist 1990
# @param r0 radius in [pc]
# @param gp global parameters
# @return sig_r^2 in [(km/s)^2]

def sigr2(r0, gp):
    if gp.investigate == 'hern':
        return sigr2_hern(r0, gp)
    else:
        gh.LOG(1, 'ga.sigr2 not defined')
        return 0*r0-1
## \fn sigr2(rad, gp)
# analytic total mass surface density profile for all investigations
# @param rad radii in pc
# @param gp global parameters
# @return sigr2

def sig_los_hern(r0, gp):
    return np.sqrt(Sig_sig_los_2_hern(r0, gp)/Sig_hern(r0, gp)[1])
## \fn sig_los_hern(r0, gp)
# sig_los determined from analytic surfden*sig2 and surfden
# @param r0 radius in [pc]
# @param gp global parameters
# @return sig_LOS for pop  in [km/s]

def sig_los(r0, gp):
    if gp.investigate == 'hern':
        return sig_los_hern(r0, gp)
    else:
        gh.LOG(1, 'ga.sig_los not defined')
        return 0*r0-1
## \fn sig_los(r0, gp)
# sigma_LOS for all investigations
# @param r0 radii in pc
# @param gp global parameters

def kappa_hern(r0):
    return 3.*np.ones(len(r0))          # [1]
## \fn kappa_hern(r0)
# analyitc value for the fourth order momentum of the LOS velocity
# @param r0 radius in [pc]
# @return 3 for Gaussian velocity distribution

def Sig_sig_los_2_hern(r0, gp):
    # \sigma_p = \sigma_projected = \sigma_{LOS}
    if gp.investigate != 'hern':
        gh.LOG(1, 'wrong investigation')
        pdb.set_trace()
    s = r0/gp.ana                            # [1]
    G1 = 1
    return G1*gp.anM**2/(12.*np.pi*gp.ana**3)*(1./(2.*(1.-s**2)**3)\
                                          *(-3.*s**2*X(s)\
                                          *(8.*s**6-28.*s**4+35.*s**2-20.)\
                                          -24.*s**6+68.*s**4-65.*s**2+6.)\
                                          -6.*np.pi*s) # [(km/s)^2 * Munit/pc^2]
## \fn Sig_sig_los_2_hern(r0, gp)
# equation 21 from Baes&Dejonghe 2002
# @param r0 radius in [pc]
# @param gp global parameters
# @return surface density * sig_LOS^2 in [(km/s)^2 * Munit/pc^2]

def Sig_sig_los_2(r0, gp):
    if gp.investigate == 'hern':
        return Sig_sig_los_2_hern(r0, gp)
    else:
        gh.LOG(1, 'ga.Sig_sig_los_2 not defined')
        return 0.*r0-1
## \fn Sig_sig_los_2(rad, gp)
# surface density time sigma_r^2
# @param rad radii in pc
# @param gp global parameters
# @return Sigma*sigma_r^2
