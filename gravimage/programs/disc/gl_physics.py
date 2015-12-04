#!/usr/bin/env ipython3

##
# @file
# calculations for velocity dispersion
# disc version

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import pdb
import numpy as np
import scipy.integrate as integrate
import scipy.constants as constants
import math
import gl_helper as gh

# New rho function:
def rho(zvec, kpar, rho0):
    logrho0 = math.log(rho0)
    logrho = np.array(logrho0 + integrate.cumtrapz(-kpar,zvec,initial=0.))
    rhovec = np.exp(logrho)
    return rhovec
# calculate density from k(z) = -d(ln(rho))/dz, rho at z=0, using trapezoidal rule
# @param kpar: kpar(z) vector, with kpar(0) as zeroth element
# @param rho0: density at z=0
# @param zvec: vector of z, at which k(z) is given (requirs z[0] != 0.)

# New function for calculating the surface density
def Sig(zvec, rhovec):
    Sigvec = 2.*np.array(integrate.cumtrapz(rhovec,zvec,initial=0.))
    return Sigvec
# calculate (total) surface density from the density, using trapezoidal rule
# @param Sigvec: surface density vector, Sigma[0] = 0.
# @param rhovec: rho(z) vector, with rho(0) as zeroth element.
# OBS: rho(z) is total density: rho = rho_DM + rho_bary (if return total Sig)
# @param zvec: z-vector, at which rho/Sig(z) is given (requires z[0] != 0.)

# New function for calculating the z-dir velocity dispersion
def sigz2(zvec,Sigvec,tiltvec,nuvec,C):
    G1 = 4.299e-6  # Newton's constant in (km)^2*kpc/(Msun*s^2)
    Kzvec = -2.*constants.pi*G1*Sigvec
    integral = integrate.cumtrapz(nuvec*(Kzvec-tiltvec),zvec,initial=0.) + C
    sig2 = integral/nuvec
    if (any(sig2<0.)):
        gh.LOG(1,'Negative sig2 in phys.sigz')
        raise ValueError('negative value in sig2 array')
        return
    #sigvec = np.sqrt(sig2)
    return sig2
# calculate z velocity dispersion using eq. 5 in 'almost' paper
# @param sigvec: velocity dispersion vector at the locations of zvec
# @param zvec: z-vector, at which nu/Sig(z) is given (assuming z[0] = 0.)
# @param nuvec: tracer number density at the locations of zvec
# @param C: Integration constant C in eq (5) in the 'almost' paper
# all arrays (zvec,Sigvec,nuvec) are required to be numpy arrays.
# outputs sigvec at z = zvec[1:], eg discards first z point


def tilt(zipol, param, gp):
    tilttmp = 0.
    for i in range(gp.nbeta):
        tilttmp += param[i]*(zipol/max(zipol))**i
    return tilttmp
## \fn tilt(zipol, param, gp)
# return sum of polynomials for tilt as fct of radius
# TODO: get tilt size from Garbari+2011
# @param zipol [pc]
# @param param n_beta parameters
# @param gp global parameters


def kappa(xipol, Kz):
    z0 = xipol
    kzpar = -np.hstack([Kz[0]/z0[0], (Kz[1:]-Kz[:-1])/(z0[1:]-z0[:-1])])
    return kzpar
## \fn kappa(xipol, Kz)
# calculate vertical force from Kz
# @param xipol height above midplane [pc]
# @param Kz acceleration in z direction


def sig_rz(z, zpars, tpars):
    # Mirror prior
    tparsu = abs(tpars)

    # dz = zpars[2]-zpars[1]
    # sig_Rz = np.zeros(gp.nipol)
    # sig_Rz[0] = tparsu[0] * dz / 2.
    # for i in range(1,gp.nipol):
    #   sig_Rz[i] = sig_Rz[i-1] + tparsu[i] * dz
    # f = gp.ipol(zpars,sig_Rz,z)

    # Alternative here --> don't assume monotonic!
    f = gh.ipol(zpars, tparsu, z)
    return f
## \fn sig_rz(z, zpars, tpars)
# General function to describe the tilt profile
# @param z [pc]
# @param zpars [pc] z, on which sig is defined
# @param tpars tilt parameters: Rsun, hr, hsig


def rho_baryon_simplenu(zvec, params):
    G1 = 4.299e-6 # Newton's constant in (km)^2*kpc/(Msun*s^2)
    K=params[0]
    D=params[1]
    return (1/(4*np.pi*G1)) * abs((K*(D**2)/((D**2 + zvec**2)**(1.5))))

def rho_baryon_obs(zvec, params):
    H2_Sigma = params[0]
    H2_h = params[1]
    HIcnm_Sigma = params[2]
    HIcnm_h = params[3]
    HIwnm1_Sigma = params[4]
    HIwnm1_h = params[5]
    HIwnm2_Sigma = params[6]
    HIwnm2_h = params[7]
    HII_Sigma = params[8]
    HII_h = params[9]
    MS3_Sigma = params[10]
    MS3_h = params[11]
    MS4_Sigma = params[12]
    MS4_h = params[13]
    MS5_Sigma = params[14]
    MS5_h = params[15]
    MS8_Sigma = params[16]
    MS8_h = params[17]
    MS_thick_fraction = params[18]
    MS_thick_h = params[19]
    dwarfs_Sigma = params[20]
    dwarfs_h1 = params[21]
    dwarfs_h2 = params[22]
    dwarfs_beta = params[23]

    #This calculation is performed in units of Msun and pc, while the rest of
    # the code uses kpc. Hence convert zvec to pc:
    zvec = 1000.*zvec

    H2_dens = H2_Sigma/(H2_h*math.sqrt(math.pi))*np.exp(-np.square(zvec/H2_h))
    HIcnm_dens = HIcnm_Sigma/(HIcnm_h*math.sqrt(math.pi))*np.exp(-np.square(zvec/HIcnm_h))
    HIwnm1_dens = HIwnm1_Sigma/(HIwnm1_h*math.sqrt(math.pi))*np.exp(-np.square(zvec/HIwnm1_h))
    HIwnm2_dens = HIwnm2_Sigma/(2.*HIwnm2_h)*np.exp(-zvec/HIwnm2_h)
    HII_dens = HII_Sigma/(2.*HII_h)*np.exp(-zvec/HII_h)

    MS3_dens = MS3_Sigma/(2.*MS3_h*np.square(np.cosh(zvec/MS3_h)))
    MS4_dens = MS4_Sigma/(2.*MS4_h*np.square(np.cosh(zvec/MS4_h)))
    MS5_dens = MS5_Sigma/(2.*MS5_h*np.square(np.cosh(zvec/MS3_h)))*(1.-MS_thick_fraction)
    MS8_dens = MS8_Sigma/(2.*MS8_h*np.square(np.cosh(zvec/MS8_h)))*(1.-MS_thick_fraction)
    MS_thick_dens = MS_thick_fraction*(MS5_Sigma+MS8_Sigma)/(2.*MS_thick_h)*np.exp(-zvec/MS_thick_h)

    h_dwarf = (1.-dwarfs_beta)*dwarfs_h1 + dwarfs_beta*dwarfs_h2
    dwarf_dens = dwarfs_Sigma/(2.*h_dwarf)*((1.-dwarfs_beta)*(np.cosh(zvec/dwarfs_h1))**-2 + dwarfs_beta*np.exp(-zvec/dwarfs_h2))

    star_dens = MS3_dens+MS4_dens+MS5_dens+MS8_dens+MS_thick_dens+dwarf_dens
    #star_dens = dwarf_dens
    gas_dens =  H2_dens+HIcnm_dens+HIwnm1_dens+HIwnm2_dens+HII_dens

    # Above is in units of Msun/pc3, transforming into units Msun/kpc3 below
    gas_dens = gas_dens*1e9     #Converting to Msun/kpc3
    star_dens = star_dens*1e9   # Converting to Msun/kpc3

    return gas_dens+star_dens
    #return (0.5*star_dens + 5.*gas_dens)/2. #BODGE to make the baryon dens more similar to simplenu
    # TODO FIXME !!! 
    # Divided by 2 to see if less baryons improve multinest fit.

def rho_dm_simplenu(zvec, params): # DM = const + dark disc
    G1 = 4.299e-6 # Newton's constant in (km)^2*kpc/(Msun*s^2)
    const = params[0]
    K=params[1]
    D=params[2]
    DD_component =  1/(4*np.pi*G1) * abs((K*(D**2)/((D**2 + zvec**2)**(1.5))))
    #print ('const:',const,' K:',K,' D:',D)
    return (const + DD_component)


## \fn rho_baryon_simplenu(zvec, K, D)
# Calculate baryon density from the model used to generate the simplenu mock data
# @param zvec [kpc]
# @param K [kpc s^-2], baryon force parameter
# @param D [kpc], baryon disc scale height
