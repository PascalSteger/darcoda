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
import scipy.special as special
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

def sigz2_tilt(zvec, h_tracer, A, n): # Assume tilt tem to be of the form A*z**n
    F_vec = np.zeros(len(zvec))
    for i in range(len(zvec)):
        z = zvec[i]
        F_temp = h_tracer**(n+1.)*special.gamma(n+1.)*(special.gammainc(n+1.,z/h_tracer) - 1.)
        F_vec[i] = F_temp
    #print ('F_vec:',F_vec)
    #print ('zvec:',zvec)
    #print ('h_tracer:',h_tracer)
    #print ('(...):',special.gammainc(n+1.,z/h_tracer) - 1.)

    sig2_tilt = -A*F_vec*np.exp(zvec/h_tracer)
    return sig2_tilt

def sigz2_obsbary(zvec, params, h_tracer, rho_dm): # params: bary_params
    Sig_tot_inf =params[0] ; Sig_dwarf_inf =params[1] ; rho0_MS_thick =params[2]
    h = params[3] ; h1 = params[4]  ; h2 = params[5]  ; h3 = params[6]
    x = params[7] ; Sig_HII_inf = params[8] ; h_HII = params[9]
    
    G1 = 4.299e-6  # Newton's constant in (km)^2*kpc/(Msun*s^2)

    h_eff = (1.-x)*h2 + x*h3
    beta = (h-h1)/(h_eff-h1)

    int_dwarf_thin = Sig_dwarf_inf*(1.-beta)*h1/h * (h_tracer*0.61*h1/(h_tracer+0.61*h1)*np.exp(-zvec/(0.61*h1)) - h_tracer)
    #Sig_dwarf_thin = Sig_dwarf_inf*(1.-beta)*h1/h *(1.-np.exp(-zvec/(0.61*h1)))
    
    Sig_MS_thick_inf = 2.*h_eff*rho0_MS_thick
    Sig_thick_inf = Sig_dwarf_inf*beta*h_eff/h + Sig_MS_thick_inf
    
    temp1 = h_tracer*h2/(h_tracer+h2)*np.exp(-zvec/h2) - h_tracer
    #(1.-np.exp(-zvec/h2))
    temp2 = h_tracer*h3/(h_tracer+h3)*np.exp(-zvec/h3) - h_tracer
    #(1.-np.exp(-zvec/h3))
    int_Sig_thick = Sig_thick_inf/h_eff* ((1-x)*h2*temp1 + x*h3*temp2)
    #Sig_thick = Sig_thick_inf/h_eff*((1-x)*h2*(1.-np.exp(-zvec/h2)) + x*h3*(1.-np.exp(-zvec/h3)))
    
    int_Sig_HII = Sig_HII_inf*(h_tracer*h_HII/(h_tracer+h_HII)*np.exp(-zvec/h_HII) - h_tracer)
    #Sig_HII = Sig_HII_inf*(1.-np.exp(-zvec/h_HII))
    
    Sig_unimp = Sig_tot_inf - Sig_dwarf_inf - Sig_MS_thick_inf - Sig_HII_inf
    int_Sig_unimp = -Sig_unimp*h_tracer*np.ones(len(zvec))

    int_dm = -2.*rho_dm*h_tracer*(h_tracer+zvec)
    #Sig_dm = 2.* rho_dm * zvec
    
    sig2 = -2.*constants.pi*G1*(int_dwarf_thin + int_Sig_thick + int_Sig_HII + int_Sig_unimp + int_dm)
    #print ('---------------------------')
    #print ('int_dwarf_thin:',-int_dwarf_thin/1e6)
    #print ('int_Sig_thick',-int_Sig_thick/1e6)
    #print ('int_Sig_HII:',-int_Sig_HII/1e6)
    #print ('int_Sig_unimp:',-int_Sig_unimp/1e6)
    #print ('analytic sig2:',sig2)

    return sig2
# Analytic calculation of sigz2 for the obs_bary baryon model. Without tilt.
# Here sigz2 is positive for all z by definition.


# New function for calculating the z-dir velocity dispersion
# Constant C is now defined at zvec[-1] (instead of at zvec[0])
def sigz2(zvec,Sigvec,tiltvec,nuvec,C):
    #print ('In sigz2, C=',C)  # Print TAG
    G1 = 4.299e-6  # Newton's constant in (km)^2*kpc/(Msun*s^2)
    Kzvec = -2.*constants.pi*G1*Sigvec
    #print ('Kzvec:',Kzvec)
    #print ('tiltvec:',tiltvec)
    #print ('nuvec:',nuvec)
    #print ('integrand:',nuvec*(Kzvec-tiltvec))

    #old_integral = integrate.cumtrapz(nuvec*(Kzvec-tiltvec),zvec,initial= 0.) + C

    # Could maybe use integrate.quad 
    # Would then need a function for nuvec*(Kzvec-tiltvec):
    # phys.Sig_baryon_obs, tilt function, nu function, and DM func (with DD?)
    # Then need to do something like zvec -> np.array([zvec]) if zvec is float
    #  in called funcion.


    int_vec = nuvec*(Kzvec-tiltvec)

    integral = np.zeros(len(zvec))
    for i in range(len(zvec)):
        #zprim = zvec[i:]  ;  int_prim = int_vec[i:]
        intprim = integrate.cumtrapz(int_vec[i:],zvec[i:],initial = 0.)[-1]
        integral[i] = intprim

    #print ('gl_physics C:',C)
    #print ('integral:',integral)

    # C is defined at gp.z_bincenter_vecs[pop][-1] in gl_class_cube

    sig2 = (C - integral)/nuvec
    #print ('sigz2:',sig2)
    if (any(sig2<0.)):
        #print ('NEGATIVE VALUE IN sig2: min=',np.amin(sig2))
        gh.LOG(1,'Negative sig2 in phys.sigz')
        raise ValueError('negative value in sig2 array')
        return
    #if (sig2[-1] < sig2[-2]): # Disallow falling sig2 for last two bins
    #    #print ('sig2 falling at high z end')
    #    gh.LOG(1,'sig2 falling for high z in phys.sig2')
    #    raise ValueError('sig2 falling for high z')
    #    return

    return sig2    # Returns sig_z^2, sqrt taken in loglike if fitting to sig_z
# calculate z velocity dispersion using eq. 5 in 'almost' paper
# @param sigvec: velocity dispersion vector at the locations of zvec
# @param zvec: z-vector, at which nu/Sig(z) is given
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
    return (1/(4.*np.pi*G1)) * abs((K*(D**2)/((D**2 + zvec**2)**(1.5))))

def Sig_baryon_simplenu(zvec, params):  # SS: 6 May 2016
    G1 = 4.299e-6 # Newton's constant in (km)^2*kpc/(Msun*s^2)
    K=params[0]
    D=params[1]
    Sig_bary = (1/(2.*np.pi*G1)) * K*zvec/(np.sqrt(D**2 + zvec**2))
    return Sig_bary

def rho_baryon_obs(zvec, params):
    #z = 1000.*zvec  # converting so z in pc. No more! gl_params changed
    Sig_tot = params[0] ; Sig_dwarf_inf = params[1] ; rho0_MS_thick = params[2]
    h = params[3] ; h1 = params[4]  ; h2 = params[5]  ; h3 = params[6]
    x = params[7] ; Sig_HII_inf = params[8] ; h_HII = params[9]

    h_eff = (1.-x)*h2 + x*h3
    beta = (h-h1)/(h_eff-h1)

    rho0_dwarf = Sig_dwarf_inf/(2.*h)
    rho_dwarf_thin =rho0_dwarf*(1.-beta)/np.square(np.cosh(zvec/h1))

    rho0_thick = beta*rho0_dwarf + rho0_MS_thick
    rho_thick = rho0_thick*((1.-x)*np.exp(-zvec/h2) + x*np.exp(-zvec/h3))

    rho_HII = Sig_HII_inf/(2.*h_HII)*np.exp(-zvec/h_HII)

    rho_tot = rho_dwarf_thin + rho_thick + rho_HII

    if rho_tot[0]<rho_tot[-1]:
        print ('******************************')
        #print ('rho_dwarf_thin:',rho_dwarf_thin[0],' : ',rho_dwarf_thin[-1])
        print ('rho_thick:',rho_thick[0],' : ',rho_thick[-1])
        #print ('rho_HII:',rho_HII[0],' : ',rho_HII[-1])
        print ('h2:',h2,' h3:',h3,' x:',x)
        print ('1:',(1.-x)*h2*(np.ones(len(zvec))-np.exp(-zvec/h2)))
        print ('2:',x*h3*(np.ones(len(zvec))-np.exp(-zvec/h3)))
        pdb.set_trace()
    
    if np.amin([rho_dwarf_thin,rho_thick,rho_HII])<0:
        print ('rho_bary < 0 ERROR !!!!!!')
        pdb.set_trace()

    return rho_tot

def Sig_baryon_obs(zvec, params): # zvec in kpc
    #z = 1000.*zvec  # converting so z in pc. No more! gl_params changed
    Sig_tot_inf =params[0] ; Sig_dwarf_inf =params[1] ; rho0_MS_thick =params[2]
    h = params[3] ; h1 = params[4]  ; h2 = params[5]  ; h3 = params[6]
    x = params[7] ; Sig_HII_inf = params[8] ; h_HII = params[9]

    h_eff = (1.-x)*h2 + x*h3
    beta = (h-h1)/(h_eff-h1)

    #Sig_dwarf_thin1  = Sig_dwarf_inf*(1.-beta)*h1/h * np.tanh(zvec/h1)
    Sig_dwarf_thin = Sig_dwarf_inf*(1.-beta)*h1/h *(1.-np.exp(-zvec/(0.61*h1)))
    
    #print ('Sig ratio:',Sig_dwarf_thin1/Sig_dwarf_thin)

    Sig_MS_thick_inf = 2.*h_eff*rho0_MS_thick
    Sig_thick_inf = Sig_dwarf_inf*beta*h_eff/h + Sig_MS_thick_inf

    Sig_thick = Sig_thick_inf/h_eff*((1-x)*h2*(1.-np.exp(-zvec/h2)) + x*h3*(1.-np.exp(-zvec/h3)))

    Sig_HII = Sig_HII_inf*(1.-np.exp(-zvec/h_HII))

    Sig_unimp = Sig_tot_inf - Sig_dwarf_inf - Sig_MS_thick_inf - Sig_HII_inf

    if np.amin([Sig_dwarf_thin,Sig_thick,Sig_HII])<0 or Sig_unimp<0:
        print ('Sig_bary <0 ERROR !!!!!!!!!!!')
        print ('Sig_thick:',Sig_thick)
        pdb.set_trace()

    return Sig_unimp + Sig_dwarf_thin + Sig_thick + Sig_HII

def rho_baryon_obs_old(zvec, params):
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
    #G1 = 4.299e-6 # Newton's constant in (km)^2*kpc/(Msun*s^2)
    const_rho = params[0]
    #K=params[1]
    Sig_DD_inf=params[1]
    D=params[2]
    #DD_component =  1/(4*np.pi*G1) * abs((K*(D**2)/((D**2 + zvec**2)**(1.5))))
    DD_rho =  1./2.*Sig_DD_inf*D**2/(D**2 + zvec**2)**(1.5)
    #print ('const:',const,' K:',K,' D:',D)
    return const_rho + DD_rho

def Sig_dm_simplenu(zvec, params): # DM = const + dark disc
    const_Sig = 2.*params[0]*zvec
    Sig_DD_inf = params[1]
    D = params[2]
    DD_Sig = Sig_DD_inf*zvec/np.sqrt(np.square(zvec)+D**2)
    return const_Sig + DD_Sig


## \fn rho_baryon_simplenu(zvec, K, D)
# Calculate baryon density from the model used to generate the simplenu mock data
# @param zvec [kpc]
# @param K [kpc s^-2], baryon force parameter
# @param D [kpc], baryon disc scale height
