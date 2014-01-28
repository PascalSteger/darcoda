#!/bin/env python3

##
# @file
# global params for data analysis step 1 used to generate input for spherical MCMC

# (c) 2013 ETHZ Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import pdb
import gl_params as gp
from gl_class_files import *

showplots = False
n = 200
nbins = gp.nipol


if gp.investigate == 'hern':
    repr  = 1     # choose simulation representation
    ncomp = gp.pops+1     # number of populations + 1 (for all comps together)
    # numbers of tracers. both 0: take all particles
    ntracers1 = 10000  # stars
    ntracers2 = 0      # DM      (in dual_unit_hern_*)
    
    munit = 1.                                 # [Msun]
    Rcut = 1.e10                               # [Rvir]
    
    Rmin = 0. * gp.ascale; Rmax = 3.*gp.ascale # [pc]
    rprior = 3.
    
    nit = 30             # set number of iterations
    procs = 1            # number of processes for parallel execution
    
    dir = gp.files.dir
        
    simname = gp.files.get_sim_name() # dir+'simulation/'+prename+'unit_hern_%i' %(repr)
    simpos = dir+'simulation/'+simname+'stars_pos.txt'
    simvel = dir+'simulation/'+simname+'stars_vel.txt'
    fileposcartesian = []; fileposcenter = []; filevelcartesian = [];
    fileposcartesian.append(dir+'simulation/'+simname+'pos_0.txt')
    fileposcartesian.append(dir+'simulation/'+simname+'pos_1.txt')
    fileposcenter.append(dir+'simulation/'+simname+'centeredpos_0.txt')
    fileposcenter.append(dir+'simulation/'+simname+'centeredpos_1.txt')
    filevelcartesian.append(dir+'simulation/'+simname+'vel_0.txt')
    filevelcartesian.append(dir+'simulation/'+simname+'vel_1.txt')
    
    fileposspherical =[]; filevelspherical = []
    fileposspherical.append(dir+'simulation/'+simname+'sphericalpos_0.txt')
    fileposspherical.append(dir+'simulation/'+simname+'sphericalpos_1.txt')
    filevelspherical.append(dir+'simulation/'+simname+'sphericalvel_0.txt')
    filevelspherical.append(dir+'simulation/'+simname+'sphericalvel_1.txt')
        
    filemass = []; filedenfalloff = []; filesig = []; filekappa = []
    filemass.append(dir+'enclosedmass/'+simname+'enclosedmass_0.txt')
    filemass.append(dir+'enclosedmass/'+simname+'enclosedmass_1.txt')
    filedenfalloff.append(dir+'densityfalloff/'+simname+'falloffnotnorm_0.txt')
    filedenfalloff.append(dir+'densityfalloff/'+simname+'falloffnotnorm_1.txt')
    filesig.append(dir+'siglos/'+simname+'veldisplos_0.txt')
    filesig.append(dir+'siglos/'+simname+'veldisplos_1.txt')
    filekappa.append(dir+'kappalos/'+simname+'kappalos_0.txt')
    filekappa.append(dir+'kappalos/'+simname+'kappalos_1.txt')
    
     #################### files for Walker's mock data ####################
elif gp.investigate == 'walk': # or just want to try some other generic pymc stuff:
    # prior:  analysis is done out to rprior*r_core (half-light radius)
    # if < 0: include all particles, r=max(x**2+y**2)
    rprior = 3.0
    r_DM  = 1000.
    
    def rhodm(r):
        exp = ((1.*gamma_DM-beta_DM)/alpha_DM)
        rho = rho0*(r/r_DM)**(-gamma_DM)*(1.+(r/r_DM)**alpha_DM)**exp
        return rho
    
    fi = Files()
    fi.set_walk()
    dir = fi.dir
    fil = dir+'mem2'
    
    ncomp = gp.pops+1    # do analysis for all of them, plus one for all comp's
    
    pmsplit = 0.9 # minimum probability of membership required for analysis
    # use 0 if grw_* should be called from within gravlite
    fileposcartesian = dir+'simulation/pos.txt'
    filevelcartesian = dir+'simulation/vel_my.txt'
    
elif gp.investigate == 'gaia':
    fi  = Files()
    fi.set_gaia()
    dir = fi.dir
    fil = dir+'dat'
    r_DM = 1000.
    rprior = 3.                 # take data out to rprior*R_s
    ncomp  = 2                  # TODO: search Gaia data for two-components, or find out how to split them

elif gp.investigate == 'triax':
    fi  = Files()
    dir = fi.dir
    fil = dir + fi.set_triax()+'.dat'
    r_DM = 1500.                  # [pc] TODO check
    rprior = 3.                   # [r_DM?]
    
Rerror  = 0.1      # distance error in [Rscale]
vrerror = 0.1      # [km/s] 0.01 # velocity error. only raises sig_los


import os; import os.path
def newdir(bname):
    if not os.path.exists(bname):
        os.makedirs(bname)
    return
## \fn newdir(bname)
# create new directory
# @param bname path


def get_com_file(n):
    return dir+'centeredpos_%i.txt' %(n)
## \fn get_com_file(n)
# get filename of COM file
# @param n population


def get_com_png(n):
    return dir+'centeredpos_%i.png' %(n)
## \fn get_com_png(n)
# get filename of COM png
# @param n population


def get_params_file(n):
    return dir+'scale_%i.txt' %(n)
## \fn get_params_file(n)
# get scale file
# @param n population


def get_ntracer_file(n):
    return dir+'ntracer_%i.txt' %(n)
## \fn get_ntracer_file(n)
# get file with tracer information
# @param n population


def get_pos_sphere_file(n):
    return dir+'sphericalpos_%i.txt' %(n)
## \fn get_pos_sphere_file(n)
# get data file with spherical positions
# @param n population


def get_vel_sphere_file(n):
    return dir+'sphericalvel_%i.txt' %(n)
## \fn get_vel_sphere_file(n)
# get data file with velocities in spherical coordinates
# @param n population


def get_enc_mass_file(n):
    newdir(dir+'enclosedmass/')
    return dir+'enclosedmass/enclosedmass_%i.txt' %(n)
## \fn get_enc_mass_file(n)
# get data files with enclosed mass
# @param n population


def get_dens_file(n):
    # TODO: use gl_class_files.nufiles here!
    newdir(dir+'nu/')
    return dir+'nu/nunotnorm_%i.txt' %(n)
## \fn get_dens_file(n)
# get data file with tracer density falloff
# @param n population


def get_dens_png(n):
    return dir+'nu/nunotnorm_%i.png' %(n)
## \fn get_dens_png(n)
# get png output file
# @param n population


def get_siglos_file(n):
    newdir(dir+'siglos/')
    return dir+'siglos/siglos_%i.txt' %(n)
## \fn get_siglos_file(n)
# get data file with line of sight velocity dispersion
# @param n population


def get_siglos_png(n):
    return dir+'siglos/siglos_%i.png' %(n)
## \fn get_siglos_png(n)
# get png output file for sigma
# @param n population


def get_kurtosis_file(n):
    newdir(dir+'kappalos/')
    return dir+'kappalos/kappalos_%i.txt' %(n)
## \fn get_kurtosis_file(n)
# get datafile for fourth order of the LOS velocity
# @param n population


def get_kurtosis_png(n):
    return dir+'kappalos/kappalos_%i.png' %(n)
## \fn get_kurtosis_png(n)
# get output png filename for kappa
# @param n population
