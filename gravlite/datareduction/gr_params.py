#!/bin/env python3

##
# @file
# global params for data analysis step 1 used to generate input for spherical MCMC

# (c) 2013 ETHZ Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import pdb
import gl_params
gp = gl_params.Params()
from gl_class_files import *

showplots = False
n = 200 # number of iterations in MCMC_*
# nbins = gp.nipol


if gp.investigate == 'hern':
    repr  = 1     # choose simulation representation
    ncomp = gp.pops+1     # number of populations + 1 (for all comps together)
    # numbers of tracers. both 0: take all particles
    
    munit = 1.                                 # [Msun]
    Rcut = 1.e10                               # [Rvir]
    
    Rmin = 0. * gp.ascale; Rmax = 3.*gp.ascale # [pc]
    
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
    r_DM  = 1000.
    
    def rhodm(r):
        exp = ((1.*gamma_DM-beta_DM)/alpha_DM)
        rho = rho0*(r/r_DM)**(-gamma_DM)*(1.+(r/r_DM)**alpha_DM)**exp
        return rho
    
    fi = Files(gp)
    fi.set_walk(gp)
    dir = fi.dir
    fil = dir+'mem2'
    
    ncomp = gp.pops+1    # do analysis for all pops, plus one for all comp's
    
    pmsplit = 0.9 # minimum probability of membership required for analysis
    # use 0 if grw_* should be called from within gravlite
    fileposcartesian = dir+'simulation/pos.txt'
    filevelcartesian = dir+'simulation/vel_my.txt'
    
elif gp.investigate == 'gaia':
    fi  = Files()
    fi.set_gaia()
    dir = fi.dir
    fil = dir + 'dat'
    r_DM = 1000.
    ncomp  = 2 # TODO: search Gaia corresponding two-components

elif gp.investigate == 'triax':
    fi  = Files()
    fi.set_triax(gp)
    dir = fi.dir
    fil = dir + 'dat'
    r_DM = 1500.                  # [pc] TODO check

elif gp.investigate == 'obs':
    fi = Files(gp)
    fi.set_obs(gp)
    dir = fi.dir
    fil = dir+'mem2'
    pmsplit = 0.9
    ncomp = gp.pops+1 # analysis for all components together, and each
                      # population by itself
    
Rerror  = 0.1      # distance error in [Rscale]
vrerror = 0.1      # [km/s] 0.01 # velocity error. only raises sig_los


def get_com_file(n):
    return gp.files.dir+'centeredpos_' + str(n) + '.txt'
## \fn get_com_file(n)
# get filename of COM file
# @param n population


def get_com_png(n):
    return gp.files.dir+'centeredpos_' + str(n) + '.png'
## \fn get_com_png(n)
# get filename of COM png
# @param n population


def get_pos_sphere_file(n):
    return gp.files.dir+'sphericalpos_' + str(n) + '.txt'
## \fn get_pos_sphere_file(n)
# get data file with spherical positions
# @param n population


def get_vel_sphere_file(n):
    return gp.files.dir+'sphericalvel_' + str(n) + '.txt'
## \fn get_vel_sphere_file(n)
# get data file with velocities in spherical coordinates
# @param n population


def get_dens_png(n):
    return gp.files.dir+'nu/nunotnorm_' + str(n) + '.png'
## \fn get_dens_png(n)
# get png output file
# @param n population


def get_siglos_png(n):
    return gp.files.dir+'siglos/siglos_' + str(n) + '.png'
## \fn get_siglos_png(n)
# get png output file for sigma
# @param n population


def get_kurtosis_png(n):
    return gp.files.dir+'kappalos/kappalos_' + str(n) + '.png'
## \fn get_kurtosis_png(n)
# get output png filename for kappa
# @param n population
