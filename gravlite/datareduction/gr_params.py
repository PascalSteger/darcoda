#!/bin/env python3

##
# @file
# global params for data analysis step 1 used to generate input for spherical MCMC
# (c) 2013 ETHZ Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import pdb
import gl_params as gp
from gl_class_files import *

## run all commands
def run():
    showplots = False
    n = 100
    nbins = gp.nipol
    
    
    if gp.investigate == 'hernquist':
        sim = 1           # choose simulation
        ncomp = 2         # number of populations + 1 (for all comps together)
        # number of populations/tracers. both set to 0 means: take all particles
        pops = 1
        ntracers1 = 10000
        ntracers2 = 0
        
        munit = 1.                                 # [Msun]
        Rcut = 1.e10                               # [Rvir]
        
        bins = gp.nipol                            # set binning
        Rmin = 0. * gp.ascale; Rmax = 3.*gp.ascale # [pc]
        rprior = 3.
        
        nit = 30             # set number of iterations
        procs = 1            # number of processes for parallel execution
        
        dir = gp.files.dir
        # dir = '/home/ast/read/dark/dwarf_data/data_hernquist/'
        # dir = '/home/psteger/sci/dwarf_data/data_hernquist/'
        
        if pops == 2:
            prename = 'dual_'
            
            if ntracers1 >0:
                comp = '_dm'
            else:
                comp = '_stars'
        else:
            prename = ''
            comp = ''
    
        simname = dir+'simulation/'+prename+'unit_hern_%i' %(sim)
        simpos = simname+comp+'_pos.txt'
        simvel = simname+comp+'_vel.txt'
        fileposcartesian = []; fileposcenter = []; filevelcartesian = []; 
        fileposcartesian.append(dir+'simulation/'+prename+'unit_hern_%i_pos_0.txt'%(sim))
        fileposcartesian.append(dir+'simulation/'+prename+'unit_hern_%i_pos_1.txt'%(sim))
        fileposcenter.append(dir+'simulation/'+prename+'unit_hern_%i_centeredpos_0.txt'%(sim))
        fileposcenter.append(dir+'simulation/'+prename+'unit_hern_%i_centeredpos_1.txt'%(sim))
        filevelcartesian.append(dir+'simulation/'+prename+'unit_hern_%i_vel_0.txt'%(sim))
        filevelcartesian.append(dir+'simulation/'+prename+'unit_hern_%i_vel_1.txt'%(sim))
        
        fileposspherical =[]; filevelspherical = []
        fileposspherical.append(dir+'simulation/'+prename+'unit_hern_%i_sphericalpos_0.txt'%(sim))
        fileposspherical.append(dir+'simulation/'+prename+'unit_hern_%i_sphericalpos_1.txt'%(sim))
        filevelspherical.append(dir+'simulation/'+prename+'unit_hern_%i_sphericalvel_0.txt'%(sim))
        filevelspherical.append(dir+'simulation/'+prename+'unit_hern_%i_sphericalvel_1.txt'%(sim))
        
        filemass = []; filedenfalloff = []; filesig = []; filekappa = []
        filemass.append(dir+'enclosedmass/'+prename+'unit_hern_%i_enclosedmass_0.txt'%(sim))
        filemass.append(dir+'enclosedmass/'+prename+'unit_hern_%i_enclosedmass_1.txt'%(sim))
        filedenfalloff.append(dir+'densityfalloff/'+prename+'unit_hern_%i_falloffnotnorm_0.txt' %(sim))
        filedenfalloff.append(dir+'densityfalloff/'+prename+'unit_hern_%i_falloffnotnorm_1.txt' %(sim))
        filesig.append(dir+'siglos/'+prename+'unit_hern_%i_veldisplos_0.txt'%(sim))
        filesig.append(dir+'siglos/'+prename+'unit_hern_%i_veldisplos_1.txt'%(sim))
        filekappa.append(dir+'kappalos/'+prename+'unit_hern_%i_kappalos_0.txt'%(sim))
        filekappa.append(dir+'kappalos/'+prename+'unit_hern_%i_kappalos_1.txt'%(sim))

    #################### files for Walker's mock data ####################
    elif gp.investigate == 'walker': # or just want to try some other generic pymc stuff:
        # prior:  analysis is done out to rprior*r_core (half-light radius)
        # if < 0: include all particles, r=max(x**2+y**2)
        rprior = 3.0
        r_DM  = 1000.
        
        def rhodm(r):
            exp = ((1.*gamma_DM-beta_DM)/alpha_DM)
            rho = rho0*(r/r_DM)**(-gamma_DM)*(1.+(r/r_DM)**alpha_DM)**exp
            return rho
  
        fi = Files()
        fi.set_walker()
        dir = fi.dir
        fil = dir+'mem2'
        
        ncomp = 4  # 3 possibilities: 0 (both), tracer pop 1, tracer pop 2, 
        # 3: foreground contamination
        # a value of 4 means: do analysis for all of them
        
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
        ncomp  = 2
    
    elif gp.investigate == 'triaxial':
        fi  = Files()
        dir = fi.dir
        fil = dir + fi.set_triaxial()+'.dat'
        r_DM = 1500.                  # [pc] TODO check
        rprior = 3.                   # [r_DM?]

    Rerror  = 0.1      # distance error in [Rcore]
    vrerror = 0.1      # [km/s] 0.01 # velocity error. only raises sig_los

    import os; import os.path

## create new directory
# @param bname path
def newdir(bname):
    if not os.path.exists(bname):
        os.makedirs(bname)
    return

## get filename of COM file
# @param n population
def get_com_file(n):
    return dir+'centeredpos_%i.txt' %(n)

## get filename of COM png
# @param n population
def get_com_png(n):
    return dir+'centeredpos_%i.png' %(n)

## get scale file
# @param n population
def get_params_file(n):
    return dir+'scale_%i.txt' %(n)

## get file with tracer information
# @param n population
def get_ntracer_file(n):
    return dir+'ntracer_%i.txt' %(n)

## get data file with spherical positions
# @param n population
def get_pos_sphere_file(n):
    return dir+'sphericalpos_%i.txt' %(n)

## get data file with velocities in spherical coordinates
# @param n population
def get_vel_sphere_file(n):
    return dir+'sphericalvel_%i.txt' %(n)

## get data files with enclosed mass
# @param n population
def get_enc_mass_file(n):
    newdir(dir+'enclosedmass/')
    return dir+'enclosedmass/enclosedmass_%i.txt' %(n)

## get data file with tracer density falloff
# @param n population
def get_dens_file(n):
    newdir(dir+'nu/')
    return dir+'nu/nunotnorm_%i.txt' %(n)

## get png output file
# @param n population
def get_dens_png(n):
    return dir+'nu/nunotnorm_%i.png' %(n)

## get data file with line of sight velocity dispersion
# @param n population
def get_siglos_file(n):
    newdir(dir+'siglos/')
    return dir+'siglos/siglos_%i.txt' %(n)

## get png output file for sigma
# @param n population
def get_siglos_png(n):
    return dir+'siglos/siglos_%i.png' %(n)

## get datafile for fourth order of the LOS velocity
# @param n population
def get_kurtosis_file(n):
    newdir(dir+'kappalos/')
    return dir+'kappalos/kappalos_%i.txt' %(n)

## get output png filename for kappa
# @param n population
def get_kurtosis_png(n):
    return dir+'kappalos/kappalos_%i.png' %(n)
