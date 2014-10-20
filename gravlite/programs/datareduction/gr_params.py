#!/bin/env ipython3

##
# @file
# global params for data analysis step 1 used to generate input for spherical MCMC

# (c) 2013 ETHZ Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import pdb
from pylab import *
ion()

import gl_params
gp = gl_params.Params()
from gl_class_files import Files
import gl_helper as gh


class Params():
    def __init__(gp):
        # show plots during execution of data readout?
        # set automatically if gr_MCMCbin.py is called on the command line
        showplots = False

        n = 3 # number of iterations in gr_MCMCbin
        Rerr  = 0. #0.01      # distance error in [Xscale]
        vrerr = 0. #2.0      # [km/s] 0.01 # velocity error. only raises sig_los

        if gp.investigate == 'hern':
            repr  = 1     # choose simulation representation
            Rcut = 1.e10  # [Rvir]
            Rmin = 0. # [Rscale]i
            Rmax = 3. # [Rscale]

            simname = gp.files.get_sim_name(gp) # dir+'simulation/'+prename+'unit_hern_%i' %(repr)
            if gp.pops == 1:
                simpos = gp.files.dir+'simulation/'+simname+'pos.txt'
                simvel = gp.files.dir+'simulation/'+simname+'vel.txt'
            elif gp.pops == 2:
                simpos = gp.files.dir+'simulation/'+simname+'stars_pos.txt'
                simvel = gp.files.dir+'simulation/'+simname+'stars_vel.txt'
            else:
                gh.LOG(0, 'get data for more than 2 pops in Hernquist profile')
                pdb.set_trace()

        elif gp.investigate == 'walk': # or just want to try some other generic pymc stuff:
            r_DM  = 1000.

            def rhodm(r):
                exp = ((1.*gamma_DM-beta_DM)/alpha_DM)
                rho = rho0*(r/r_DM)**(-gamma_DM)*(1.+(r/r_DM)**alpha_DM)**exp
                return rho
            ## \fn rhodm(r)
            # in walker case, calculate rho_DM profile
            # @param r radii [pc]

        fi = Files(gp)
        fi.set_walk(gp)
        dir = fi.dir
        fil = dir+'mem2'

        pmsplit = 0.9 # minimum probability of membership required for analysis
        # use 0 if grw_* should be called from within gravlite
        fileposcartesian = dir+'simulation/pos.txt'
        filevelcartesian = dir+'simulation/vel_my.txt'

    elif gp.investigate == 'gaia':
        fi  = Files(gp)
        fi.set_gaia(gp)
        dir = fi.dir
        fil = dir + 'dat'
        r_DM = 1000.

    elif gp.investigate == 'triax':
        fi  = Files(gp)
        fi.set_triax(gp)
        dir = fi.dir
        fil = dir + 'dat'
        r_DM = 1500.                  # [pc]

    elif gp.investigate == 'obs':
        fi = Files(gp)
        fi.set_obs(gp)
        dir = fi.dir
        fil = dir+'mem2'
        pmsplit = 0.9
    ## \fn __init__(gp)
    # set up common parameters for data pre-processing
    # @param gp global parameters


def determine_radius(R, Rmin, Rmax, gp):
    if gp.binning == 'linspace':
        gh.LOG(2, ' bins in linear spacings: ', gp.nipol)
        return gh.bin_r_linear(Rmin, Rmax, gp.nipol)
    elif gp.binning == 'logspace':
        gh.LOG(2, ' bins in log spacings: ', gp.nipol)
        return gh.bin_r_log(Rmax/gp.nipol, Rmax, gp.nipol)
    elif gp.binning == 'consttr':
        gh.LOG(2, ' particles per bin: ', len(R)/gp.nipol)
        return gh.bin_r_const_tracers(R, gp.nipol)
    else:
        gh.LOG(1, 'unknown gp.binning in gpr.determine_radius')
        pdb.set_trace()
## \fn determine_radius(R, Rmin, Rmax, gp)
# determine bin radii once and for all. this must not be changed between
# readout and gravlite run. if you wish to change: set gp.getnewdata =
# True in gl_params.py
# @param R radii of all tracer particles
# @param Rmin float
# @param Rmax float
# @param gp global parameters
