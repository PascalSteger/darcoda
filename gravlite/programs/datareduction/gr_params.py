#!/bin/env ipython3

##
# @file
# global params for data analysis step 1 used to generate input for spherical MCMC

# (c) 2014 ETHZ Pascal Steger, psteger@phys.ethz.ch

import pdb

import gl_helper as gh
class Params():
    def __init__(self, gp):
        from gl_class_files import Files
        # show plots during execution of data readout?
        # set automatically if gr_MCMCbin.py is called on the command line
        self.showplots = False

        self.n = 3 # number of iterations in gr_MCMCbin
        self.Rerr  = 0. # 0.01      # distance error in [Xscale]
        self.vrerr = 2.0 # [km/s] 0.01 # velocity error. only raises sig_los

        if gp.investigate == 'hern':
            self.repr  = 1     # choose simulation representation
            self.Rcut = 1.e10  # [Rvir]
            self.Rmin = 0. # [Rscale]i
            self.Rmax = 3. # [Rscale]

            self.simname = gp.files.get_sim_name(gp) # dir+'simulation/'+prename+'unit_hern_%i' %(repr)
            if gp.pops == 1:
                self.simpos = gp.files.dir+'simulation/'+self.simname+'pos.txt'
                self.simvel = gp.files.dir+'simulation/'+self.simname+'vel.txt'
            elif gp.pops == 2:
                self.simpos = gp.files.dir+'simulation/'+self.simname+'stars_pos.txt'
                self.simvel = gp.files.dir+'simulation/'+self.simname+'stars_vel.txt'
            else:
                gh.LOG(0, 'get data for more than 2 pops in Hernquist profile')
                pdb.set_trace()

        elif gp.investigate == 'walk': # or just want to try some other generic pymc stuff:
            self.r_DM  = 1000.

            def rhodm(r):
                # TODO beta_DM: define up here somewhere
                exp = ((1.*gamma_DM-beta_DM)/alpha_DM)
                rho = rho0*(r/r_DM)**(-gamma_DM)*(1.+(r/r_DM)**alpha_DM)**exp
                return rho
            ## \fn rhodm(r)
            # in walker case, calculate rho_DM profile
            # @param r radii [pc]

            self.fi = Files(gp)
            self.fi.set_walk(gp)
            self.dir = self.fi.dir
            self.fil = self.dir+'mem2'

            self.pmsplit = 0.9 # minimum probability of membership required for analysis
            # use 0 if grw_* should be called from within gravlite
            self.fileposcartesian = self.dir+'simulation/pos.txt'
            self.filevelcartesian = self.dir+'simulation/vel_my.txt'

        elif gp.investigate == 'gaia':
            self.fi  = Files(gp)
            self.fi.set_gaia(gp)
            self.dir = self.fi.dir
            self.fil = self.dir + 'dat'
            self.r_DM = 1000.

        elif gp.investigate == 'triax':
            self.fi  = Files(gp)
            self.fi.set_triax(gp)
            self.dir = self.fi.dir
            self.fil = self.dir + 'dat'
            self.r_DM = 1500.                  # [pc]

        elif gp.investigate == 'obs':
            self.fi = Files(gp)
            self.fi.set_obs(gp)
            self.dir = self.fi.dir
            self.fil = self.dir+'mem2'
            self.pmsplit = 0.9
    ## \fn __init__(gp)
    # set up common parameters for data pre-processing
    # @param gp global parameters
