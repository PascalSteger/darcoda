#!/bin/env ipython3

##
# @file
# global params for data analysis step 1 used to generate input for spherical MCMC

# (c) GPL v3 2015 ETHZ Pascal Steger, pascal@steger.aero

import pdb

import gi_helper as gh
class grParams():
    def __init__(self, gp):
        from gi_class_files import Files
        # show plots during execution of data readout?
        # set automatically if gr_MCMCbin.py is called on the command line
        self.showplots = False
        self.n = 300 # number of iterations in gr_MCMCbin
        self.Rerr  = 0. # 0.01      # distance error in [Xscale]
        self.vrerr = 2.0 # [km/s] 0.01 # velocity error. raises sig_los, errors in it
        if gp.investigate == 'hern':
            self.repr  = 1     # choose simulation representation
            self.Rcut = 1.e10  # [Rvir]
            self.Rmin = 0. # [Rscale]i
            self.Rmax = 10. # [Rscale]
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
            self.fi = Files(gp)
            self.fi.set_walk(gp)
            self.dir = self.fi.dir
            self.fil = self.dir+'mem2'

            self.pmsplit = 0.9 # minimum probability of membership required for analysis
            # use 0 if grw_* should be called from within gravimage
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

        elif gp.investigate == 'coll':
            self.fi = Files(gp)
            self.fi.set_coll(gp)
            self.dir = self.fi.dir
            self.fil = self.dir+'dat'

    ## \fn __init__(gp)
    # set up common parameters for data pre-processing
    # @param gp global parameters
