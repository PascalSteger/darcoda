#!/usr/bin/env ipython3

##
# @file
# all parameters for the gravlite MCMC, gaia investigation

# (c) 2013 ETHZ Pascal S.P. Steger

import numpy as np
import os, ipdb, logging, socket, getpass

def sanitize_investigate(inv):
    if inv == 'walk':     return True
    if inv == 'gaia':     return True
    if inv == 'obs':      return True
    if inv == 'hern':     return True
    if inv == 'triax':    return True
    raise Exception('wrong investigative case in gl_params')
    return False
## \fn sanitize_investigate(inv)
# check whether there is a valid investigation chosen
# @param inv string


class Params():
    def __init__(self, timestamp = ''):
        # basic setup
        # -----------------------------------------------------------
        self.investigate  = 'gaia' # determine which data set to work on
                                  # 'hern': check simple Hernquist prof. from simwiki
                                  # 'walk': check with full obs. cont. data from Walker
                                  # 'gaia': 6D data (x,y,z,vx,vy,vz) from gaia
                                  #         challenge, 1 pop only
                                  # 'obs': real data from Fornax dwarf galaxy
        sanitize_investigate(self.investigate)
        self.case = 3 # gaia models (1..8) Walker (0..2,4,5; use 1, 2)
                      # triax (1-4:core, 5-8:cusp)
        self.pops = 1 # number of stellar tracer populations
        if self.investigate == 'hern':
            self.pops = 1
        self.ntracer = [1e6] # number of tracer stars pop1 (Hernquist case), pop2, ...



        # data options
        # ------------------------------------------------------------
        self.getnewdata = True     # get new data computed from
                                    # observations before burn-in
        self.getnewpos  = True      # read in the positions and v_LOS again
        if self.getnewdata == False: self.getnewpos = False
        self.binning    = 'consttr' # linspace, logspace, consttr
        self.metalpop   = False     # split metallicities with a separate MCMC
        self.usekappa   = False # switch to turn on (True) or off the
                                # calculation of kappa
        self.usezeta    = False # switch to turn on (True) or off the
                                # calculation of virial parameters zeta_a,b
        self.walker3D   = False     # for walker mock data: use 3D models
        self.hern_sim_pops = 1 # use hernquist model with 1 or 2 particle
                               # types. do not use second type (DM) as population
        if self.pops == 1: self.hern_sim_pops = 1
        self.maxR       = 5.        # [Xscale], max range in radial bins



        ########## MultiNest options
        # ----------------------------------------------------------------------
        self.chi2_Sig_converged = False # set to False to first converge on Sig
        self.chi2_switch = 10.          # if 10chi^2>chi2_switch, switch sig calc on
        # Set number of terms for enclosedmass+tracer+anisotropy bins
        # = model parameters:
        self.nipol = 12   # IF CHANGED => set getnewdata = True to run
                         # data readout again
        self.nexp  = 3    # more fudge parameters at r<rmin and r>rmax
        self.nepol = self.nipol + 2*self.nexp     # number of parameters for
                                                # direct mapping of nu(r)
        self.nrho = self.nipol + 2*self.nexp + 3 # +3 means 1 more
                                                # parameter for the
                                                # density at half
                                                # light radius, 1 more
                                                # parameter for the
                                                # asymptote to 0, 1
                                                # more parameter for
                                                # the asymptote to
                                                # \infty
        self.nfine = 100 #2*self.nipol  # number of entries in integral lookup table
                         # gives no. log spaced points
        self.rinfty = 10 # interpolate from last slope to slope at
                          # 10*max(xipol), where asymptote to \infty
                          # is reached, must be >= 11
        self.nbeta = 4   # number of parameters for beta
        # next: # live points, > ndim, < 2^ndim, about number of
        # ellipsoids in phase space to be found
        self.geom = 'sphere'
        if self.investigate == 'obs':
            N_nu = (self.pops+1)*self.nrho+1
        else:
            N_nu = self.pops*self.nrho
        self.ndim = self.nrho + N_nu + self.pops*self.nbeta
        self.nlive = 2*self.ndim
        self.err = 1e300    # chi^2 for models which are impossible


        ########## spherical case and both cases
        # ----------------------------------------------------------------------
        self.rhohalf = -1.        # prior density for rho at half-light radius of tracers
                                  # calculated in gl_data
        self.log10rhospread = 6.       # with this spread, [dex] in log space
        self.nuspread = 1.0       # analog for nu profile
        self.rlimnr = 1 # scale below which range of
                         # n(r<rlimnr*r_half)<maxrhoslope/2, in multiples of r_half.
                         # calculated to values in [pc] in gl_data.read_nu;
                         # if set to -1 here, use maxrhoslope everywhere
        self.rlimnr_nu = 1 # same for nu, using same rhalf
        self.maxrhoslope  = 5    # maximum slope (change if
                                 # monotonicity prior used) of rho
        self.maxrhoslope_nu = 5
        # prior (max +/- range) for dn(r)/dlog(r)
        #   determine how far nr can wander with the max allowed nr slope
        #    on from min(gp.xipol) to max(gp.xipol)
        self.nrtol  = 2*self.maxrhoslope
        self.nrtol_nu = 2*self.maxrhoslope_nu # same for nu profile
        self.maxbetaslope = 1.5   # linear (and 2nd..order) max slope
                                  # of beta* in polynomial representation
        self.minbetastar = -0.99  # clipping for beta, default: -0.99
        self.maxbetastar = 0.99    # clipping for beta, default:  1.00
        self.beta00prior = False  # prior beta(r=0) = 0
        self.MtoLmin = 0.8        # boundaries for flat prior on constant mass-to-light ratio
        self.MtoLmax = 3.
        self.monotonic = False    # monotonicity-prior on nr_rho(x)
        self.monotonic_nu = False # monotonicity-prior on nr_nu(x)


        ########## debug options
        # ----------------------------------------------------------------------
        self.checksig   = False  # check sigma calculation routine with 'walk'
        self.stopstep   = 1     # stop after step number ..., enter debugger


        # unitsXS
        # ----------------------------------------------------------------------
        self.G1  = 6.67398e-11                # [m^3 kg^-1 s^-2]
        self.pc_in_m  = 3.08567758e16              # [m]
        self.msun= 1.981e30                   # [kg]
        self.km_in_m  = 1000.                      # [m]
        self.kpc = 1000.                      # [pc]
        self.G1  = self.G1*self.msun/self.km_in_m**2/self.pc_in_m
        # [pc msun^-1 (km/s)^2]
        if self.investigate == 'hern':
            self.G1 = 1.            # as defined by Justin, so we can rescale model
            self.ana        = 1.    # scale radius of Hernquist profile in [pc]
            self.anM        = 1.    # total mass of Hernquist profile in [Msun]


        ########## filesystem-related
        # ----------------------------------------------------------------------
        host_name = socket.gethostname()
        user_name = getpass.getuser()
        if 'darkside' in host_name:
            self.machine = 'darkside'
        elif 'pstgnt332' in host_name:
            self.machine = 'pstgnt332'
        elif ('lisa' in host_name) and ('hsilverw' in user_name):
            self.machine = 'lisa_HS'
        elif ('lisa' in host_name) and ('sofia' in user_name):
            self.machine = 'lisa_SS'

        import import_path as ip
        ip.set_geometry(self.geom, self.machine) # load spherical or
                                                 # disc version
                                                 # of the code
        import gl_class_files
        self.files = gl_class_files.Files(self, timestamp)
        from gl_data import Datafile
        self.dat = Datafile()


        ########## global arrays
        # ----------------------------------------------------------------------
        self.xipol = np.array([]) # [pc] hold radius bin centers
        self.xepol = np.array([]) # [pc] extended by 3 fudge bins
        self.xfine = np.array([]) # [pc] radii for lookup tables,
                                  #      gp.nfine long
        # scaling: Xscale in [pc], surfdens_central (=Sig0) in
        # in [Munit/pc^2], and totmass_tracers
        # [Munit], and max(sigma_LOS) in [km/s]
        self.rscale=[];        self.nu0pc=[]
        self.Xscale=[];        self.Sig0pc=[]
        self.totmass_tracers=[];       self.maxsiglos=[]
        # for investigations without data:
        if self.investigate != 'walk' and\
           self.investigate != 'triax' and\
           self.investigate != 'gaia' and\
           self.investigate != 'hern' and\
           self.investigate != 'discmock':
            # each is set for all components and first component by
            # default
            self.rscale.append(1.);           self.rscale.append(1.)
            self.Xscale.append(1.);           self.Xscale.append(1.)
            self.nu0pc.append(1.);            self.nu0pc.append(1.)
            self.Sig0pc.append(1.);           self.Sig0pc.append(1.)
            self.totmass_tracers.append(1.);          self.totmass_tracers.append(1.)
            self.maxsiglos.append(1.);        self.maxsiglos.append(1.)

    ## \fn __init__(self, timestamp = '')
    # set up all parameters used in the course of the MultiNest run,
    # and the analysis routines in plot_multinest @param timestamp =
    # '', for output


    def __repr__(self):
        return "Params: "+self.files.dir
    ## \fn __repr__(self)
    # string representation for ipython
