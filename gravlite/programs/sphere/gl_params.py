#!/usr/bin/env ipython3

##
# @file
# all parameters for the gravlite MCMC, gaia investigation

# (c) 2013 ETHZ Pascal S.P. Steger

import numpy as np
import os, pdb, logging, socket, getpass

def check_investigate(inv):
    if inv == 'walk':     return True
    if inv == 'gaia':     return True
    if inv == 'obs':      return True
    if inv == 'hern':     return True
    if inv == 'triax':    return True
    raise Exception('wrong investigative case in gl_params')
    return False
## \fn check_investigate(inv)
# check whether there is a valid investigation chosen
# @param inv string


class Params():
    def __init__(self, timestamp = ''):

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

        self.investigate  = 'gaia' # determine which data set to work on
                                  # 'hern': check simple Hernquist prof. from simwiki
                                  # 'walk': check with full obs. cont. data from Walker
                                  # 'gaia': 6D data (x,y,z,vx,vy,vz) from gaia
                                  #         challenge, 1 pop only
                                  # 'obs': real data from Fornax dwarf galaxy
        check_investigate(self.investigate)
        self.case = 5 # gaia models (1..8) Walker (0..2,4,5; use 1, 2)
                      # triax (1-4:core, 5-8:cusp)
        self.pops = 1 # number of stellar tracer populations
                      # if changed: set getnewdata=True!
        # Set number of tracer stars to look at take all particles #
        # case 0 want to set ntracer = 3e3 # case 1 ntracer = 1e4 #
        # case 2
        self.ntracer = [1e6, 1e6] # pop0, pop1, pop2, ..., pop_N

        # data options
        # ----------------------------------------------------------------------
        self.getnewdata = True # get new data computed from observations before burn-in
        self.getnewpos  = True # get new positions of particles, important for Hernquist runs
        self.binning = 'consttr' # 'linspace', 'logspace', 'consttr': binning of particles
        self.metalpop   = False # split metallicities with a separate
                                # MCMC
        self.walker3D = False # for walker mock data: use 3D models
        self.hern_dual = 2 # use hernquist model with 1 or 2 particle
                     # types. do not use second type (DM) as
                     # population
        self.maxR = 5.            # [Xscale], max range in radial bins

        # debug options
        # ----------------------------------------------------------------------
        self.debug = False
        self.checksig = False
        self.stopstep = 1

        # MultiNest options
        # ----------------------------------------------------------------------
        self.chi2_Sig_converged = False # set to False to first converge on Sig
        self.chi2_switch = 10.
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
        self.nfine = 2*self.nipol  # number of entries in integral lookup table
                         # gives no. log spaced points
        self.rinfty = 5. # interpolate from last slope to slope at
                          # 10*max(xipol), where asymptote to \infty
                          # is reached, must be >= 11
        self.nbeta = 4    # number of parameters for beta, in sum of
                         # polynomials
        self.x0turn = -1 # [pc] pinch radius for beta polynomial, set in data readin
        # next: # live points, > ndim, < 2^ndim, about number of
        # ellipsoids in phase space to be found
        self.geom = 'sphere'
        if self.investigate == 'obs':
            N_nu = (self.pops+1)*self.nrho+1
        else:
            N_nu = self.pops*self.nrho
        self.ndim = self.nrho + N_nu + self.pops*self.nbeta
        self.nlive = 10*self.ndim
        self.err = 1e300    # chi^2 for models which are impossible

        # spherical case and both cases
        # ----------------------------------------------------------------------
        self.rhohalf = -1.    # prior density for rho at
                                  # half-light radius of tracers
                                  # calculated in gl_data
        self.log10rhospread = 1.       # with this spread, [dex] in log space
        self.log10nuspread = 0.5
        self.rlimnr = 1       # radius in [Rhalf] below which n(r) is bounded by maxrhoslope/2
        self.rlimnr_nu = 1    # same for nrnu
        self.nrtol  = 1./(8./self.nipol) # prior (max +/- range) for dn(r)/dlog(r); 8 is log(3000[pc])
        self.nrtol_nu = 1./(8./self.nipol) # max change in dn(r)/d log(r)
        self.maxrhoslope  = 4    # maximum slope (change if
                                 # monotonicity prior used) of rho
        self.maxnuslope = 5      # same for nrnu
        self.beta00prior = False  # beta(r=0)=0
        self.minbetastar = -0.99  # clipping for beta, default: -0.99
        self.maxbetastar = 0.99    # clipping for beta, default:  1.00
        self.MtoLmin = 0.8
        self.MtoLmax = 3.


        # disc case
        # ----------------------------------------------------------------------
        # norm1 = 17.**2 # offset of sig[0]/nu[0], from int starting
        # at zmin instead of 0 norm2 = 10.**2 # and for the second
        # component, if there is one
        self.quadratic = False    # linear or quad interpol.
        self.monotonic = False    # mono-prior on rho(x)
        self.monotonic_nu = False # mono-prior on nu(x)
        self.adddarkdisc = False  # for disc mock case: add a dark disc?
        self.baryonmodel = 'sim' # read in surface density from
                                 # corresponding surfden file
                                 # 'silvia', 'sim', 'simple'

        # integration options
        # ----------------------------------------------------------------------
        self.even = 'avg' # for simps integration (everywhere): 'avg',
                          # 'first', 'last'
        self.usekappa   = False # switch to turn on (True) or off the
                                # calculation of kappa
        self.usezeta    = False # switch to turn on (True) or off the
                                # calculation of virial parameters zeta_a,b

        # filesystem-related
        # ----------------------------------------------------------------------
        import import_path as ip
        ip.set_geometry(self.geom, self.machine) # load spherical or
                                                 # disc version
                                                 # of the code
        import gl_class_files
        self.files = gl_class_files.Files(self, timestamp)
        from gl_data import Datafile
        self.dat = Datafile()


        # global arrays
        # ----------------------------------------------------------------------
        self.xipol = np.array([]) # [pc] hold radius bin centers
        self.xepol = np.array([]) # [pc] extended by 3 fudge bins
        self.xfine = np.array([]) # [pc] radii for lookup tables,
                                  #      gp.nfine long
        # scaling: Xscale in [pc], surfdens_central (=Sig0) in
        # in [Munit/pc^2], and totmass
        # [Munit], and max(v_LOS) in [km/s]
        self.rscale=[];        self.nu0pc=[]
        self.Xscale=[];        self.Sig0pc=[]
        self.totmass_tracers=[];       self.maxsiglos=[]
        # for investigations without data:
        if self.investigate != 'walk' and\
           self.investigate != 'triax' and\
           self.investigate != 'gaia' and\
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
