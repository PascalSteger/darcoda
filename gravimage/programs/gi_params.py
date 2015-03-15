#!/usr/bin/env ipython3

##
# @file
# all parameters for the gravimage MCMC, gaia investigation

# (c) GPL v3 2015 ETHZ Pascal S.P. Steger, pascal@steger.aero

import numpy as np
import pdb
import gi_base as gb

def check_investigate(inv):
    if inv == 'walk':     return True
    if inv == 'gaia':     return True
    if inv == 'obs':      return True
    if inv == 'hern':     return True
    if inv == 'triax':    return True
    if inv == 'coll':     return True
    raise Exception('wrong investigative case in gi_params')
    return False
## \fn check_investigate(inv)
# check whether there is a valid investigation chosen
# @param inv string

class Params():
    def __init__(self, timestamp = '-1', investigate = '', case = -1):
        if timestamp == '-1':
            self.restart = False
        else:
            self.restart = True

        if investigate != '':
            self.investigate = investigate
        else:
            self.investigate  = 'obs' # determine which data set to work on
                                      # 'hern': simple Hernquist prof. from simwiki
                                      # 'walk': with full obs. cont. data from Walker
                                      # 'gaia': 1pop 6D data, gaia challenge
                                      # 'coll': collisional system
                                      # 'obs': real data from Fornax dwarf galaxy
        check_investigate(self.investigate)
        if case != -1:
            self.case = case
            #os.system('sed -i "s/case = 1/case = '+str(case)+'/"')
        else:
            self.case = 5 # gaia models (1..8) Walker (0..2,4,5; use 1, 2)
                          # triax (1-4:core, 5-8:cusp), obs (1:for,car,scl,sex,dra)

        self.pops = 2 # number of stellar tracer populations # if changed: set getnewdata=True!
        # Set number of tracer stars to look at
        self.ntracer = [1e6, 1e6] # pop0, pop1, pop2, ..., pop_N

        # data options
        self.getnewdata = True # new data computed from observations before burn-in
        if self.restart: self.getnewdata = False
        self.selfconsistentnu = False # tracer star density profile for dSph?
        self.binning = 'consttr' # linspace, logspace, consttr: binning of particles
        self.metalpop = True # split metallicities with a separate MCMC
        self.Rdiff = 'median' # median, min, max for median,
        self.walker3D = False # for walker mock data: use 3D models
        self.hern_dual = 2 # use hernquist model with 1 or 2 particle
                     # types. do not use second type (DM) as population
        self.maxR = 5. # [Xscale], max range in radial bins

        # MultiNest options
        self.getSigdata = False # get previously stored parameters for nu,  after a Sig convergence run
        self.chi2_switch = 80   # turn on sig calculation if chi2 < chi2_switch
        self.chi2_Sig_converged = 1000 # how many times to be below that threshold?
        # Set number of terms for enclosedmass&tracer&anisotropy bins = model parameters:
        self.nipol = 12  # set getnewdata = True to run data readout again if this value changes
        self.nexp  = 3    # more fudge parameters at r<rmin and r>rmax
        self.nepol = self.nipol + 2*self.nexp     # number of parameters for
                                                # direct mapping of nu(r)
        self.nrho = self.nipol + 2*self.nexp + 3 # +3 means 1 more
                                                # parameter for the density at half
                                                # light radius, 1 more parameter for the
                                                # asymptote to 0, 1 more parameter for
                                                # the asymptote to \infty
        self.nfine = 10*self.nipol  # number of entries in integral lookup table
                         # gives no. log spaced points
        self.rinfty = 5. # interpolate from last slope to slope at
                         # 10*max(xipol), where asymptote to \infty
                         # is reached, must be >= 11
        self.nbeta = 4   # number of parameters for beta, in sum of
                         # polynomials
        self.x0turn = -1 # [pc] pinch radius for beta polynomial, set in data readin
        # next: # live points, > ndim, < 2^ndim, about number of
        # ellipsoids in phase space to be found
        self.geom = 'sphere'
        # N_nu = number of populations to work with
        if self.investigate == 'obs':
            N_nu = (self.pops+1)*self.nrho+1 #last +1 for MtoL param
        else:
            N_nu = self.pops*self.nrho
        self.ndim = self.nrho + N_nu + self.pops*self.nbeta
        self.nlive = 10*self.ndim
        self.err = 1e300    # chi^2 for models which are impossible

        # parameter spaces
        # --------------------------------------------------------
        self.rhohalf = -1.    # prior density for rho at half-light radius of tracers
                              # calculated in gi_data, in linear space
        self.log10rhospread = 1.  # with this spread, [dex] in log space
        self.log10nuspread = 0.3  # same for nu
        self.rlimnr = 1       # radius in [Rhalf] below which n(r) is bounded by maxrhoslope/2
        self.rlimnr_nu = 1    # same for nrnu
        self.innerslope = 2.999
        self.maxnuslope = 5   # same for nrnu
        self.nrtol  = 1.5     # scale of change of dn/dr
        self.nrtol_nu = self.maxnuslope*2 # scale of change for nr_nu
        self.nupar_min = np.zeros(self.nrho)  # ranges to be sampled
        self.nupar_max = np.ones(self.nrho)*self.nrtol_nu
        self.beta00prior = False  # beta(r=0)=0
        self.minbetastar_0 = -0.99  # clipping for beta*, default: -0.99
        self.maxbetastar_0 = 1.00  # clipping for beta*, default:  1.00
        self.minbetastar_inf = -0.99
        self.maxbetastar_inf = 1.00
        self.betalogrs = 5.86803451 # for checkbeta, fitted value for Gaia02
        self.MtoLmin = 0.8
        self.MtoLmax = 3.
        self.monotonic = False    # monotonicity prior on n(x) for rho(x)
        self.monotonic_nu = False # monotonicity prior on n(x) for nu(x)

        # integration options
        # --------------------------------------------------
        self.usekappa   = False # switch to turn on (True) or off the
                                # calculation of kappa
        self.usezeta    = False # switch to turn on (True) or off the
                                # calculation of virial parameters zeta_a,b

        # automatic plotting options
        # --------------------------------------------------
        self.last_plot = -1    # timestamp of last automatic plot, set to -1
        self.plot_after = 36000  # [s] to elapse before automatic plotting called again

        # filesystem-related
        # --------------------------------------------------
        self.machine = gb.get_machine()
        import import_path as ip
        ip.set_geometry(self.geom, self.machine) # load spherical or disc version of the code
        import gi_class_files
        self.files = gi_class_files.Files(self, timestamp)
        from gi_data import Datafile
        self.dat = Datafile()

        # debug options
        # --------------------------------------------------
        self.debug = False # enable calling debug routines. Turn off in production runs!
        self.checkbeta = False # check that if right r_s and beta(r_infty) is set,
                              # we get the right profiles back
        self.checksig = False # check sigma_LOS calculation steps in gi_int
        self.stopstep = 5 # step to stop at by default

        # global arrays
        # --------------------------------------------------
        self.xipol = np.array([]) # [pc] hold radius bin centers
        self.xepol = np.array([]) # [pc] extended by 3 fudge bins
        self.xfine = np.array([]) # [pc] radii for lookup tables,
                                  #      gp.nfine long
        # scaling: Xscale in [pc], surfdens_central (=Sig0) in
        # in [Munit/pc^2], and totmass
        # [Munit], and max(v_LOS) in [km/s]
        self.rscale = []
        self.Xscale = []
        self.Sig0pc = []
        self.nu0pc = []
        self.totmass_tracers = []
        self.maxsiglos = []
        self.ana = 1000.
        self.anM = 1.
        if self.investigate == 'hern':
            # each is set for all components and first component by
            # default
            self.rscale.append(1.);           self.rscale.append(1.)
            self.Xscale.append(1.);           self.Xscale.append(1.)
            self.nu0pc.append(1.);            self.nu0pc.append(1.)
            self.Sig0pc.append(1.);           self.Sig0pc.append(1.)
            self.totmass_tracers.append(1.);          self.totmass_tracers.append(1.)
            self.maxsiglos.append(1.);        self.maxsiglos.append(1.)
    ## \fn __init__(self, timestamp, investigate, case)
    # set up all parameters used in the course of the MultiNest run,
    # and the analysis routines in plot_multinest
    # @param timestamp = '', for new output. set to some existing output dir for restart, plotting
    # @param investigate if set, string of mode
    # @param case if set, number of case

    def __repr__(self):
        return "Params: "+self.files.dir
    ## \fn __repr__(self)
    # string representation for ipython
