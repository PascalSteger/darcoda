#!/usr/bin/env ipython3

##
# @file
# all parameters for the gravlite MCMC, gaia investigation

# (c) 2013 ETHZ Pascal S.P. Steger

import numpy as np
import os, pdb, logging, socket, getpass

def check_investigate(inv):
    if inv == 'discmock': return True
    if inv == 'discsim': return True
    raise Exception('wrong investigative case in gl_params')
    return False

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

        self.geom = 'disc'
        self.investigate  = 'discmock' # determine which data set to work on
                                  # 'discmock': set up simple model for disc
                                  # 'discsim': read in disc simulation
        check_investigate(self.investigate)

        self.case = 0 # used in spherical case
        self.pops = 2 # number of stellar tracer populations
                      # if changed: set getnewdata=True!
        self.checksig = False # debug sig calculation?
        
        # Set number of tracer stars to look at take all particles #
        # case 0 want to set ntracer = 3e3 # case 1 ntracer = 1e4 #
        # case 2
        self.ntracer = [10000, 10000] # pop1, pop2, ..., pop_N (and
                                      # take sum for all tracers)

        # unitsXS
        self.G1  = 6.67398e-11                # [m^3 kg^-1 s^-2]
        self.pc  = 3.08567758e16              # [m]
        self.msun= 1.981e30                   # [kg]
        self.km  = 1000.                      # [m]
        self.G1  = self.G1*self.msun/self.km**2/self.pc # [pc msun^-1
                                                        # (km/s)^2]

        ########## data options
        self.getnewdata = True # get new data computed from
                                # observations before burn-in
        self.consttr    = True  # set radial bin by constant number of
                                # tracer particles
        self.metalpop   = False # split metallicities with a separate
                                # MCMC
        self.walker3D = False # for walker mock data: use 3D models
        self.hern_dual = 2 # use hernquist model with 1 or 2 particle
                     # types. do not use second type (DM) as
                     # population
        self.maxR = 5.            # [Xscale], max range in radial bins



        ########## MultiNest options
        # Set number of terms for enclosedmass+tracer+anisotropy bins
        # = model parameters:
        self.chi2_Sig_converged = False # first converge on Sig if set to False
        self.chi2_switch = 100. # if chi2<chi2_switch, add chi2_sig 

        self.nipol = 15   # IF CHANGED => set getnewdata = True to run
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
        self.nfine = 30  # number of entries in integral lookup table
                         # gives no. log spaced points
        self.rinfty = 5. # interpolate from last slope to slope at
                          # 10*max(xipol), where asymptote to \infty
                          # is reached, must be >= 11
        self.nbeta = 1   # number of parameters for beta, in sum of
                         # polynomials


        # TODO: if not using rhostar, subtract gp.nrho:
        self.ndim = 1 + 2*self.nrho + 1 + self.pops*(self.nrho + self.nbeta)
        # +1 for norm, rho, rhostar, MtoL, [nu_i, tilt_i]

        # live points, > ndim, < 2^ndim, about number of
        # ellipsoids in phase space to be found
        self.nlive = 2*self.ndim

        self.err = 1e300    # chi^2 for models which are impossible


        ########## spherical case
        self.rhohalf = 10**-1.    # prior density for rho at
                                  # half-light radius of tracers
                                  # calculated in gl_data
        self.rhospread = 1.       # with this spread, [dex] in log space
        self.nuspread = 0.5
        self.iscale = -1 # scale below which range of
                                         # n(r)<2. instead of
                                         # maxrhoslope; is adapted in
                                         # gl_data.read_nu; if set to
                                         # -1 here, use maxrhoslope
                                         # everywhere
        self.iscale_nu = -1

        self.nrtol  = 1.5/(8./self.nipol) # prior (max +/- range) for dn(r)/dlog(r); 8 is log(3000[pc])
        self.nrtol_nu = 1.5/(8./self.nipol)

        self.maxrhoslope  = 0.5    # maximum slope (change if
                                  # monotonicity prior used) of rho
        self.maxrhoslope_nu = 0.5
        self.maxlog10nu = 4.     # direct sampling of nu: min value
        self.minlog10nu = 0.     # direct sampling of nu: max value
        self.maxbetaslope = 0.2   # linear (and 2nd..order) max slope
                                  # of beta*
        self.MtoLmin = 0.8
        self.MtoLmax = 3.


        ########## disc case
        # norm1 = 17.**2 # offset of sig[0]/nu[0], from int starting
        # at zmin instead of 0 norm2 = 10.**2 # and for the second
        # component, if there is one
        self.quadratic = False    # linear or quad interpol.
        self.monotonic = False    # mono-prior on rho(z)
        self.monotonic_nu = False # mono-prior on nu(z)
        self.adddarkdisc = False  # for disc mock case: add a dark disc?

        self.baryonmodel = 'sim' # read in surface density from
                                 # corresponding surfden file
                                 # 'silvia', 'sim', 'simple'

        ########## integration options
        self.even = 'avg' # for simps integration (everywhere): 'avg',
                          # 'first', 'last'
        self.usekappa   = False # switch to turn on (True) or off the
                                # calculation of kappa
        self.usezeta    = True  # switch to turn on (True) or off the
                                # calculation of virial parameters zeta_a,b


        ########## filesystem-related
        import import_path as ip
        ip.set_geometry(self.geom, self.machine) # load spherical or
                                                 # disc version
                                                 # of the code
        import gl_class_files
        self.files = gl_class_files.Files(self, timestamp)
        from gl_data import Datafile
        self.dat = Datafile()


        ########## global arrays
        self.xipol = np.array([]) # [pc] hold radius bin centers
        self.xepol = np.array([]) # [pc] extended by 3 fudge bins
        self.xfine = np.array([]) # [pc] radii for lookup tables, 
                                  #      gp.nfine long

        # scaling: Xscale in [pc], surfdens_central (=Sig0) in
        # in [Munit/pc^2], and totmass
        # [Munit], and max(v_LOS) in [km/s]
        self.rscale=[];        self.nu0pc=[]
        self.Xscale=[];        self.Sig0pc=[]
        self.totmass=[];       self.maxsiglos=[]

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
            self.totmass.append(1.);          self.totmass.append(1.)
            self.maxsiglos.append(1.);        self.maxsiglos.append(1.)

    ## \fn __init__(self, timestamp = '')
    # set up all parameters used in the course of the MultiNest run,
    # and the analysis routines in plot_multinest @param timestamp =
    # '', for output
    

    def __repr__(self):
        return "Params: "+self.files.dir
    ## \fn __repr__(self)
    # string representation for ipython

