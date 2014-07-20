#!/usr/bin/env ipython3

##
# @file
# all parameters for the gravlite MCMC, gaia investigation

# (c) 2013 ETHZ Pascal S.P. Steger

import numpy as np
import os, pdb, logging

class Params():
    def __init__(self, timestamp = ''):
        self.machine = 'darkside'
        if not os.path.exists("/home/ast/read"):
            self.machine = 'local'    # machine: 'local' or 'darkside'

        self.investigate  = 'discmock' # determine which data set to work on
                                  # 'discmock': set up simple model for disc
                                  # 'discsim': read in disc simulation
                                  # 'hern': check simple Hernquist prof. from simwiki
                                  # 'walk': check with full obs. cont. data from Walker

                        # 'gaia': 6D data (x,y,z,vx,vy,vz) from gaia
                        #         challenge, 1 pop only 'obs': real
                        #         data from Fornax dwarf galaxy

        self.case = 0 # gaia models (1..8) Walker (0..2,4,5; use 1 or 2)
                      # triax (0:core, 1:cusp)

        import import_path as ip
        ip.set_geometry(self.investigate, self.machine) # load
                                                        # spherical or
                                                        # disc version
                                                        # of the code
        
        # Set number of tracer stars to look at take all particles #
        # case 0 want to set ntracer = 3e3 # case 1 ntracer = 1e4 #
        # case 2
        self.ntracer = [10000, 10000] # pop1, pop2, ..., pop_N (and
                                      # take sum for all tracers)

        self.getnewdata = True # get new data computed from
                                # observations before burn-in
        self.metalpop   = False # split metallicities with a separate
                                # MCMC
        self.consttr = True     # set radial bin by constant number of
                                # tracer particles

        # unitsXS
        self.G1  = 6.67398e-11                # [m^3 kg^-1 s^-2]
        self.pc  = 3.08567758e16              # [m]
        self.msun= 1.981e30                   # [kg]
        self.km  = 1000.                      # [m]
        self.G1  = self.G1*self.msun/self.km**2/self.pc # [pc msun^-1
                                                        # (km/s)^2]
        self.geom = 'sphere'
        if self.investigate == 'discmock' or self.investigate == 'discsim':
            self.geom = 'disc'

        self.err = 1e30    # chi^2 for models which are impossible
        self.model = False # for walker mock data: plot model set beta
                           # to Osipkov-Merrit profile
        if self.investigate != 'walk':
            self.model = False

        ########## integration options
        self.densint = False # use integral of dens to find mass, not
                             # binned sum
        self.even = 'avg' # for simps integration (everywhere): 'avg',
                          # 'first', 'last'
        self.usekappa   = False # switch to turn on (True) or off the
                                # calculation of kappa
        self.usezeta    = True  # switch to turn on (True) or off the
                                # calculation of virial parameters zeta_a,b
        self.sim = 2 # use hernquist model with 1 or 2 particle
                     # types. do not use second type (DM) as
                     # population
        self.pops = 2

        # Set number of terms for enclosedmass+tracer+anisotropy bins
        # = model parameters:
        self.nipol = 10   # IF CHANGED => set getnewdata = True to run
                         # data readout again

        self.nexp  = 3    # 3 more parameters for n(r) at 2, 4, and
                          # 8*rmax
        self.nepol = self.nipol + self.nexp + 3 # +3 means 1 more
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
        self.iscale = (self.nepol+0.5)/2 # scale below which range of
                                         # n(r)<2. instead of
                                         # maxrhoslope; is adapted in
                                         # gl_data.read_nu; if set to
                                         # -1 here, use maxrhoslope
                                         # everywhere
        self.rinfty = 30. # interpolate from last slope to slope at
                          # 10*max(xipol), where asymptote to \infty
                          # is reached, must be >= 11
        self.nbeta = 4   # number of parameters for beta, in sum of
                         # polynomials

        # next: # live points, > ndim, < 2^ndim, about number of
        # ellipsoids in phase space to be found
        self.ndim = self.nepol + (self.pops+1)*self.nepol + self.pops*self.nbeta
        self.nlive = 2*self.ndim


        ########## priors
        self.bprior   = False    # Baryon minimum surfden prior:
                                 # rho_*=sum_i nu_i only for mock!
        self.Mmodel = np.zeros(self.nipol) # for disc simulation case
        self.baryonmodel = 'sim' # read in surface density from
                                 # corresponding surfden file

        self.sprior  = False      # rising sig_z
        if self.geom=='sphere':   #  needed in disc case only
            self.sprior = False

        self.rhohalf = -1.        # prior density for rho at
                                  # half-light radius of tracers
                                  # calculated in gl_data
        self.rhospread = 1.       # with this spread, [dex] in log space
        self.nuspread = 0.5*np.ones(self.pops+1) # [dex] in log10 space

        self.nrtol  = 1.5/(8./self.nipol) # prior (max +/- range) for dn(r)/dlog(r); 8 is log(3000[pc])
        # norm1 = 17.**2 # offset of sig[0]/nu[0], from int starting
        # at zmin instead of 0 norm2 = 10.**2 # and for the second
        # component, if there is one
        self.quadratic = False    # linear or quad interpol.
        self.monotonic = False    # mono-prior on rho(z)
        self.monotonic_nu = True # and nu(z)
        self.adddarkdisc = False  # for disc mock case: add a dark disc?

        self.xipol = np.array([]) # [pc] hold radius bin centers
        self.xepol = np.array([]) # [pc] extended by 3 fudge bins

        self.maxR = 5.            # [Rscale], max range in radial bins
        self.maxrhoslope  = 4.    # maximum slope (change if
                                  # monotonicity prior used) of rho
        self.maxlognu = 3.        # direct sampling of nu: min value
        self.minlognu = -8.       # direct sampling of nu: max value
        self.maxbetaslope = 0.2   # linear (and 2nd..order) max slope
                                  # of beta*

        from gl_data import Datafile
        self.dat = Datafile()

        # Units: Rscale in [pc], surfdens_central (=Sig0) in
        # in [Munit/pc^2], and totmass
        # [Munit], and max(v_LOS) in [km/s]

        self.rscale=[];        self.nu0pc=[]
        self.Rscale=[];        self.Sig0pc=[]
        self.totmass=[];       self.maxvlos=[]
        
        

        # for investigations without data:
        if self.investigate != 'walk' and\
           self.investigate != 'triax' and\
           self.investigate != 'gaia' and\
           self.investigate != 'discmock':
            # each is set for all components and first component by
            # default
            self.rscale.append(1.);           self.rscale.append(1.)
            self.Rscale.append(1.);           self.Rscale.append(1.)
            self.nu0pc.append(1.);            self.nu0pc.append(1.)
            self.Sig0pc.append(1.);           self.Sig0pc.append(1.)
            self.totmass.append(1.);          self.totmass.append(1.)
            self.maxvlos.append(1.);          self.maxvlos.append(1.)

        import gl_class_files
        self.files = gl_class_files.Files(self, timestamp)
    ## \fn __init__(self, timestamp = '')
    # set up all parameters used in the course of the MultiNest run,
    # and the analysis routines in plot_multinest @param timestamp =
    # '', for output
    

    def __repr__(self):
        return "Params: "+self.files.dir
    ## \fn __repr__(self)
    # string representation for ipython

