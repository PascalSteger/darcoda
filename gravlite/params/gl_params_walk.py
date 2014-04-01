#!/usr/bin/env python3

##
# @file
# all parameters for the gravlite MCMC, gaia investigation

# (c) 2013 ETHZ Pascal S.P. Steger

import numpy as np
import os, pdb, logging

class Params():
    def __init__(self):
        self.loglevel = logging.WARNING
        # import sys
        # logging.basicConfig(stream=sys.stderr, level=self.loglevel)
        # self.LOG = logging.getLogger(__name__)
        self.machine = 'darkside'
        if not os.path.exists("/home/ast/read"):
            self.machine = 'local'    # machine: 'local' or 'darkside'

        self.investigate  = 'walk' # determine which data set to work on
                        # 'discmock': set up simple model for disc
                        # 'discsim': read in disc simulation
                        # 'checkdwarf': check int for analytic dwarf, sig_LOS
                        # 'hern': check simple Hernquist prof. from simwiki
                        # 'walk': check with full obs. cont. data from Walker
                        # 'gaia': 6D data (x,y,z,vx,vy,vz) from gaia challenge, 1 pop only
                        # 'fornax': real data from Fornax dwarf galaxy

        self.case = 0 # choose gaia models (1..8) or Walker (0..2,4,5) models, or triax models (0:core, 1:cusp)
        self.projcase = 4 # (1:X, 2:Y, 3:Z, 4:intermediate)

        import import_path as ip
        ip.set_geometry(self.investigate, self.machine)   # load spherical or disc version of the code
        
        # Set number of tracer stars to look at
        # take all particles                       # case 0
        # want to set ntracer = 3e3              # case 1
        #             ntracer = 1e4              # case 2
        self.ntracer = [10000, 10000] # pop1, pop2, ..., pop_N (and take sum for all tracers)

        self.getnewdata = True # get new data computed from observations before burn-in
        self.metalpop   = False # split metallicities with a separate MCMC
        self.lograd  = False # log steps for radial bin in readout, show x-axis in log scale
        self.consttr = True # set radial bin by constant number of tracer particles

        # units
        self.G1  = 6.67398e-11                # [m^3 kg^-1 s^-2]
        self.pc  = 3.08567758e16              # [m]
        self.msun= 1.981e30                   # [kg]
        self.km  = 1000.                      # [m]
        self.G1  = self.G1*self.msun/self.km**2/self.pc # [pc msun^-1 km^2 s^-2]

        self.geom = 'sphere'
        self.err = -1e30
        self.model = False # for walker mock data: plot model. set beta to Osipkov-Merrit profile
        if self.investigate != 'walk':
            self.model = False

        ########## integration options
        self.densint = False # use integral of dens to find mass, not binned sum
        self.even = 'avg' # for simps integration (everywhere): 'avg', 'first', 'last'
        self.usekappa   = False # switch to turn on (True) or off the calculation of kappa

        self.sim = 2 # use hernquist model with 1 or 2 particle types. do not use second type (DM) as population
        self.pops      = 1

        self.rstarhalf = 1000.

        # Set number of terms for enclosedmass+tracer+anisotropy bins = model parameters:
        self.nipol = 7   # IF CHANGED => set getnewdata = True to run data readout again
        self.nexp  = 3    # 3 more parameters for n(r) at 2, 4, and 8*rmax
        self.nepol = self.nipol + self.nexp + 3 # +3 means 1 more parameter for the density at half light radius,
                  #          1 more parameter for the asymptote to 0,
                  #          1 more parameter for the asymptote to \infty
        self.rinfty = 30. # interpolate from last slope to slope at 10*max(xipol), where asymptote to \infty is reached, must be >= 11
        self.nbeta = 2   # number of parameters for beta, in sum of polynomials

        # next: # live points, > ndim, < 2^ndim, about number of ellipsoids in phase space to be found
        self.nlive = 2*(self.nepol + self.pops*self.nepol + self.pops*self.nbeta)
        # if more than nlive points trigger the same penalty value => Multinest finishes

        ########## priors
        self.bprior   = True    # Baryon minimum surfden prior: rho_*=sum_i nu_i only for mock!
        self.Mmodel = np.zeros(self.nipol)
        self.baryonmodel = 'sim' # read in surface density from corresponding surfden file

        self.delta0 = np.zeros(self.nipol)
        # sigmaprior, default: -1
        self.sigmaprior1 = 0.3
        self.sigmaprior2 = 0.3

        self.sprior  = False            # rising sig_z needed in disc case
        if self.geom=='sphere':
            self.sprior = False
        self.rprior  = True    # regularize Nuz
        self.rhotol  = 1.0
        # norm1   = 17.**2 # offset of sigma[0]/nu[0], from int starting at zmin instead of 0
        # norm2   = 10.**2 # and for the second component, if there is one
        self.quadratic = False           # linear or quad interpol. 
        self.monotonic = False           # mono-prior on nu(z)

        self.xipol = np.array([])
        self.xepol = np.array([])

        self.maxrhoslope  = 6.
        self.maxnuslope   = 6.
        self.maxbetaslope = 3. # linear (and 2nd, 3rd...) max slope of beta*

        self.nuerrcorr = 1. # 2.   # scale error from grw_dens by this amount
        self.sigerrcorr = 1.  # 1.5  # best done directly in the corresponding grw_dens, grw_siglos
                     # thus, both times 1. is best
        self.kaperrcorr = 1.

        ########## global variables, not set to value here
        from gl_data import Datafile
        self.dat = Datafile()
        self.dof = 1

        # Units: 
        # rscale in [pc], surfdens_central (=Nu0) in [munit/rscale**2],
        # and in [munit/pc**2], and totmass [munit], and max(v_LOS) in [km/s]

        self.rscale=[]
        self.nu0rscale=[]
        self.nu0pc=[]
        self.totmass=[]
        self.maxvlos=[] # unit system
        self.Rscale=[]
        self.Nu0rscale=[]
        self.Nu0pc=[]
        if self.investigate != 'walk' and self.investigate != 'triax' and self.investigate != 'gaia':
            # each is set for all components and first component by default
            self.rscale.append(1.)
            self.Rscale.append(1.)
            self.rscale.append(1.)
            self.Rscale.append(1.)
            self.nu0rscale.append(1.)
            self.Nu0rscale.append(1.)
            self.nu0rscale.append(1.)
            self.Nu0rscale.append(1.)
            self.nu0pc.append(1.)
            self.Nu0pc.append(1.)
            self.nu0pc.append(1.)
            self.Nu0pc.append(1.)
            self.totmass.append(1.)
            self.totmass.append(1.)
            self.maxvlos.append(1.)
            self.maxvlos.append(1.)

        import gl_class_files
        self.files = gl_class_files.Files(self)

    def len(self):
        return 1
