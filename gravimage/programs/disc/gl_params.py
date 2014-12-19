#!/usr/bin/env python3

##
# @file
# all parameters for the gravimage MCMC, gaia investigation

# (c) GPL v3 2014 ETHZ Pascal S.P. Steger

import numpy as np
import ipdb
import gl_base as gb

def check_investigate(inv):
    if inv == 'discmock': return True
    if inv == 'discsim': return True
    raise Exception('wrong investigative case in gl_params')
    return False

class Params():
    def __init__(self, timestamp = '', investigate = '', case = -1):

        # Set machine and user variables
        # ----------------------------------------------------------------------
        self.machine = gb.get_machine()

        # Set investigation and geometry
        if investigate != '':
            self.investigate = investigate
        else:
            self.investigate  = 'discmock' # determine which data set to work on
                                  # 'discmock': set up simple model for disc
                                  # 'discsim': read in disc simulation
            self.geom = 'disc'

        check_investigate(self.investigate)
        self.case = 0 # used in spherical case
        self.pops = 2 # number of stellar tracer populations
                      # if changed: set getnewdata=True!

        # debug options
        # ----------------------------------------------------------------------
        self.checksig = False # debug sig calculation?
        self.debug = False # stop at wrong sanitazion?

        # Set number of tracer stars to look at take all particles #
        # case 0 want to set ntracer = 3e3 # case 1 ntracer = 1e4 #
        # case 2
        self.ntracer = [10000, 10000] # pop1, pop2, ..., pop_N (and
                                      # take sum for all tracers)

        # data options
        # ----------------------------------------------------------------------
        self.getnewdata = True  # get new data computed from
                                # observations before burn-in
        self.getnewpos  = True  # redo the first data conversion step
        self.binning = 'consttr' # 'linspace', 'logspace', 'consttr': binning of particles
        self.metalpop   = False # split metallicities with a separate
                                # MCMC
        self.maxR = 5.            # [Xscale], max range in radial bins


        # MultiNest options
        # ----------------------------------------------------------------------
        # Set number of terms for enclosedmass+tracer+anisotropy bins
        # = model parameters:
        self.chi2_nu_converged = False # first converge on Sig if set to False
        self.chi2_switch = 10. # if chi2<chi2_switch, add chi2_sig
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
        self.nlive = 100*self.ndim
        self.err = 1e300    # chi^2 for models which are impossible


        # parameter spaces
        # ----------------------------------------------------------------------
        self.rhohalf = 10**-1.    # prior density for rho at
                                  # half-light radius of tracers
                                  # calculated in gl_data
        self.log10rhospread = 2.       # with this spread, [dex] in log space
        self.rlimnr = -1 # scale below which range of
                         # n(r<rlimnr*r_half)<maxrhoslope/2, in multiples of r_half.
                         # calculated to values in [pc] in gl_data.read_nu;
                         # if set to -1 here, use maxrhoslope everywhere
        self.rlimnr_nu = -1 # same for nu, using same rhalf
        self.log10nuspread = 2.
        self.maxrhoslope  = 5.    # maximum slope (change if
                                  # monotonicity prior used) of rho
        self.maxnuslope = 5.
        # nztol: prior (max +/- range) for dn(r)/dlog(r)
        #   determine how far nr can wander with the max allowed nr slope
        #    from min(gp.xipol) to max(gp.xipol)
        self.nztol  = self.maxrhoslope
        self.nztol_nu = 2*self.maxnuslope # same for nu profile
        self.maxbetaslope = 0.2   # linear (and 2nd..order) max slope
                                  # of beta*
        self.beta00prior = False  # prior beta(r=0) = 0
        self.MtoLmin = 0.8
        self.MtoLmax = 3.


        # disc case
        # ----------------------------------------------------------------------
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

        # integration options
        # ----------------------------------------------------------------------
        self.usekappa   = False # switch to turn on (True) or off the
                                # calculation of kappa
        self.usezeta    = False  # switch to turn on (True) or off the
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
        # in [Munit/pc^2], and totmass_tracers
        # [Munit], and max(sigma_LOS) in [km/s]
        self.rscale=[];        self.nu0pc=[]
        self.Xscale=[];        self.Sig0pc=[]
        self.totmass_tracers=[];       self.maxsiglos=[]
        # for investigations without data:
        if self.investigate != 'discmock':
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
