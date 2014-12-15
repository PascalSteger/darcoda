#!/usr/bin/env ipython3

##
# @file
# all parameters for the gravimage MCMC, gaia investigation

# (c) 2013 ETHZ Pascal S.P. Steger

import numpy as np
import pdb
import socket
import getpass

def check_investigate(inv):
    if inv == 'discmock': return True
    if inv == 'discsim': return True
    raise Exception('wrong investigative case in gl_params')
    return False

class Params():
    def __init__(self, timestamp = '', investigate = '', case = -1):

        # Set machine and user variables
        # ----------------------------------------------------------------------
        host_name = socket.gethostname()
        user_name = getpass.getuser()
        if 'darkside' in host_name:
            self.machine = 'darkside'
        elif 'pstgnt332' in host_name:
            self.machine = 'pstgnt332'
        elif ('lisa' in host_name) and ('login' in host_name) and ('hsilverw' in user_name):
            self.machine = 'lisa_HS_login'
        elif ('lisa' in host_name) and ('login' not in host_name) and ('hsilverw' in user_name):
            self.machine = 'lisa_HS_batch'
        elif ('lisa' in host_name) and ('login' in host_name) and ('sofia' in user_name):
            self.machine = 'lisa_SS_login'
        elif ('lisa' in host_name) and ('login' not in host_name) and ('sofia' in user_name):
            self.machine = 'lisa_SS_batch'

        # Set investigation and geometry
        if investigate != '':
            self.investigate = investigate
        else:
            self.investigate  = 'simplenu' # determine which data set to work on
                                  # 'discmock': set up simple model for disc
                                  # 'discsim': read in disc simulation
        self.geom = 'disc'

        #check_investigate(self.investigate)
        self.case = 0 # used in spherical case
        self.ntracer_pops = 1 # number of stellar tracer populations
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

        # ----------------------------------------------------------------------
        #HS Working Line
        #Everthing above this line is old, and hasn't been considered for keeping
        #deletion, or modification
        # ----------------------------------------------------------------------

        self.nbins=15 # Number of bins to split tracer stars into
        self.nrhonu = self.nbins + 2 # Number of points where rho and nu parameters will be set,
                                   # e.g. bin centres, plus zC=0, and zLS (position
                                   # of Last Star)
        self.nbaryon_pops = 0 # Number of baryon populations to look at
                               # =0 if doing simple mass model (eg DM profile describes
                               # all mmass)
        self.nbaryon_params = 0 # Number of parameters to describe baryon population
                                 #  Holmberg & Flynn = 15
                                 #  with baryon observational information = nrho
        self.ndim = 1 + 2*(self.nrhonu + 1) + self.nbaryon_pops*self.nbaryon_params
            # Constant C from sigma_z calculation, nrho + 1 params for rho (nrho
            # points for kz_rho, plus central density of rho, eg rho_C), similarly
            # nrho +1  params for nu, plus the number of params for all baryon pops

        #Limits for central densities (z=0)
#        self.rho_C_max = 0.2 #Msun pc^-3, for either DM or baryons (cf rho_b = 0.0914 Msun pc^-3, Flynn+ 2006)
        self.rho_C_max = 0.2*1.E9  #Msun kpc^-3
        self.rho_C_min = 0.0 #Msun pc^-3
        self.nu_C_max = 0.0 # no. stars pc^-3, full value calculated in external_data
        self.nu_C_min = 10.0 # no. stars pc^-3

        #Limits for central kz values (z=0)
        self.kz_rho_C_max = 10.0
#        self.kz_rho_C_min = 0.0
        self.kz_rho_C_min = -1.0 #SS
        self.kz_nu_C_max = 10.0
#        self.kz_nu_C_min = 0.0
        self.kz_nu_C_min = -1.0 #SS

        #Maximum kz_slope (=dk/dz)
        self.max_kz_slope = 5.0

        #Limits for sigz central value
        self.sigz_C_max = 30.
        self.sigz_C_min = 10.







        #HS Working Line
        #Everthing below this line is old, and hasn't been considered for keeping
        #deletion, or modification
        # ----------------------------------------------------------------------


        # MultiNest options
        # ----------------------------------------------------------------------
        # Set number of terms for enclosedmass+tracer+anisotropy bins
        # = model parameters:
        self.chi2_nu_converged = False # first converge on Sig if set to False
        self.chi2_switch = 100. # if chi2*10 < chi2_switch, add chi2_sig
        self.chi2_switch_mincount = 250. # demand that this number of profiles with
                                        # chi2<chi2_switch are found before adding chi2_sig
        self.chi2_switch_counter = 0. # start the counter at 0

        self.nipol = 15   # IF CHANGED => set getnewdata = True to run
                         # data readout again
        self.nexp  = 3    # more fudge parameters at r<rmin and r>rmax
        self.nepol = self.nipol + 2*self.nexp     # number of parameters for
                                                # direct mapping of nu(r)

        self.nfine = 30  # number of entries in integral lookup table
                         # gives no. log spaced points
        self.rinfty = 5. # interpolate from last slope to slope at
                          # 10*max(xipol), where asymptote to \infty
                          # is reached, must be >= 11
        self.nbeta = 0   # number of parameters for beta, in sum ofgre
                         # polynomials
        # TODO: if not using rhostar, subtract gp.nrho:

        # +1 for norm, rho, rhostar, MtoL, [nu_i, tilt_i]
        # live points, > ndim, < 2^ndim, about number of
        # ellipsoids in phase space to be found
        self.nlive = 100*self.ndim
        self.err = 1e300    # chi^2 for models which are impossible

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
        self.z_bincenters = np.array([]) # [pc] holds the bin centers, H Silverwood 21/11/14
        self.z_binmins = np.array([])
        self.z_binmaxs = np.array([])
        self.z_all_pts = np.array([]) # [pc] holds [zC = 0, z_bin_centers, zLS]


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
