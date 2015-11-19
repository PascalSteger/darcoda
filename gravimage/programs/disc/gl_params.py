#!/usr/bin/env ipython3

##
# @file
# all parameters for the gravimage MCMC, gaia investigation

# (c) 2013 ETHZ Pascal S.P. Steger

import numpy as np
import pdb
import gl_helper as gh
from mpi4py import MPI
#import socket
#import getpass

def check_investigate(inv):
    if inv == 'discmock': return True
    if inv == 'discsim': return True
    raise Exception('wrong investigative case in gl_params')
    return False

class Params():
    def __init__(self, timestamp = '', investigate = '', case = -1):

        myrank = MPI.COMM_WORLD.Get_rank()
        nprocs = MPI.COMM_WORLD.Get_size()
        procnm = MPI.Get_processor_name()

        print('P', myrank, ': Initiating Params()')

        # Set investigation and geometry
        # ----------------------------------------------------------------------
        if investigate != '':
            self.investigate = investigate
        else:
            self.investigate  = 'simplenu' # determine which data set to work on
                                  # 'discmock': set up simple model for disc
                                  # 'discsim': read in disc simulation
        self.geom = 'disc'

        #check_investigate(self.investigate)
        self.case = case # used in spherical case

        # debug options
        # ----------------------------------------------------------------------
        self.checksig = False # debug sig calculation?
        self.debug = False # stop at wrong sanitazion?


        # data and analysis options
        # ----------------------------------------------------------------------
        self.getnewdata = True  # get new data computed from
                                # observations before burn-in
        self.getnewpos  = True  # redo the first data conversion step

        if self.investigate == 'simplenu':
            self.external_data_file, self.external_data_file_tilt = gh.ext_file_selector_simplenu([2], '1e6', '', False)
                #ext_file_selector_simplenu(pops, sampling, darkdisk, tilt)
                #pops = [x,x], eg which population to use
                #sampling = 1e4, 1e5, 1e6,
                #darkdisk= '', 'dd', 'bdd'
                #tilt = True or False

        elif self.investigate == 'disc_nbody' and self.case == 1: ######## Hunt & Kawata Data sets ########
            self.external_data_file = ['/Hunt_Kawata_GCD/HK_cut_1kpc_GalC_z_vz.dat']
            self.external_data_file_tilt = ['/Hunt_Kawata_GCD/HK_cut_1kpc_GalC_z_vz.dat']

        elif self.investigate == 'disc_nbody' and self.case == 2: ######## Garbari et al. Wedges ########
            Gar_Root = '/Garbari_Nbody/wedges/'
            Gar_Simulation = 'mwhr' #'unevolved', 'mwhr', 'LLMC10', 'LLMC60'
            Gar_Radius = 'r8500' #'r8500', 'r10500', 'r12500'
            Gar_Angle = 'ang0' #'ang45', 'ang90', 'ang135', 'ang180', 'ang225', 'ang270', 'ang315'
            self.external_data_file = [Gar_Root + Gar_Simulation + '_' + Gar_Radius + '_' + Gar_Angle + '_stars_z_vz_vR.dat']

        #Check case agains data file
        if self.investigate == 'simplenu' and 'simplenu' not in self.external_data_file[0]:
                print('Investigate = ', self.investigate, ', Case = ', self.case, ', data files = ', self.external_data_file[0])
                raise Exception('Case and data file mismatch')
        if self.investigate == 'disc_nbody' and self.case == 1 and 'Hunt_Kawata_GCD' not in self.external_data_file[0]:
                print('Investigate = ', self.investigate, ', Case = ', self.case, ', data files = ', self.external_data_file[0])
                raise Exception('Case and data file mismatch')
        if self.investigate == 'disc_nbody' and self.case == 2 and 'Garbari_Nbody' not in self.external_data_file[0]:
                print('Investigate = ', self.investigate, ', Case = ', self.case, ', data files = ', self.external_data_file[0])
                raise Exception('Case and data file mismatch')

        #Count number of tracer populations
        self.ntracer_pops = len(self.external_data_file) # number of stellar tracer populations

        #self.data_z_cut = 1.2  # [kpz] only use (& bin) data up to this z limit
        self.data_z_cut = [2.4]  # (set > data z_max to use all avaiable data)
        if len(self.data_z_cut) != self.ntracer_pops:
            raise Exception('Incorrect data_z_cut vector length')

        self.tilt = False   # If also modelling the tilt

        self.darkmattermodel = 'const_dm' # const_dm = const DM dens in z
        #self.darkmattermodel = 'kz_dm'  # kz_dm = kz parameterization of DM CURRENTLY NOT WORKING
        #self.darkmattermodel = 'ConstPlusDD' # constant DM + DM disc component

        #self.dd_data = False # if we are to plot dd analytics or not

        self.binning = 'consttr' # 'linspace', 'logspace', 'consttr': binning of particles
        #self.binning = 'linspace' # 'linspace', 'logspace', 'consttr': binning of particles
        self.nbins = [20]   # Number of bins to split each population of tracer stars into
        if len(self.nbins) != self.ntracer_pops:
            raise Exception('Number of bins vector is of incorrect length')
        self.nrhonu = sum(self.nbins) + 1 # Number of points where rho and nu parameters will be set,
                                   # e.g. bin centres, plus zC=0

        #Tracer Density description
        #self.nu_model = 'kz_nu'
        self.nu_model = 'gaussian_data'

        #------------------------
        self.TheoryData = False   # If true using theoretical bin values as indata
        self.hyperparams = False  # Use hyperparameters, 2 params, range (0.1->10)*meanerr
        #Dark matter options
        self.adddarkdisc = False  # for disc mock case: add a dark disc?
        # self.adddarkdisc is currently not used !
        self.param_headers = []



        #Baryon options
        self.baryonmodel = 'simplenu_baryon' #set baryon model
                                             # none = all mass in DM
                                    # setting baryonmodel == 'none' not possible without major rewrite
                                    # simplenu_baryon = model used to generate simplenu mock data
        self.nbaryon_pops = 1 # Number of baryon populations to look at
                                    # =0 if doing simple mass model (eg DM profile describes
                                    # all mmass)
        self.nbaryon_params = 2 # Number of parameters to describe baryon population
                                    #  simplenu_baryon = 2
                                    #  Holmberg & Flynn = ?
                                    #  with baryon observational information = nrho
        self.ntilt_params = 3  # Number of parameters used to describe tilt

        self.scan_rhonu_space = False #Search directly in rho or nu space, i.e. no kr parametrization
        if self.scan_rhonu_space:
            self.nrhonu = self.nbins #Param count will be rho_C + rho_bins

        # Parameters escribing the underlying baryonic model and the mock data
        # This is assuming baryonmodel = simplenu_baryon above, rewrite required if this is changed.

        #Total dimensions count
        if self.darkmattermodel == 'kz_dm':
            self.ndim = self.ntracer_pops + (self.nrhonu + 1) + (sum(self.nbins) + self.ntracer_pops) + self.nbaryon_pops*self.nbaryon_params
                # Constant C from sigma_z calculations (one for each tracer population)
                # DM kz_rho defined at all bin centres and zC=0, plus rho_C
                # kz_nu definitions at bin centres, plus at zC = 0 for each population
                # baryon parameters

        elif self.darkmattermodel == 'const_dm':
            self.ndim = self.ntracer_pops + 1 + (sum(self.nbins) + 2*self.ntracer_pops) + self.nbaryon_pops*self.nbaryon_params
                # Constant C from sigma_z calculations (one for each tracer population)
                # Constant DM density
                # kz_nu definitions at bin centres, plus at zC = 0 for each population
                # baryon parameters

        elif self.darkmattermodel == 'ConstPlusDD':
            self.ndim = self.ntracer_pops + 3 + (sum(self.nbins) + 2*self.ntracer_pops) + self.nbaryon_pops*self.nbaryon_params
                # Constant C from sigma_z calculations (one for each tracer population)
                # Constant DM density, plus two params for a dark disc
                # kz_nu definitions at bin centres, plus at zC = 0, and nu_C, for each population
                # baryon parameters

        if self.nu_model == 'gaussian_data':
            self.ndim -= self.ntracer_pops

        if self.hyperparams:
            self.ndim += 2

        if self.tilt:
            self.ndim += self.ntilt_params * self.ntracer_pops

        self.z_err_measurement = 0.05 # Measurement error on z, fraction, eg 0.05 = 5%
        self.vz_SDerr_meas = 5.  # Measurement error on vz, [km s^-1]
        self.mc_err_N_iters = int(100) #Number of iterations to perform when doing MC error estimation


        # Priors
        # ----------------------------------------------------------------------
        # Limits for central densities (z=0)
        self.rho_C_max = 1.0E8  #Msun kpc^-3 for either DM or baryons (cf rho_b = 0.0914 Msun pc^-3, Flynn+ 2006)
        self.rho_C_min = 1.0E6 #Msun pc^-3
        self.rho_C_prior_type = 'log' #log, linear, gaussian
        self.nu_C_max = 0.0 # no. stars pc^-3, full value calculated in external_data
        self.nu_C_min = 0.0 # no. stars pc^-3
        self.nu_C_prior_type = 'log'

        # Limits for central kz values (z=0)
        self.kz_rho_C_max = 5.0
        self.kz_rho_C_min = -5.0 #SS
        self.kz_nu_C_max = 5.#20.0
        self.kz_nu_C_min = -5.0 #SS

        # Maximum kz_slope (=dk/dz)
        #self.max_kz_slope = 10.0
        self.max_kz_slope = 5.0

        # Limits for sigz central value
        self.sigz_C_max = 0. #set from data in gr_external_data
        self.sigz_C_min = 0.

        # Monotonicity priors
        self.monotonic_rho = False    # mono-prior on rho(z)
        self.monotonic_nu = False # mono-prior on nu(z)

        # kz selection scheme
        self.kz_rho_selection = 'gaussian'
        self.kz_nu_selection = 'tophat'

        # Log or linear priors for rhonu scanning
        self.prior_type_rho = 'gaussian' # 'log' or 'linear'
        self.prior_type_nu = 'gaussian' # 'log' or 'linear'

        # Simplenu Baryon model priors

        #self.simplenu_baryon_K_max = 1650  #JR model has K = 1500.
        #self.simplenu_baryon_K_min = 1350
        #self.simplenu_baryon_D_max = 0.2  #JR model has D = 0.18
        #self.simplenu_baryon_D_min = 0.16

        self.simplenu_baryon_K_max = 2000  #1600 #1700. #JR model has K = 1500.
        self.simplenu_baryon_K_min = 1000  #1400 #1300.
        self.simplenu_baryon_D_max = 0.24  #0.20 #0.5 #JR model has D = 0.18
        self.simplenu_baryon_D_min = 0.12  #0.16 #0.05

        # Simplenu DM disc model priors
        self.simplenu_dm_K_max = 1500.  #JR model has K = 300.
        self.simplenu_dm_K_min = 0.
        self.simplenu_dm_D_max = 3.5  #JR model has D = 2.5
        self.simplenu_dm_D_min = 1.5

        # Tilt priors
        #self.tilt_A_max = -0.005  # simple2 mock has A = -0.0087
        #self.tilt_A_min = -0.012
        #self.tilt_n_max = 1.9    # simple2 mock has n = 1.44
        #self.tilt_n_min = 1.
        #self.tilt_R_max = 3.5     # simple2 mock has R = 2.5
        #self.tilt_R_min = 1.5

        #self.tilt_A_max = 0.0  # simple2 mock has A = -8.7 -0.0087
        #self.tilt_A_min = -1.0
        self.tilt_A_max = 0.0  # simple2 mock has A = -181.77 for z[kpc]
        self.tilt_A_min = -200.00

        self.tilt_n_max = 3.0    # simple2 mock has n = 1.44
        self.tilt_n_min = 0.1
        self.tilt_R_max = 4.5     # simple2 mock has R = 2.5
        self.tilt_R_min = 0.5


        # MultiNest options
        # ----------------------------------------------------------------------
        self.map_priors = False

        # live points, > ndim, < 2^ndim, about number of
        # ellipsoids in phase space to be found
        self.nlive = 100*self.ndim

        if self.map_priors:
            self.nlive = int(1E4)*self.ndim

        self.err = 1e300    # chi^2 for models which are impossible

        #fraction of profiles to save, set <0 for no profile saving
        self.save_fraction = -1.0

        #Plotting flag
        self.plotting_flag = False

        # filesystem-related
        # ----------------------------------------------------------------------
        import import_path as ip
        ip.set_geometry(self.geom) # load spherical or
                                    # disc version
                                    # of the code

        print('P', myrank, ': Instantiating Files class')
        import gl_class_files
        self.files = gl_class_files.Files(self, timestamp)
        from gl_data import Datafile
        self.dat = Datafile()

        # global arrays
        # ----------------------------------------------------------------------
        self.z_bincenter_vecs = [] # [(np.zeros(self.nbins[ii])) for ii in range(0, self.ntracer_pops)]
        self.z_binmin_vecs = [] #[(np.zeros(self.nbins[ii])) for ii in range(0, self.ntracer_pops)]
        self.z_binmax_vecs = [] #[(np.zeros(self.nbins[ii])) for ii in range(0, self.ntracer_pops)]
        self.z_all_pts_unsort = np.array([]) #Simple concatenation of bincentres and zC = 0
        self.z_all_pts_sorted = [] #Sorted list of bincentres and zC = 0
        self.z_vec_masks = [] #[[].append(None) for ii in range(0, self.ntracer_pops)]

        # Global constants
        self.Rsun = 8.  # [kpc]  Sun's distance to galactic center

    ## \fn __init__(self, timestamp = '')
    # set up all parameters used in the course of the MultiNest run,
    # and the analysis routines in plot_multinest @param timestamp =
    # '', for output


    def __repr__(self):
        return "Params: "+self.files.dir
    ## \fn __repr__(self)
    # string representation for ipython
