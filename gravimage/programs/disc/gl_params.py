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

        # Set machine and user variables
        # ----------------------------------------------------------------------
        self.machine, dummy = gh.detect_machine()

        # Set global constants
        #self.Msun = 1.9891e30 # Solar mass in kg
        #self.pc = 3.0857e16   # pc in m 
        #self.Mp = 1.6726231e-27 # proton mass in kg


        # Set investigation and geometry
        # ----------------------------------------------------------------------
        if investigate != '':
            self.investigate = investigate
        else:
            self.investigate  = 'simplenu' # determine which data set to work on
            #self.investigate  = 'obsbary' # determine which data set to work on
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


        # data and analysis options
        # ----------------------------------------------------------------------
        self.getnewdata = True  # get new data computed from
                                # observations before burn-in
        self.getnewpos  = True  # redo the first data conversion step

        ######## New mock data: obsbary ##########
        #self.external_data_file = ['/obsbary/obs_32e3_raw.dat']
        #self.external_data_file = ['/simplenu/obs_32e3_raw.dat']
        self.external_data_file = ['/simplenu/obs_32e3_tilt_raw.dat']
        self.external_data_file_tilt = ['/simplenu/obs_32e3_sigRz.dat']
        #self.external_data_file = ['/obsbary/simple2_1e4nu_sigz_raw.dat']

        ######## 1E4 Sample Size ########
        ##Popn 1 1E4
        #self.external_data_file= ['/simplenu/simple_1e4nu_sigz_raw.dat']
        ###Popn 2 1E4
        #self.external_data_file= ['/simplenu/simple2_1e4nu_sigz_raw.dat']
        ##Popn 1 & 2 1E4
        #self.external_data_file= ['/simplenu/simple_1e4nu_sigz_raw.dat', '/simplenu/simple2_1e4nu_sigz_raw.dat']

        ##Popn 1 1E4 DD
        #self.external_data_file= ['/simplenu/simple_dd_1e4nu_sigz_raw.dat']
        ##Popn 2 1E4 DD
        #self.external_data_file= ['/simplenu/simple2_dd_1e4nu_sigz_raw.dat']
        ##Popn 1 & 2 1E4 DD
        #self.external_data_file= ['/simplenu/simple_dd_1e4nu_sigz_raw.dat', '/simplenu/simple2_dd_1e4nu_sigz_raw.dat']

        ##Popn 2 1E4 Big DD
        #self.external_data_file= ['/simplenu/simple2_bdd_1e4nu_sigz_raw.dat']

        ##Popn 2 1E4 Big DD with Tilt
        #self.external_data_file= ['/simplenu/simple2_bdd_tilt_1e4nu_sigz_bin.dat']
        #self.external_data_file_tilt= ['/simplenu/simple2_bdd_tilt_1e4nu_sigRz_raw.dat']


        ######## 1E5 Sample Size ########
        ##Popn 1 1E5
        #self.external_data_file= ['/simplenu/simple_1e5nu_sigz_raw.dat']
        ##Popn 2 1E5
        #self.external_data_file= ['/simplenu/simple2_1e5nu_sigz_raw.dat']
        ##Popn 1 & 2 1E5
        #self.external_data_file= ['/simplenu/simple_1e5nu_sigz_raw.dat', '/simplenu/simple2_1e5nu_sigz_raw.dat']

        ##Popn 1 1E5 DD
        #self.external_data_file= ['/simplenu/simple_dd_1e5nu_sigz_raw.dat']
        ##Popn 2 1E5 DD
        #self.external_data_file= ['/simplenu/simple2_dd_1e5nu_sigz_raw.dat']
        ##Popn 1 & 2 1E5 DD
        #self.external_data_file= ['/simplenu/simple_dd_1e5nu_sigz_raw.dat', '/simplenu/simple2_dd_1e5nu_sigz_raw.dat']

        ##Popn 2 1E5 Big DD
        #self.external_data_file= ['/simplenu/simple2_bdd_1e5nu_sigz_raw.dat']

        ##Popn 2 1E5 Big DD with Tilt
        #self.external_data_file= ['/simplenu/simple2_bdd_tilt_1e5nu_sigz_bin.dat']
        #self.external_data_file_tilt= ['/simplenu/simple2_bdd_tilt_1e5nu_sigRz_raw.dat']


        ####### 1E6 Sample Size ########
        ##Popn 1 1E6
        #self.external_data_file= ['/simplenu/simple_1e6nu_sigz_raw.dat']
        #Popn 2 1E6
        #self.external_data_file= ['/simplenu/simple2_1e6nu_sigz_raw.dat']
        ##Popn 1 & 2 1E6
        #self.external_data_file= ['/simplenu/simple_1e6nu_sigz_raw.dat', '/simplenu/simple2_1e6nu_sigz_raw.dat']

        ##Popn 1 1E6 DD
        #self.external_data_file= ['/simplenu/simple_dd_1e6nu_sigz_raw.dat']
        ##Popn 2 1E6 DD
        #self.external_data_file= ['/simplenu/simple2_dd_1e6nu_sigz_raw.dat']
        ##Popn 1 & 2 1E6 DD
        #self.external_data_file= ['/simplenu/simple_dd_1e6nu_sigz_raw.dat', '/simplenu/simple2_dd_1e6nu_sigz_raw.dat']

        ##Popn 2 1E6 Big DD
        #self.external_data_file= ['/simplenu/simple2_dd_1e6nu_sigz_raw.dat']

        ##Popn 1 1E6 Tilt
        #self.external_data_file= ['/simplenu/simple_tilt_1e6nu_sigz_raw.dat']
        #self.external_data_file_tilt= ['/simplenu/simple_tilt_1e6nu_sigRz_raw.dat'] #MISSING
        ##Popn 2 1E6 Tilt
        #self.external_data_file= ['/simplenu/simple2_tilt_1e6nu_sigz_raw.dat']
        #self.external_data_file_tilt= ['/simplenu/simple2_tilt_1e6nu_sigRz_raw.dat']
        ##Popn 1 & 2 1E6 Tilt MISSING ELEMENTS
        #self.external_data_file= ['/simplenu/simple_tilt_1e6nu_sigz_raw.dat', '/simplenu/simple2_tilt_1e6nu_sigz_raw.dat']
        #self.external_data_file= ['/simplenu/simple_tilt_1e6nu_sigRz_raw.dat', '/simplenu/simple2_tilt_1e6nu_sigRz_raw.dat'] #MISSING 0

        ##Popn 2 1E6 Big DD and Tilt
        #self.external_data_file= ['/simplenu/simple2_bdd_tilt_1e6nu_sigz_raw.dat']
        #self.external_data_file_tilt= ['/simplenu/simple_bdd_tilt_1e6nu_sigRz_raw.dat']

        ######## Hunt & Kawata Data sets ########
        #self.external_data_file = ['/Hunt_Kawata_GCD/HK_cut_1kpc_z_vz.dat']
        #self.external_data_file_tilt = ['/Hunt_Kawata_GCD/HK_cut_1kpc_z_vz.dat']



        self.dd_data = False # if we are to plot dd analytics or not

        #self.data_z_cut = 1.2  # [kpz] only use (& bin) data up to this z limit
        self.data_z_cut = [2.4]  # (set > data z_max to use all avaiable data)

        self.tilt = True   # If also modelling the tilt

        self.darkmattermodel = 'const_dm' # const_dm = const DM dens in z
        #self.darkmattermodel = 'kz_dm'  # kz_dm = kz parameterization of DM
        #self.darkmattermodel = 'ConstPlusDD' # constant DM + DM disc component

        self.binning = 'consttr' # 'linspace', 'logspace', 'consttr': binning of particles
        #self.binning = 'linspace' # 'linspace', 'logspace', 'consttr': binning of particles
        self.nbins = [10]   # Number of bins to split each population of tracer stars into
        self.nrhonu = sum(self.nbins) + 1 # Number of points where rho and nu parameters will be set,
                                   # e.g. bin centres, plus zC=0

        #------------------------
        self.fit_outer_z_half = False  # Use only the outer half of the bins for chi2 fit
        self.TheoryData = False   # If true using theoretical bin values as indata
        self.hyperparams = False  # Use hyperparameters, 2 params, range (0.1->10)*meanerr
        #Dark matter options
        self.adddarkdisc = False  # for disc mock case: add a dark disc?
        # self.adddarkdisc is currently not used !

        #Baryon options
        #self.baryonmodel = 'simplenu_baryon' #set baryon model
        self.baryonmodel = 'obs_baryon' #more complex baryon model based on observations 
                                             # none = all mass in DM
                                    # setting baryonmodel == 'none' not possible without major rewrite
                                    # simplenu_baryon = model used to generate simplenu mock data
        self.nbaryon_pops = 1 # Number of baryon populations to look at
                                    # =0 if doing simple mass model (eg DM profile describes
                                    # all mmass)
        if self.baryonmodel == 'obs_baryon':
                    #  obs_baryon: gas: 5*2, stars: 5*2 + 4 => total 24
            self.nbaryon_params = 24 # Number of parameters to describe bary pop
        else:   # simplenu baryons
            self.nbaryon_params = 2
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
        self.rho_C_max = 1.0E10  #Msun kpc^-3 for either DM or baryons (cf rho_b = 0.0914 Msun pc^-3, Flynn+ 2006)
        self.rho_C_min = 1.0E4 #Msun pc^-3  (Lower end od DM density prior)
        self.rho_C_prior_type = 'log' #log, linear, gaussian
        self.nu_C_max = 0.0 # no. stars pc^-3, full value calculated in external_data
        self.nu_C_min = 10.0 # no. stars pc^-3
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
        self.sigz_C_max = 500. #was 5 and 50
        self.sigz_C_min = 0.5

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
        self.simplenu_baryon_K_max = 1650  #JR model has K = 1500.
        self.simplenu_baryon_K_min = 1350
        self.simplenu_baryon_D_max = 0.2  #JR model has D = 0.18
        self.simplenu_baryon_D_min = 0.16
        #self.simplenu_baryon_K_max = 2000  #1600 #1700. #JR model has K = 1500.
        #self.simplenu_baryon_K_min = 1000  #1400 #1300.
        #self.simplenu_baryon_D_max = 0.24  #0.20 #0.5 #JR model has D = 0.18
        #self.simplenu_baryon_D_min = 0.12  #0.16 #0.05

        # Observational baryon priors
        self.obs_bary_gasH2_Sigma = 1.0   # Surface density of H2 gas in Msun/pc2
        self.obs_bary_gasH2_Sigma_err = 0.3  # Estimated error on the above surface density
        self.obs_bary_gasH2_h = 105.  # H2 scale height in pc
        self.obs_bary_gasH2_h_err = 0.2*self.obs_bary_gasH2_h  # Estimated error on above h. UNKNOWN

        self.obs_bary_gasHIcnm_Sigma = 6.21   # Surface density of H2 gas in Msun/pc2
        self.obs_bary_gasHIcnm_Sigma_err = 0.15*self.obs_bary_gasHIcnm_Sigma
        self.obs_bary_gasHIcnm_h = 127.  # H2 scale height in pc
        self.obs_bary_gasHIcnm_h_err =  0.2*self.obs_bary_gasHIcnm_h 

        self.obs_bary_gasHIwnm1_Sigma = 2.51   # Surface density of H2 gas in Msun/pc2
        self.obs_bary_gasHIwnm1_Sigma_err = 0.15*self.obs_bary_gasHIwnm1_Sigma 
        self.obs_bary_gasHIwnm1_h = 318.  # H2 scale height in pc
        self.obs_bary_gasHIwnm1_h_err = 0.2*self.obs_bary_gasHIwnm1_h

        self.obs_bary_gasHIwnm2_Sigma = 2.14   # Surface density of H2 gas in Msun/pc2
        self.obs_bary_gasHIwnm2_Sigma_err = 0.15*self.obs_bary_gasHIwnm2_Sigma 
        self.obs_bary_gasHIwnm2_h = 403.  # H2 scale height in pc
        self.obs_bary_gasHIwnm2_h_err = 0.2*self.obs_bary_gasHIwnm2_h

        self.obs_bary_gasHII_Sigma = 1.8   # Surface density of H2 gas in Msun/pc2
        self.obs_bary_gasHII_Sigma_err = 0.1
        self.obs_bary_gasHII_h = 1590.  # H2 scale height in pc
        self.obs_bary_gasHII_h_err = 0.2*self.obs_bary_gasHII_h 
        
        # STARS
        self.obs_bary_MS3_Sigma = 0.5   # MS stars up to magnitude 3
        self.obs_bary_MS3_Sigma_err = 0.1*self.obs_bary_MS3_Sigma
        self.obs_bary_MS3_h = 140.
        self.obs_bary_MS3_h_err = 0.2*self.obs_bary_MS3_h

        self.obs_bary_MS4_Sigma = 0.8   # MS stars (from 3) up to magnitude 4
        self.obs_bary_MS4_Sigma_err = 0.1*self.obs_bary_MS4_Sigma
        self.obs_bary_MS4_h = 236.
        self.obs_bary_MS4_h_err = 0.2*self.obs_bary_MS4_h

        self.obs_bary_MS5_Sigma = 2.2+0.4  # including giant stars
        self.obs_bary_MS5_Sigma_err = 0.1*self.obs_bary_MS5_Sigma
        self.obs_bary_MS5_h = 384.
        self.obs_bary_MS5_h_err = 0.2*self.obs_bary_MS5_h

        self.obs_bary_MS8_Sigma = 5.8
        self.obs_bary_MS8_Sigma_err = 0.1*self.obs_bary_MS8_Sigma
        self.obs_bary_MS8_h = 400.
        self.obs_bary_MS8_h_err = 0.2*self.obs_bary_MS8_h

            # Thick disk fraction of the dimmer MS stars: MS5 and MS8
        self.obs_bary_MS_thick_fraction = 0.1
        self.obs_bary_MS_thick_fraction_err = 0.5*self.obs_bary_MS_thick_fraction
        self.obs_bary_MS_thick_h = 1000.
        self.obs_bary_MS_thick_h_err = 0.2*self.obs_bary_MS_thick_h

        self.obs_bary_dwarfs_Sigma = 17.3+1.2+4.9+0.3  # M, brown and white dwarfs
        self.obs_bary_dwarfs_Sigma_err = 0.1*self.obs_bary_dwarfs_Sigma
        self.obs_bary_dwarfs_h1 = 332.
        self.obs_bary_dwarfs_h1_err = 0.2*self.obs_bary_dwarfs_h1
        self.obs_bary_dwarfs_h2 = 609.
        self.obs_bary_dwarfs_h2_err = 0.2*self.obs_bary_dwarfs_h2
        self.obs_bary_dwarfs_beta = 0.244
        self.obs_bary_dwarfs_beta_err = 0.2*self.obs_bary_dwarfs_beta


        # Simplenu DM disc model priors
        self.simplenu_dm_K_max = 1500.  #JR model has K = 300.
        self.simplenu_dm_K_min = 0.
        self.simplenu_dm_D_max = 3.5  #JR model has D = 2.5
        self.simplenu_dm_D_min = 1.5

        # Tilt priors
        self.tilt_A_max = 270.  # 'Right answer': A = 181.77 
        self.tilt_A_min = 90
        self.tilt_n_max = 1.9    # simple2 mock has n = 1.44
        self.tilt_n_min = 1.
        self.tilt_R_max = 3.5     # simple2 mock has R = 2.5
        self.tilt_R_min = 1.5

        # MultiNest options
        # ----------------------------------------------------------------------
        self.map_priors = False

        # live points, > ndim, < 2^ndim, about number of
        # ellipsoids in phase space to be found
        self.nlive = 100*self.ndim
        self.err = 1e300    # chi^2 for models which are impossible

        #fraction of profiles to save, set <0 for no profile saving
        self.save_fraction = -1.0

        #Plotting flag
        self.plotting_flag = False

        # filesystem-related
        # ----------------------------------------------------------------------
        import import_path as ip
        ip.set_geometry(self.geom, self.machine) # load spherical or
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

        #self.z_bincenters = np.array([]) # [pc] holds the bin centers, H Silverwood 21/11/14
        #self.z_binmins = np.array([])
        #self.z_binmaxs = np.array([])
        #self.z_all_pts = np.array([]) # [pc] holds [zC = 0, z_bin_centers]


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
