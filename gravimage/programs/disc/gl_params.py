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
    def __init__(self, timestamp = '', investigate = '', case = -1, mock_suffix = 0):

        myrank = MPI.COMM_WORLD.Get_rank()
        nprocs = MPI.COMM_WORLD.Get_size()
        procnm = MPI.Get_processor_name()

        #print('P', myrank, ': Initiating Params()')

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
        self.mock_suffix = mock_suffix

        # debug options
        # ----------------------------------------------------------------------
        self.checksig = False # debug sig calculation?
        self.debug = False # stop at wrong sanitazion?


        # data and analysis options
        # ----------------------------------------------------------------------
        self.getnewdata = True  # get new data computed from
                                # observations before burn-in
        self.getnewpos  = True  # redo the first data conversion step

        no_pops = 1   # Number of pops, can currently only take the value 1 or 2
        
        switch_pop_order = False   # Switch which pop is which by toggling this flag
        # If False: young comes first (i.e. young data used if no_pops = 1)
     

        self.raw_data = False
        self.fit_to_sigz2 = False 
        # set if chi2 is calculated on sigz2 or sigz, correlate with data folder
        # will then save sigz vectors (but will still be referred to as sigz2 !)
        # Can only be false if raw_data = False !!

        if self.raw_data == True:
            print ('Binning raw data read in from file')
            #self.data_z_cut = [2.4] * no_pops  # [kpz] only use (& bin) data up to this z
            self.data_z_cut = [1.3] * no_pops
            if self.investigate == 'simplenu':
                ext_data_file_1, ext_data_file_tilt_1 = \
                    gh.ext_file_selector_simplenu([2], '1e4', '', True, True, 'simplenu_baryon', 0)   # First population  (running without tilt)
            #    gh.ext_file_selector_simplenu([2], '1e5', '', False, True, 'simplenu_baryon', mock_suffix)   # First population  (running without tilt)
                if no_pops >= 2:
                    ext_data_file_2, ext_data_file_tilt_2 = \
                        gh.ext_file_selector_simplenu([2], '1e4', '', True, True, 'simplenu_baryon', 1)   # Second population
                if no_pops >= 3:
                    ext_data_file_3, ext_data_file_tilt_3 = \
                        gh.ext_file_selector_simplenu([2], '1e4', '', True, True, 'simplenu_baryon', 2)   # Third population
                if no_pops ==1:
                    self.external_data_file = ext_data_file_1  # (1 pop)
                    self.external_data_file_tilt = ext_data_file_tilt_1  # (1 pop)
                elif no_pops == 2:
                    self.external_data_file = ext_data_file_1 + ext_data_file_2 
                    self.external_data_file_tilt = ext_data_file_tilt_1 + ext_data_file_tilt_2
                elif no_pops == 3:
                    self.external_data_file = ext_data_file_1 + ext_data_file_2 + ext_data_file_3 
                    self.external_data_file_tilt = ext_data_file_tilt_1 + ext_data_file_tilt_2 + ext_data_file_tilt_3

           # self.external_data_file, self.external_data_file_tilt = \
           #     gh.ext_file_selector_simplenu([2], '1e6', '', True, True, 'simplenu_baryon', mock_suffix) #  # Running the code with tilt

                #ext_file_selector_simplenu(pops, sampling, darkdisk, tilt, pos_sigRz, baryon_model)
                #pops = [x,x], eg which population to use  ?
                #sampling = 1e4, 1e5, 1e6,
                #darkdisk= '', 'dd', 'bdd'
                #tilt = True or False
                #pos_sigRz = True or False, to pull in data set with old or new (positive) sigRz
                #baryon_model = simplenu_baryon or obs_baryon
                #suffix


            elif self.investigate == 'disc_nbody' and self.case == 1: ######## Hunt & Kawata Data sets ########
                self.external_data_file = ['/Hunt_Kawata_GCD/HK_cut_1kpc_GalC_z_vz.dat']
                self.external_data_file_tilt = ['/Hunt_Kawata_GCD/HK_cut_1kpc_GalC_z_vz.dat']

            elif self.investigate == 'disc_nbody' and self.case == 2: ######## Garbari et al. Wedges ########
                Gar_Root = '/Garbari_Nbody/wedges/'
                Gar_Simulation = 'mwhr' #'unevolved', 'mwhr', 'LLMC10', 'LLMC60'
                Gar_Radius = 'r8500' #'r8500', 'r10500', 'r12500'
                Gar_Angle = 'ang0' #'ang45', 'ang90', 'ang135', 'ang180', 'ang225', 'ang270', 'ang315'
                self.external_data_file = [Gar_Root + Gar_Simulation + '_' + Gar_Radius + '_' + Gar_Angle + '_stars_z_vz_vR.dat']

        else:  # if raw_data == False:
            # Assuming 1 or 2 populations:
            # Mock data extracted binned data, preparing for Budenbender data
            print ('Using values extracted from already binned data')

            self.ext_z_vec = [None] * no_pops
            self.ext_nu_z_vec = [None] * no_pops # z where nu data is defined
            self.ext_nu_vec = [None] * no_pops  ; self.ext_nu_err_vec = [None] * no_pops
            self.ext_sigz2_vec = [None] * no_pops  ; self.ext_sigz2_err_vec = [None] * no_pops
            self.ext_sigRz2_vec = [None] * no_pops  ; self.ext_sigRz2_err_vec = [None] * no_pops
            binned_data_path = '/home/sofia/darcoda/Data_Sets//binned_data/'
            if self.fit_to_sigz2 == True: 
                binned_data_set = 'simplenu_mock_sigz2/'
            else:     
                binned_data_set = 'buden_meas/' # _0.dat: alpha-young, and _1.dat: alpha-old


            path_0 = binned_data_path + binned_data_set + 'binned_data_0_short.dat'
            nu_path_0 = binned_data_path + binned_data_set + 'binned_nu_data_0.dat' #or short
            bins_0 = [20,9]
            path_1 = binned_data_path + binned_data_set + 'binned_data_1.dat'
            nu_path_1 = binned_data_path + binned_data_set + 'binned_nu_data_1.dat'
            bins_1 = [20,10]
            #no_bins_0 =  [20,10]  # If above is: binned_nu_data_0.dat & binned_data_0.dat
            #no_bins_0 =  [9,9]    # If using: binned_data_0_short.dat & binned_nu_data_0_short.dat above
            if switch_pop_order:
                ext_binned_file_1 = path_0 ; ext_binned_nu_file_1 = nu_path_0 ; no_bins_1 = bins_0
                ext_binned_file_0 = path_1 ; ext_binned_nu_file_0 = nu_path_1 ; no_bins_0 = bins_1
            else:
                ext_binned_file_0 = path_0 ; ext_binned_nu_file_0 = nu_path_0 ; no_bins_0 = bins_0
                ext_binned_file_1 = path_1 ; ext_binned_nu_file_1 = nu_path_1 ; no_bins_1 = bins_1


            binned_data_0 = np.loadtxt(ext_binned_file_0, unpack=True, skiprows=1) 
            binned_data_1 = np.loadtxt(ext_binned_file_1, unpack=True, skiprows=1) 
            binned_nu_data_0 = np.loadtxt(ext_binned_nu_file_0, unpack=True, skiprows=1) 
            binned_nu_data_1 = np.loadtxt(ext_binned_nu_file_1, unpack=True, skiprows=1) 
            print ('binned_data_0:',binned_data_0)
            #print ('b_d_0[0]:',binned_data_0[0])
            print ('binned_data_1:',binned_data_1)
            #print ('binned_nu_data_0:',binned_data_0)
            #print ('binned_nu_data_1:',binned_data_1)
            #pdb.set_trace()
            #Continue job reading binned data from file

            self.ext_z_vec[0] = binned_data_0[0]
            self.ext_nu_z_vec[0] = binned_nu_data_0[0]
            self.ext_nu_vec[0] = binned_nu_data_0[1]
            self.ext_nu_err_vec[0] = binned_nu_data_0[2]

            self.ext_sigz2_vec[0] = binned_data_0[1]
            self.ext_sigz2_err_vec[0] = binned_data_0[2]
            self.ext_sigRz2_vec[0] = binned_data_0[3]
            self.ext_sigRz2_err_vec[0] = binned_data_0[4]

            if no_pops == 2:
                self.ext_z_vec[1] = binned_data_1[0]
                self.ext_nu_z_vec[1] = binned_nu_data_1[0]
                self.ext_nu_vec[1] = binned_nu_data_1[1]
                self.ext_nu_err_vec[1] = binned_nu_data_1[2]

                self.ext_sigz2_vec[1] = binned_data_1[1]
                self.ext_sigz2_err_vec[1] = binned_data_1[2]
                self.ext_sigRz2_vec[1] = binned_data_1[3]
                self.ext_sigRz2_err_vec[1] = binned_data_1[4]

        #Check case agains data file
        if self.raw_data:
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
        if self.raw_data:
            self.ntracer_pops = len(self.external_data_file) # number of stellar tracer populations
            print ('self.external_data_file:',self.external_data_file)
            if len(self.data_z_cut) != self.ntracer_pops:
                raise Exception('Incorrect data_z_cut vector length')
        else:
            self.ntracer_pops = no_pops
            
        print ('self.ntracer_pops:',self.ntracer_pops)
        #pdb.set_trace()

        #self.tilt = False   # If also modelling the tilt
        self.tilt = True   # If also modelling the tilt

        self.darkmattermodel = 'const_dm' # const_dm = const DM dens in z
        #self.darkmattermodel = 'ConstPlusDD' # constant DM + DM disc component
        #self.darkmattermodel = 'gaussian_per_bin' #Gaussian in each bin w/ median & SD equal in each bin
        #self.darkmattermodel = 'kz_dm'  # kz_dm = kz parameterization of DM CURRENTLY NOT WORKING

        #self.dd_data = False # if we are to plot dd analytics or not

        self.binning = 'consttr' # 'linspace', 'logspace', 'consttr': binning of particles
        #self.binning = 'linspace' # 'linspace', 'logspace', 'consttr': binning of particl
        self.no_Sigrho_bins = 20  # number of bins to store rho & Sig profiles in
        if self.raw_data:  
            no_bins = 18
            self.nbins = [[no_bins,no_bins]] * no_pops #(no bins MUST be same for nu & sigz)
        else:
            #self.nbins = [[20,10]] * no_pops # first nu-bins then sigz-bins
            self.nbins = []
            self.nbins.append(no_bins_0)    
            if no_pops == 2:
                self.nbins.append(no_bins_1)
        # we allow differnt no bins for nu and sigz (assume sigRz to have same bins as sigz)

        if len(self.nbins) != self.ntracer_pops:
            raise Exception('Number of bins vector is of incorrect length')
        #self.nrhonu = sum(self.nbins) + 1 # Number of points where rho and nu parameters
        #will be set, e.g. bin centres, plus zC=0  # SS: no longer used (outside gl_params)

        #Tracer Density description
        #self.nu_model = 'kz_nu'
        #self.nu_model = 'gaussian_data'
        self.nu_model = 'exponential_sum'
        self.N_nu_model_exps = 1  # Can only take values 1 or 2. If 2: flattening nu for low z

        #sigRz sign correction
        self.positive_sigRz_data_sign = True
        self.positive_sigRz_model_sign = True

        #------------------------
        self.TheoryData = False   # If true using theoretical bin values as indata
        self.hyperparams = False  # Use hyperparameters, 2 params, range (0.1->10)*meanerr
        #Dark matter options
        self.adddarkdisc = False  # for disc mock case: add a dark disc?
        # self.adddarkdisc is currently not used !
        self.param_headers = []



        #Baryon options
        self.baryonmodel = 'obs_baryon' # Hopefully realistic baryon model, based on observations
        self.analytic_sigz2 = True   # Calc sigz2 without the use of the C integration constant
        # If True this requires: baryonmodel = obs_baryon and N_nu_model_exps = 1
        # and darkmattermodel = const_dm

        #self.baryonmodel = 'trivial_baryon' # Only baryons inside innermost bin
        #self.baryonmodel = 'simplenu_baryon' # model for simplenu mock data
                                    # simplenu_baryon_gaussian Gaussian rho
        self.nbaryon_pops = 1 # Number of baryon populations to look at
                                    # =0 if doing simple mass model (eg DM profile describes
                                    # all mmass)
        # ASSUMING nbaryons_pops = 1   major rewrite to change this

        # Number of parameters to describe baryon population
        if self.baryonmodel == 'simplenu_baryon':
            self.nbaryon_params = 2
        elif self.baryonmodel == 'trivial_baryon':
            self.nbaryon_params = 1
        elif self.baryonmodel == 'obs_baryon':
            self.nbaryon_params = 10
        elif self.baryonmodel == 'simplenu_baryon_gaussian': # Not supported
            self.nbaryon_params = sum(self.nbins)+1

        if self.darkmattermodel == 'const_dm':
            self.nDM_params = 1
        elif self.darkmattermodel == 'ConstPlusDD':
            self.nDM_params = 3

        self.ntilt_params = 3  # Number of parameters used to describe tilt

        self.scan_rhonu_space = False #Search directly in rho or nu space, i.e. no kr parametrization
        if self.scan_rhonu_space: # Currently not supported
            self.nrhonu = self.nbins #Param count will be rho_C + rho_bins

        # Parameters escribing the underlying baryonic model and the mock data
        # This is assuming baryonmodel = simplenu_baryon above, rewrite required if this is changed.

        #Total dimensions count
        if self.darkmattermodel == 'kz_dm': # not supported
            self.ndim = self.ntracer_pops + (self.nrhonu + 1) + (sum(self.nbins) + self.ntracer_pops) + self.nbaryon_pops*self.nbaryon_params
                # Constant C from sigma_z calculations (one for each tracer population)
                # DM kz_rho defined at all bin centres and zC=0, plus rho_C
                # kz_nu definitions at bin centres, plus at zC = 0 for each population
                # baryon parameters

        elif self.darkmattermodel in ['const_dm', 'ConstPlusDD']:
            #self.ndim = self.ntracer_pops + 1 + (sum(self.nbins) + 2*self.ntracer_pops) + self.nbaryon_pops*self.nbaryon_params
            self.ndim = self.nDM_params + self.nbaryon_pops*self.nbaryon_params
                # DM params
                # baryon parameters

        elif self.darkmattermodel == 'gaussian_per_bin': # not supported
            self.ndim = self.ntracer_pops + (self.nrhonu) + (sum(self.nbins) + 2*self.ntracer_pops) + self.nbaryon_pops*self.nbaryon_params
                # Constant C from sigma_z calculations (one for each tracer population)
                # DM density in bin centres plus zC = 0
                # kz_nu definitions at bin centres, plus at zC = 0 for each population
                # baryon parameters

        if self.analytic_sigz2 == False:
            self.ndim += self.ntracer_pops
            # Constant C from sigma_z calculations (one for each tracer population)

        if self.nu_model == 'gaussian_data':  # Not supported
            self.ndim -= self.ntracer_pops
        elif self.nu_model == 'exponential_sum':
            #self.ndim = self.ndim - (sum(self.nbins) + 2*self.ntracer_pops) + 2*self.N_nu_model_exps*self.ntracer_pops   # SS added last ntracer_pops 25 apr -16  (one nu-model per pop)
            self.ndim = self.ndim + 2*self.N_nu_model_exps*self.ntracer_pops   # SS added last ntracer_pops 25 apr -16  (one nu-model per pop)
        # else: not supproted, would give wrong ndim when using only the above

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
        # DM density prior (for constant DM dens):
        self.rho_C_max = 40.0E6  # WAS 20 !  #Msun kpc^-3 
        self.rho_C_min = 2.0E6 #Msun pc^-3 #standard is 1e6 to 1E8
        self.rho_C_prior_type = 'linear' #log, linear, gaussian

        self.Sigma_DD_max = 40.E6  # Not (yet?) used
        self.Sigma_DD_min = 0.
        self.h_dd_min = 0.5
        self.h_DD_max = 2.

        # Below 3 lines are currently not used:
        self.nu_C_max = 0.0 # NOT USED # no. stars pc^-3, set in external_data
        self.nu_C_min = 0.0 # NOT USED # no. stars pc^-3
        self.nu_C_prior_type = 'log' # NOT USED

        #Tracer density exponential sum prior
        #  For each exponential: [nuC_median, nuC_SD, nu_h, nu_h_SD]
        self.nu_exp_sum_priors=[]
        for pop in range(self.ntracer_pops): # SS: One prior array per pop
            nu_exp_sum_priors_pop = []
            #for ii in range(0, self.N_nu_model_exps): nu_exp_sum_priors_pop.append([0, 0, 0.9, 0.5]) #nuC and nuC_SD set in external data
            # Below: no longer set in external_data:
            for ii in range(0, self.N_nu_model_exps): nu_exp_sum_priors_pop.append([-1., 2., 0.6, 0.5]) #log10(nuC) median, log10(nuC) range, h median, h range
            self.nu_exp_sum_priors.append(nu_exp_sum_priors_pop)
        #print ('nu_exp_sum_priors:',self.nu_exp_sum_priors)
        #if not len(self.nu_exp_sum_priors) == self.N_nu_model_exps:
        #    raise Exception('Wrong number of priors for nu exponentials')
        # SS: No longer true, but check not needed since given explicitly above


        #Priors for Gaussian DM
        self.rho_DM_gaussian_med = 10.0E6
        self.rho_DM_gaussian_SD = 5.0E6

        # Limits for central kz values (z=0)
        self.kz_rho_C_max = 5.0
        self.kz_rho_C_min = -5.0 #SS
        self.kz_nu_C_max = 5.#20.0
        self.kz_nu_C_min = -5.0 #SS

        # Maximum kz_slope (=dk/dz)
        #self.max_kz_slope = 10.0
        self.max_kz_slope = 5.0

        # Limits for sigz central value
        self.sigz_C_max = [None]*self.ntracer_pops # set from data in gr_external_data
        self.sigz_C_min = [None]*self.ntracer_pops # set separately for each pop

        # Monotonicity priors
        self.monotonic_rho = False    # mono-prior on rho(z)
        self.monotonic_nu = False # mono-prior on nu(z)

        # kz selection scheme
        self.kz_rho_selection = 'gaussian'
        self.kz_nu_selection = 'tophat'

        # Log or linear priors for rhonu scanning
        self.prior_type_rho = 'gaussian' # NOT USED  # 'log' or 'linear'
        self.prior_type_nu = 'gaussian' # NOT UDSED  # 'log' or 'linear'

        # Simplenu Baryon model priors
        #self.prior_type_simplenu_baryon = 'gaussian' # gaussian or 'linear'
        self.prior_type_simplenu_baryon = 'linear' # gaussian or 'linear'
        # See if this switch solves infinity-problems

        self.gaussian_rho_baryon_mid_vector = [] #def values in gl_data.read_nu
        self.gaussian_rho_baryon_SD_vector = [] #define in gl_data.read_nu
        self.gaussian_baryon_prior_SD = 0.1 #

        #       Linear limits  # widened the range
        self.simplenu_baryon_K_max = 2000 #1650  #JR model has K = 1500.
        self.simplenu_baryon_K_min = 1000 #1350
        self.simplenu_baryon_D_max = 0.3  #0.2  #JR model has D = 0.18
        self.simplenu_baryon_D_min = 0.1 #0.16

        #self.simplenu_baryon_K_max = 3000  #1600 #1700. #JR model has K = 1500.
        #self.simplenu_baryon_K_min = 100  #1400 #1300.
        #self.simplenu_baryon_D_max = 0.30  #0.20 #0.5 #JR model has D = 0.18
        #self.simplenu_baryon_D_min = 0.12  #0.16 #0.05

        #       Gaussian SD
        self.simplenu_baryon_K_mid = 1500
        self.simplenu_baryon_K_sd = 150
        self.simplenu_baryon_D_mid = 0.18
        self.simplenu_baryon_D_sd = 0.02

        #self.simplenu_baryon_K_sd = 75
        #self.simplenu_baryon_D_sd = 0.01
        #self.simplenu_baryon_K_sd = 0.15
        #self.simplenu_baryon_D_sd = 0.00002

        # trivial baryonic profile: all baryons inside innermost bin
        #self.trivial_baryon_Sig_min = 43.5E6 # Tight prior, McKee level
        #self.trivial_baryon_Sig_max = 50.3E6 
        self.trivial_baryon_Sig_min = 15.E6  # Trying wide bary prior
        self.trivial_baryon_Sig_max = 80.E6 

        # Observational baryon model priors: -----------------------------
        self.obs_baryon_Sig_tot       = 46.85E6  # [Msun/kpc2] # 10 parameters for obs_baryon
        self.wide_bary_range = True
        if self.wide_bary_range:
            self.obs_baryon_Sig_tot_err   = 0.7
        else:
            self.obs_baryon_Sig_tot_err   = 0.13  # Fraction 

        self.obs_baryon_Sig_dwarf     = 23.7E6  # McKee, F06 value: 23.0
        self.obs_baryon_Sig_dwarf_err = 0.2
        self.obs_baryon_rho0_MS       = (0.0029+0.0072+0.0006)*1e9  # [Msun/kpc3]
        self.obs_baryon_MS_beta_max   = 0.3  # beta in range 0 to beta_max 

        self.obs_baryon_h      = 400./1000. # dwarf h [kpc]
        self.obs_baryon_h_err  = 0.05 # Error in percentage of central value
        self.obs_baryon_h1     = 337./1000. #[kpc] # dwarf thin disk
        self.obs_baryon_h1_err = 0.1
        self.obs_baryon_h2     = 609./1000. 
        self.obs_baryon_h2_err = 0.2
        self.obs_baryon_h3     = 1000./1000.
        self.obs_baryon_h3_err = 0.2
        self.obs_baryon_x      = 0.5  # Fraction of thick disk modelled with h1 or h2
        self.obs_baryon_x_err  = 0.5 #   can take all values btwn 0 and 1

        self.obs_baryon_Sig_HII     = 1.7E6  # 1.8 if including Gum nebula. Include also halo stars?? 
        self.obs_baryon_Sig_HII_err = 0.2
        self.obs_baryon_h_HII       = 1590./1000.
        self.obs_baryon_h_HII_err   = 0.2
        # -------------------------------------------------------------------


        # Simplenu DM disc model priors
        #self.simplenu_dm_K_max = 1500.  #JR model has K = 300.
        #self.simplenu_dm_K_min = 0.
        self.simplenu_DD_Sig_inf_max = 40.E6 # Very big dark disk, adds to const DM so can give
        self.simplenu_DD_Sig_inf_min = 0.    #  very large DM densities
        self.simplenu_dm_D_max = 3.5  #JR model has D = 2.5
        self.simplenu_dm_D_min = 0.5  # Was 1.5

        # Tilt priors
        #self.tilt_A_max = -0.005  # simple2 mock has A = -0.0087
        #self.tilt_A_min = -0.012
        #self.tilt_n_max = 1.9    # simple2 mock has n = 1.44
        #self.tilt_n_min = 1.
        #self.tilt_R_max = 3.5     # simple2 mock has R = 2.5
        #self.tilt_R_min = 1.5

        if self.positive_sigRz_model_sign:  # True
            #self.tilt_A_max = 200.00 #A 180.08335058 
            self.tilt_A_max = 400.00 # Was 200 above before 
            self.tilt_A_min = 0.0
        else:
            self.tilt_A_max = 0.0  # old simple2 mock has A = -181.77 for z[kpc]
            self.tilt_A_min = -200.00

        #self.tilt_n_max = 3.0    # simple2 mock has n = 1.43829774
        self.tilt_n_max = 1.5   
        #self.tilt_n_min = 0.1
        self.tilt_n_min = 0.5 

        #self.prior_type_tilt_R = 'gaussian' # do 'linear' instead???
        self.prior_type_tilt_R = 'linear' # Was 'gaussian' above before !!!!
        self.tilt_R_max = 4.5    # loglike uses (1/Rsun - 2/R)
        self.tilt_R_min = 0.5    # if linear tilt_prior ( simple2 mock has R = 2.5) 


        #self.tilt_k_min = np.array([-1.3,-0.5])  # young: k: -1.3 -> 1  , old: k: -0.5 -> 1.5
        #self.tilt_k_max = np.array([1.,1.5])     # linear prior.
        # k = k_1 + k_2, functions modelled as exp(-k*R). Loglike uses (1/Rsun - k)
        # young: -0.3 < k_0 < 0 , old: k_0 = 0.5
        # both pop: -1 < k_1 < 1:
        young_k_min = -1.3  ; young_k_max = 1.
        old_k_min   = -0.5  ; old_k_max = 1.5
        # If instead -2 < k_1 < 2:
        #young_k_min = -2.3 ; young_k_max = 2.
        #old_k_min  =  -1.5 ; old_k_max   = 2.5
        
        if no_pops == 1:
            if switch_pop_order:
                self.tilt_k_min = np.array([old_k_min])
                self.tilt_k_max = np.array([old_k_max])
            else:
                self.tilt_k_min = np.array([young_k_min])
                self.tilt_k_max = np.array([young_k_max])
        else:   # no_pops = 2, code currently does not allow no_pops larger than 2
            self.tilt_k_min = np.array([young_k_min,old_k_min])
            self.tilt_k_max = np.array([young_k_max,old_k_max])
            if switch_pop_order:
                self.tilt_k_min = self.tilt_k_min[::-1]
                self.tilt_k_max = self.tilt_k_max[::-1]                    

        self.tilt_R_med = 2.5  # If gaussian tilt_prior
        self.tilt_R_sd = 0.5


        # MultiNest options
        # ----------------------------------------------------------------------
        self.map_priors = False

        # live points, > ndim, < 2^ndim, about number of
        # ellipsoids in phase space to be found
        self.nlive = 100*self.ndim

        if self.map_priors:
            self.nlive = int(1E4)

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

        gh.LOG(1,'Instantiating Files class')

        import gl_class_files
        self.files = gl_class_files.Files(self, timestamp)
        from gl_data import Datafile
        self.dat = Datafile()

        # global arrays
        # ----------------------------------------------------------------------
        self.z_bincenter_vecs = [] # [(np.zeros(self.nbins[ii])) for ii in range(0, self.ntracer_pops)]
        self.nu_z_bincenter_vecs = []  # z vectors for the nu data
        self.z_vec_Sigrho = []  # z values where Sigma and rho are given
        self.extend_z_binc_vecs = [] # z_bincenter_vecs with added zeros in the front
        self.extend_nu_z_binc_vecs = [] # z_bincenter_vecs with added zeros in the front
        self.z_binmin_vecs = [] #[(np.zeros(self.nbins[ii])) for ii in range(0, self.ntracer_pops)]
        self.z_binmax_vecs = [] #[(np.zeros(self.nbins[ii])) for ii in range(0, self.ntracer_pops)]
        self.z_all_pts_unsort = np.array([]) #Simple concatenation of bincentres and zC = 0
        self.z_all_pts_sorted = [] #Sorted list of bincentres and zC = 0
        self.z_vec_masks = [] #[[].append(None) for ii in range(0, self.ntracer_pops)]

        # Global constants
        self.Rsun = 8.  # [kpc]  Sun's distance to galactic center

        self.print_flag = False  # If the code should print things for testing or not

    ## \fn __init__(self, timestamp = '')
    # set up all parameters used in the course of the MultiNest run,
    # and the analysis routines in plot_multinest @param timestamp =
    # '', for output


    def __repr__(self):
        return "Params: "+self.files.dir
    ## \fn __repr__(self)
    # string representation for ipython
