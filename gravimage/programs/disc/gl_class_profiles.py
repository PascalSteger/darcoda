#!/usr/bin/env ipython3

##
# @file
# store profiles

# (c) 2013 ETHZ psteger@phys.ethz.ch

import pdb
import numpy as np
import gl_project
import gl_physics as phys
import gl_helper as gh

class Profiles:
    def __init__(self, ntracer_pops, nbins):#, nrhonu,  nbaryon_pops, nbaryon_params):
        self.ntracer_pops = ntracer_pops
        self.nbins = nbins

        #z-profile points
        self.z_C = 0.0
        self.z_vecs = [(np.zeros(nbins[ii])) for ii in range(0, ntracer_pops)]
        self.z_vecs_comb_w0 = np.zeros(sum(nbins)+1)
        self.binmins = [(np.zeros(nbins[ii])) for ii in range(0, ntracer_pops)] #?
        self.binmaxs = [(np.zeros(nbins[ii])) for ii in range(0, ntracer_pops)] #?
        self.z_vec_masks = [[].append(None) for ii in range(0, ntracer_pops)]

        #Dark Matter profile parameters and derived mass density
        self.rho_DM_C      = 0.0               #Multinest
        self.kz_rho_DM_C   = 0.0               #Multinest
        self.kz_rho_DM_vec = np.zeros(sum(nbins))   #Multinest
        self.rho_DM_vec    = np.zeros(sum(nbins))   #Derived from phys
        self.Sig_DM_C      = 0.0               #always zero if z_C = 0
        self.Sig_DM_vec    = np.zeros(sum(nbins))   #Derived from phys

        #Baryon mass density and parameters
        #self.baryon_params  = np.zeros(nbaryon_pops*nbaryon_params)
        self.rho_baryon_C   = 0.0
        self.rho_baryon_vec = np.zeros(sum(nbins))
        self.Sig_baryon_C = 0.0                 #always zero if z_C = 0
        self.Sig_baryon_vec = np.zeros(sum(nbins))

        #Total Mass density
        self.rho_total_C = 0.0
        self.rho_total_vec = np.zeros(sum(nbins))
        self.Sig_total_C = 0.0                 #always zero if z_C = 0
        self.Sig_total_vec = np.zeros(sum(nbins))

        #Tracer profile parameters and derived mass density
        self.norm_C    = np.zeros(ntracer_pops)
        self.nu_C      = np.zeros(ntracer_pops)            #Multinest
        self.kz_nu_C   = np.zeros(ntracer_pops)            #Multinest

        self.kz_nu_vecs = [(np.zeros(nbins[ii])) for ii in range(0, ntracer_pops)] #Multinest
        self.nu_vecs    = [(np.zeros(nbins[ii])) for ii in range(0, ntracer_pops)] #Derived from phys
        self.sigz2_vecs = [(np.zeros(nbins[ii])) for ii in range(0, ntracer_pops)] #Derived from phys

        #chi2 of profile
        self.chi2 = 0.0
        self.chi2_nu_vecs     = [(np.zeros(nbins[ii])) for ii in range(0, ntracer_pops)]
        self.chi2_sigz2_vecs  = [(np.zeros(nbins[ii])) for ii in range(0, ntracer_pops)]
        self.chi2_sigRz2_vecs = [(np.zeros(nbins[ii])) for ii in range(0, ntracer_pops)]

        #weight of profile (used in multinest analysis)
        self.mn_weight = 0.0

        #Hyperparameters
        self.hyper_nu = 0.
        self.hyper_sigz2 = 0.

        #Tilt
        self.tilt_vecs     = [(np.zeros(nbins[ii])) for ii in range(0, ntracer_pops)]
        self.sigmaRz2_vecs  = [(np.zeros(nbins[ii])) for ii in range(0, ntracer_pops)]

    ## \fn __init__(self, pops, nipol)
    # constructor
    # @param ntracer_pops = number of tracer populations
    # @param nbins = number of vertical bins

    # set_prof just for the vectors, so that the length can be checked

    def set_prof(self, prof, vec, pop, gp):
        if prof in ['kz_rho_DM_vec', 'rho_DM_vec', 'Sig_DM_vec', 'rho_baryon_vec', 'Sig_baryon_vec',
            'rho_total_vec', 'Sig_total_vec']:
            proper_vec_length = gp.nrhonu - 1
        elif prof in ['z_vecs_comb_w0']:
            proper_vec_length = gp.nrhonu
        elif prof in ['z_vecs', 'kz_nu_vecs', 'nu_vecs', 'sigz2_vecs', 'tilt_vecs', 'sigmaRz2_vecs',\
                        'chi2_nu_vecs', 'chi2_sigz2_vecs', 'chi2_sigRz2_vecs']:
            proper_vec_length = gp.nbins[pop]

        gh.sanitize_vector(vec, proper_vec_length, -1e30, 1e30, gp.debug)

        #z-profile points
        if prof == 'z_vecs':
            self.z_vecs[pop] = vec
        elif prof == 'z_vecs_comb_w0':
            self.z_vecs_comb_w0 = vec

        #Dark matter
        elif prof == 'kz_rho_DM_vec':
            self.kz_rho_DM_vec = vec
        elif prof == 'rho_DM_vec':
            self.rho_DM_vec = vec
        elif prof == 'Sig_DM_vec':
            self.Sig_DM_vec = vec

        #Baryons
        elif prof == 'baryon_params':
            print("Still working on baryons, might be a tad dodgy")
            self.baryon_params = vec
        elif prof == 'rho_baryon_vec':
            self.rho_baryon_vec = vec
        elif prof == 'Sig_baryon_vec':
            self.Sig_baryon_vec = vec

        #Total mass
        elif prof == 'rho_total_vec':
            self.rho_total_vec = vec
        elif prof == 'Sig_total_vec':
            self.Sig_total_vec = vec

        #Tracer stars
        elif prof == 'kz_nu_vecs':
            self.kz_nu_vecs[pop] = vec
        elif prof == 'nu_vecs':
            self.nu_vecs[pop]  = vec
        elif prof == 'sigz2_vecs':
            self.sigz2_vecs[pop] = vec


        #Tilt
        elif prof == 'tilt_vecs':
            self.tilt_vecs[pop] = vec
        elif prof == 'sigmaRz2_vecs':
            self.sigmaRz2_vecs[pop] = vec

        #Chi2 profiles
        elif prof == 'chi2_nu_vecs':
            self.chi2_nu_vecs[pop] = vec
        elif prof == 'chi2_sigz2_vecs':
            self.chi2_sigz2_vecs[pop] = vec
        elif prof == 'chi2_sigRz2_vecs':
            self.chi2_sigRz2_vecs[pop] = vec



    ## \fn set_prof(self, prof, vec, pop, gp)
    # store density vector
    # @param prof profile identifier
    # @param vec array of floats
    # @param pop population, if applicable
    # @param gp global parameters, used for gp.xepol radii


    def get_prof(self, prof, pop):
        #z-profile points
        if prof == 'z_vec':
            return self.z_vec

        #Dark matter
        elif prof == 'kz_rho_DM_vec':
            return self.kz_rho_DM_vec
        elif prof == 'rho_DM_vec':
            return self.rho_DM_vec
        elif prof == 'Sig_DM_vec':
            return self.Sig_DM_vec

        #Baryons
        elif prof == 'baryon_params':
            print("Still working on baryons, might be a tad dodgy")
            return self.baryon_params
        elif prof == 'rho_baryon_vec':
            return self.rho_baryon_vec
        elif prof == 'Sig_baryon_vec':
            return self.Sig_baryon_vec

        #Total mass
        elif prof == 'rho_total_vec':
            return self.rho_total_vec
        elif prof == 'Sig_total_vec':
            return self.Sig_total_vec

        #Tracer stars
        elif prof == 'kz_nu_vecs':
            return self.kz_nu_vecs[pop]
        elif prof == 'nu_vecs':
            return self.nu_vecs[pop]
        elif prof == 'sigz2_vecs':
            return self.sigz2_vecs[pop]

        #Tilt
        elif prof == 'tilt_vecs':
            return self.tilt_vecs[pop]
        elif prof == 'sigmaRz2_vecs':
            return self.sigmaRz2_vecs[pop]

        #Chi2 profile
        elif prof == 'chi2_nu_vecs':
            return self.chi2_nu_vecs[pop]
        elif prof == 'chi2_sigz2_vecs':
            return self.chi2_sigz2_vecs[pop]
        elif prof == 'chi2_sigRz2_vecs':
            return self.chi2_sigRz2_vecs[pop]

    ## \fn get_prof(self, prof, pop)
    # return density array
    # @param prof of this profile
    # @param pop and of this population (if applicable)

    def __repr__(self):
        return "Profiles (disc): "+str(self.chi2)
    ## \fn __repr__(self)
    # string representation for ipython


## \class Profiles
# class for storing all profiles
