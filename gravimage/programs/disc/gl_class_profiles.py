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
        self.z_vec = np.zeros(nbins)
        self.binmin = np.zeros(nbins) #?
        self.binmax = np.zeros(nbins) #?
        self.xbins = np.zeros(nbins+1) #to mesh with general code
        self.x0 = np.zeros(nbins) #to mesh with general code

        #Dark Matter profile parameters and derived mass density
        self.rho_DM_C      = 0.0               #Multinest
        self.kz_rho_DM_C   = 0.0               #Multinest
        self.kz_rho_DM_vec = np.zeros(nbins)   #Multinest
        self.rho_DM_vec    = np.zeros(nbins)   #Derived from phys
        self.Sig_DM_C      = 0.0               #always zero if z_C = 0
        self.Sig_DM_vec    = np.zeros(nbins)   #Derived from phys

        #Baryon mass density and parameters
        #self.baryon_params  = np.zeros(nbaryon_pops*nbaryon_params)
        self.rho_baryon_C   = 0.0
        self.rho_baryon_vec = np.zeros(nbins)
        self.Sig_baryon_C = 0.0                 #always zero if z_C = 0
        self.Sig_baryon_vec = np.zeros(nbins)

        #Total Mass density
        self.rho_total_C = 0.0
        self.rho_total_vec = np.zeros(nbins)
        self.Sig_total_C = 0.0                 #always zero if z_C = 0
        self.Sig_total_vec = np.zeros(nbins)

        #Tracer profile parameters and derived mass density
        self.nu_C      = 0.0                               #Multinest
        self.kz_nu_C   = 0.0                               #Multinest
        self.kz_nu_vec = np.zeros(ntracer_pops * nbins)    #Multinest
        self.nu_vec    = np.zeros(ntracer_pops * nbins)    #Derived from phys
        self.sigz2_vec = np.zeros(ntracer_pops * nbins)    #Derived from phys

        #chi2 of profile
        self.chi2 = 0.0

        #Hyperparameters
        self.hyper_nu = 0.
        self.hyper_sigz2 = 0.
    ## \fn __init__(self, pops, nipol)
    # constructor
    # @param ntracer_pops = number of tracer populations
    # @param nbins = number of vertical bins

    # set_prof just for the vectors, so that the length can be checked

    def set_prof(self, prof, vec, pop, gp):
        gh.sanitize_vector(vec, len(gp.z_bincenters), -1e30, 1e30, gp.debug)

        #z-profile points
        if prof == 'z_vec':
            self.z_vec = vec

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
        elif prof == 'kz_nu_vec':
            self.kz_nu_vec = vec
        elif prof == 'nu_vec':
            self.nu_vec  = vec
        elif prof == 'sigz2_vec':
            self.sigz2_vec = vec

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
        elif prof == 'kz_nu_vec':
            return self.kz_nu_vec
        elif prof == 'nu_vec':
            return self.nu_vec
        elif prof == 'sigz2_vec':
            return self.sigz2_vec

        ##chi2 of profile
        #elif prof == 'chi2'
        #    return self.chi2

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
