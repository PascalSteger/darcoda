#!/usr/bin/env ipython3

##
# @file
# store profiles

# (c) 2013 ETHZ psteger@phys.ethz.ch

import pdb
import numpy as np
from gl_project import rho_param_INT_Rho, rho_INTIPOL_Rho
import gl_physics as phys
import gl_helper as gh

class Profiles:
    def __init__(self, pops, nipol):
        self.pops = pops
        self.nipol= nipol
        self.x0   = np.zeros(nipol)
        self.chi2 = 0.0
        self.rho  = np.zeros(nipol)
        self.nr   = np.zeros(nipol)
        self.M    = np.zeros(nipol)
        self.betastar = np.zeros((pops+1)*nipol)
        self.beta = np.zeros((pops+1)*nipol) # (pops+1) for overall (rho*), 1, 2, ...
        self.nu   = np.zeros((pops+1)*nipol) # dito
        self.nrnu = np.zeros((pops+1)*nipol)
        self.Sig  = np.zeros((pops+1)*nipol)
        self.sig  = np.zeros((pops+1)*nipol)
        self.kap  = np.zeros((pops+1)*nipol)
        self.zetaa = np.zeros((pops+1)*nipol)
        self.zetab = np.zeros((pops+1)*nipol)
    ## \fn __init__(self, pops, nipol)
    # constructor
    # @param pops number of populations
    # @param nipol number of radial bins


    def set_prof(self, prof, arr, pop, gp):
        if len(arr) != len(self.x0):
            raise Exception('wrong array to be stored')
        if prof == 'rho':
            self.rho = arr
        elif prof == 'nr':
            self.nr = arr
        elif prof == 'M':
            self.M = arr
        elif prof == 'nu':
            self.nu[pop*self.nipol:(pop+1)*self.nipol] = arr # arr has nrho-3 entries
        elif prof == 'nrnu':
            self.nrnu[pop*self.nipol:(pop+1)*self.nipol] = arr
        elif prof == 'betastar':
            self.betastar[pop*self.nipol:(pop+1)*self.nipol] = arr
        elif prof == 'beta':
            self.beta[pop*self.nipol:(pop+1)*self.nipol] = arr
        elif prof == 'Sig':
            self.Sig[pop*self.nipol:(pop+1)*self.nipol] = arr
        elif prof == 'sig':
            self.sig[pop*self.nipol:(pop+1)*self.nipol] = arr
        elif prof == 'kap':
            self.kap[pop*self.nipol:(pop+1)*self.nipol] = arr
        elif prof == 'zetaa':
            self.zetaa[pop*self.nipol:(pop+1)*self.nipol] = arr
        elif prof == 'zetab':
            self.zetab[pop*self.nipol:(pop+1)*self.nipol] = arr
    ## \fn set_prof(self, prof, arr, pop, gp)
    # store density array
    # @param prof profile identifier
    # @param arr array of floats
    # @param pop population, if applicable
    # @param gp global parameters

        
    def get_prof(self, prof, pop):
        if prof == 'rho':
            return self.rho
        elif prof == 'nr':
            return self.nr
        elif prof == 'M':
            return self.M
        elif prof == 'Sig':
            return self.Sig[pop*self.nipol:(pop+1)*self.nipol]
        elif prof == 'nu':
            return self.nu[pop*self.nipol:(pop+1)*self.nipol]
        elif prof == 'nrnu':
            return self.nrnu[pop*self.nipol:(pop+1)*self.nipol]
        elif prof == 'betastar':
            return self.betastar[pop*self.nipol:(pop+1)*self.nipol]
        elif prof == 'beta':
            return self.beta[pop*self.nipol:(pop+1)*self.nipol]
        elif prof == 'sig':
            return self.sig[pop*self.nipol:(pop+1)*self.nipol]
        elif prof == 'kap':
            return self.kap[pop*self.nipol:(pop+1)*self.nipol]
        elif prof == 'zetaa':
            return self.zetaa[pop*self.nipol:(pop+1)*self.nipol]
        elif prof == 'zetab':
            return self.zetab[pop*self.nipol:(pop+1)*self.nipol]
        return self.rho
    ## \fn get_prof(self, prof, pop)
    # return density array
    # @param prof of this profile
    # @param pop and of this population (if applicable)


    def __repr__(self):
        return "Profiles: "+str(self.pops)+" populations, "+str(self.nipol)+" nipol, "+" chi2 = "+str(self.chi2)
    ## \fn __repr__(self)
    # string representation for ipython


## \class Profiles
# class for storing all profiles

