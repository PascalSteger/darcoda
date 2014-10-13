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
    def __init__(self, pops, nipol):
        self.pops = pops
        self.nipol= nipol
        self.x0   = np.zeros(nipol)
        self.chi2 = 0.0
        self.rho  = np.zeros(nipol)
        self.nr   = np.zeros(nipol) #n(r) for mass density
        self.nrnu = np.zeros(nipol) #n(r) for tracer density
        self.M    = np.zeros(nipol)
        self.tilt = np.zeros((pops+1)*nipol) # (pops+1) for overall, 1, 2, ...
        self.nu   = np.zeros((pops+1)*nipol) # dito
        self.sig  = np.zeros((pops+1)*nipol)
        self.Sig  = np.zeros((pops+1)*nipol)
        self.kap  = np.zeros((pops+1)*nipol)
    ## \fn __init__(self, pops, nipol)
    # constructor
    # @param pops number of populations
    # @param nipol number of radial bins


    def set_prof(self, prof, vec, pop, gp):
        gh.sanitize_vector(vec, len(self.x0), -1e30, 1e30, gp.debug)
        if prof == 'rho':
            self.rho = vec
        elif prof == 'nr':
            self.nr = vec
        elif prof == 'M':
            self.M = vec
        elif prof == 'nu':
            self.nu[pop*self.nipol:(pop+1)*self.nipol] = vec
        elif prof == 'Sig':
            self.Sig[pop*self.nipol:(pop+1)*self.nipol] = vec
        elif prof == 'tilt':
            self.tilt[pop*self.nipol:(pop+1)*self.nipol] = vec
        elif prof == 'sig':
            self.sig[pop*self.nipol:(pop+1)*self.nipol] = vec
        elif prof == 'kap':
            self.kap[pop*self.nipol:(pop+1)*self.nipol] = vec
    ## \fn set_prof(self, prof, vec, pop, gp)
    # store density vector
    # @param prof profile identifier
    # @param vec array of floats
    # @param pop population, if applicable
    # @param gp global parameters, used for gp.xepol radii


    def get_prof(self, prof, pop):
        if prof == 'rho':
            return self.rho
        elif prof == 'nr':
            return self.nr
        elif prof == 'M':
            return self.M
        elif prof == 'nu':
            return self.nu[pop*self.nipol:(pop+1)*self.nipol]
        elif prof == 'Sig':
            return self.Sig[pop*self.nipol:(pop+1)*self.nipol]
        elif prof == 'tilt':
            return self.tilt[pop*self.nipol:(pop+1)*self.nipol]
        elif prof == 'sig':
            return self.sig[pop*self.nipol:(pop+1)*self.nipol]
        elif prof == 'kap':
            return self.kap[pop*self.nipol:(pop+1)*self.nipol]
        return self.rho
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
