#!/usr/bin/env python3

##
# @file
# store profiles

# (c) GPL v3 2014 ETHZ psteger@phys.ethz.ch

import pdb
import numpy as np
import gl_helper as gh

class Profiles:
    def __init__(self, pops, nipol):
        self.pops = pops
        self.nipol= nipol
        self.x0   = np.zeros(nipol)
        self.xbins = np.zeros(nipol)
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
        self.zetaa = 0.
        self.zetab = 0.
    ## \fn __init__(self, pops, nipol)
    # constructor
    # @param pops number of populations
    # @param nipol number of radial bins

    def set_prof(self, prof, vec, pop, gp):
        gh.sanitize_vector(vec, len(self.x0), -200, 1e30, gp.debug)
        if prof == 'rho':
            self.rho = vec
        elif prof == 'nr':
            self.nr = vec
        elif prof == 'M':
            self.M = vec
        elif prof == 'nu':
            self.nu[pop*self.nipol:(pop+1)*self.nipol] = vec # vec has nrho-3 entries
        elif prof == 'nrnu':
            self.nrnu[pop*self.nipol:(pop+1)*self.nipol] = vec
        elif prof == 'betastar':
            self.betastar[pop*self.nipol:(pop+1)*self.nipol] = vec
        elif prof == 'beta':
            self.beta[pop*self.nipol:(pop+1)*self.nipol] = vec
        elif prof == 'Sig':
            self.Sig[pop*self.nipol:(pop+1)*self.nipol] = vec
        elif prof == 'sig':
            self.sig[pop*self.nipol:(pop+1)*self.nipol] = vec
        elif prof == 'kap':
            self.kap[pop*self.nipol:(pop+1)*self.nipol] = vec
        else:
            raise Exception('unknown profile to be set in gl_class_profiles.set_prof')
    ## \fn set_prof(self, prof, vec, pop, gp)
    # store density array
    # @param prof profile identifier
    # @param vec array of floats
    # @param pop population, if applicable
    # @param gp global parameters


    def set_zeta(self, zetaa, zetab, pop):
        self.zetaa = zetaa
        self.zetab = zetab
    ## \fn set_zeta(self, zetaa, zetab, pop)
    # store zeta scalars, RichardsonFairbairn 2014
    # @param zetaa first virial parameter
    # @param zetab second virial parameter
    # @param pop 1,2,...


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
        else:
            raise Exception('specify prof!')
    ## \fn get_prof(self, prof, pop)
    # return density array
    # @param prof of this profile
    # @param pop and of this population (if applicable)


    def get_zeta(self, pop):
        return self.zetaa, self.zetab
    ## \fn get_zeta(self, pop)
    # return zeta parameters, scalar
    # @param pop int 1,2,... for population


    def __repr__(self):
        return "Profiles: "+str(self.pops)+" populations, "+str(self.nipol)+" nipol, "+" chi2 = "+str(self.chi2)
    ## \fn __repr__(self)
    # string representation for ipython


## \class Profiles
# class for storing all profiles
