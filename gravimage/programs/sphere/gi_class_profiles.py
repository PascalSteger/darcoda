#!/usr/bin/env ipython3

##
# @file
# store profiles

# (c) GPL v3 2014 ETHZ psteger@phys.ethz.ch

import pdb
import numpy as np
import gi_helper as gh

class Profiles:
    def __init__(self, pops, nepol):
        self.pops = pops
        self.nepol= nepol
        self.x0   = np.zeros(nepol)
        self.xbins = np.zeros(nepol)
        self.chi2 = 0.0
        self.rho  = np.zeros(nepol)
        self.nr   = np.zeros(nepol)
        self.J = np.zeros(nepol)
        self.M    = np.zeros(nepol)
        self.MtoL = -1.0
        self.betastar = np.zeros((pops+1)*nepol)
        self.beta = np.zeros((pops+1)*nepol) # (pops+1) for overall (rho*), 1, 2, ...
        self.nu   = np.zeros((pops+1)*nepol) # dito
        self.nrnu = np.zeros((pops+1)*nepol)
        self.Sig  = np.zeros((pops+1)*nepol)
        self.sig  = np.zeros((pops+1)*nepol)
        self.kap  = np.zeros((pops+1)*nepol)
        self.zetaa = 0.
        self.zetab = 0.
    ## \fn __init__(self, pops, nepol)
    # constructor
    # @param pops number of populations
    # @param nepol number of radial bins

    def set_prof(self, prof, vec, pop, gp):
        gh.sanitize_vector(vec, len(self.x0), -200, 1e30, gp.debug)
        if prof == 'rho':
            self.rho = vec
        elif prof == 'nr':
            self.nr = vec
        elif prof == 'J':
            self.J = vec
        elif prof == 'M':
            self.M = vec
        elif prof == 'nu':
            self.nu[pop*self.nepol:(pop+1)*self.nepol] = vec # vec has nrho-3 entries
         elif prof == 'nrnu':
            self.nrnu[pop*self.nepol:(pop+1)*self.nepol] = vec
        elif prof == 'betastar':
            self.betastar[pop*self.nepol:(pop+1)*self.nepol] = vec
        elif prof == 'beta':
            self.beta[pop*self.nepol:(pop+1)*self.nepol] = vec
        elif prof == 'Sig':
            self.Sig[pop*self.nepol:(pop+1)*self.nepol] = vec
        elif prof == 'sig':
            self.sig[pop*self.nepol:(pop+1)*self.nepol] = vec
        elif prof == 'kap':
            self.kap[pop*self.nepol:(pop+1)*self.nepol] = vec
        else:
            raise Exception('unknown profile to be set in gi_class_profiles.set_prof')
    ## \fn set_prof(self, prof, vec, pop, gp)
    # store density array
    # @param prof profile identifier
    # @param vec array of floats
    # @param pop population, if applicable
    # @param gp global parameters

    def set_MtoL(self, MtoL):
        self.MtoL = MtoL
        return
    ## \fn set_MtoL(self, MtoL)
    # store mass-to-light ratio
    # @param MtoL float

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
        elif prof == 'J':
            return self.J
        elif prof == 'Sig':
            return self.Sig[pop*self.nepol:(pop+1)*self.nepol]
        elif prof == 'nu':
            return self.nu[pop*self.nepol:(pop+1)*self.nepol]
        elif prof == 'nrnu':
            return self.nrnu[pop*self.nepol:(pop+1)*self.nepol]
        elif prof == 'betastar':
            return self.betastar[pop*self.nepol:(pop+1)*self.nepol]
        elif prof == 'beta':
            return self.beta[pop*self.nepol:(pop+1)*self.nepol]
        elif prof == 'sig':
            return self.sig[pop*self.nepol:(pop+1)*self.nepol]
        elif prof == 'kap':
            return self.kap[pop*self.nepol:(pop+1)*self.nepol]
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
        return "Profiles: "+str(self.pops)+" populations, "+str(self.nepol)+" nepol, "+" chi2 = "+str(self.chi2)
    ## \fn __repr__(self)
    # string representation for ipython


## \class Profiles
# class for storing all profiles
