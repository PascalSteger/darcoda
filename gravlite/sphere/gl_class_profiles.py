#!/usr/bin/env python3

##
# @file
# store profiles

# (c) 2013 ETHZ psteger@phys.ethz.ch

class Profiles:
    def __init__(self, pops, nipol):
        self.pops = pops
        self.nipol= nipol
        import numpy as np
        self.rho  = np.zeros(nipol)
        self.M    = np.zeros(nipol)
        self.beta = np.zeros((pops+1)*nipol) # (pops+1) for overall, 1, 2, ...
        self.nu   = np.zeros((pops+1)*nipol) # dito
        self.sig  = np.zeros((pops+1)*nipol)
        self.kap  = np.zeros((pops+1)*nipol)
    ## \fn __init__(self, pops, nipol)
    # constructor
    # @param pops number of populations
    # @param nipol number of radial bins

        
    def set_rho(self, arr):
        self.rho = arr
    ## \fn set_rho(self, arr)
    # store density array
    # @param arr array of floats

        
    def get_rho(self):
        return self.rho
    ## \fn get_rho(self)
    # return density array
        
        
    def set_M(self, arr):
        self.M = arr
    ## \fn set_M(self, arr)
    # store mass array
    # @param arr float array


    def get_M(self):
        return self.M
    ## \fn get_M(self)
    # return mass profile in [Msun/pc^3]

    
    def set_nu(self, pop, arr):
        self.nu[pop*self.nipol:(pop+1)*self.nipol] = arr
    ## \fn set_nu(self, pop, arr)
    # store tracer density array, at right position
    # @param pop population number, starts at 0 for all comp, 1 for first, ...
    # @param arr float of profile

    
    def get_nu(self, pop):
        return self.nu[pop*self.nipol:(pop+1)*self.nipol]
    ## \fn get_nu(self, pop):
    # return tracer density
    # @param pop population number

    
    def set_beta(self, pop, arr):
        self.beta[pop*self.nipol:(pop+1)*self.nipol] = arr
    ## \fn set_beta(self, pop, arr)
    # store velocity anisotropy array
    # @param pop population number, starts at 0
    # @param arr array of floats


    def get_beta(self, pop):
        return self.beta[pop*self.nipol:(pop+1)*self.nipol]
    ## \fn get_beta(self, pop)
    # get stored beta profile
    # @param pop population number

    
    def set_sig_kap(self, pop, sigarr, kaparr):
        self.sig[pop*self.nipol:(pop+1)*self.nipol] = sigarr
        self.kap[pop*self.nipol:(pop+1)*self.nipol] = kaparr
    ## \fn set_sig_kap(self, pop, sigarr, kaparr)
    # store velocity dispersion, fourth order moment profile
    # @param pop population number, starts at 0
    # @param sigarr array of floats
    # @param kaparr array of floats

    
    def get_sig(self, pop):
        return self.sig[pop*self.nipol:(pop+1)*self.nipol]
    ## \fn get_sig(self, pop)
    # return sigma of a given population
    # @param pop population number, start at 0 for first component

    
    def get_kap(self, pop):
        return self.kap[pop*self.nipol:(pop+1)*self.nipol]
    ## \fn get_kap(self, pop)
    # return kappa of a given population
    # @param pop population number, start at 0 for first component

## \class Profiles
# class for storing all profiles
