#!/usr/bin/env python3

##
# @file
# parameters for MultiNest.
# cube = [0,1]^ndim,  ndim = nipol + 2*pops*nipol, 
# Holds representations for overall density,
# [tracer density, anisotropy]_{for each population}
# Gives access to density, population profiles,
# and calculate physical values from [0,1]
# populations are counted from 0 = first component,
# 1 = first additional component

# (c) 2013 ETHZ Pascal S.P. Steger

import pdb
import numpy as np
import numpy.random as npr
import random
import gl_params as gp
import gl_physics as phys
import gl_helper as gh

def nu_offset(pop):
    return gp.nipol + pop*gp.nipol + pop*gp.nbeta
## \fn nu_offset(pop)
# determine offset for tracer density profiles in data cube
# @param pop which population (0=first, 1=first additional, ...)

    
def beta_offset(pop):
    return gp.nipol + pop*gp.nipol + pop*gp.nbeta + gp.nipol
## \fn beta_offset(pop)
# determine offset for anisotropy profiles in data cube
# @param pop which population (0=first, 1=first additional, ...)
# @param nipol number of radial bins


def mapping(arr, minrange, maxrange):
    tmp = arr * (np.log10(maxrange)-np.log10(minrange))
    tmp = tmp + np.log10(minrange)
    tmp = 10.**tmp
    return tmp
## \fn mapping(arr, minrange, maxrange)
# scale [0,1] to [minrange, maxrange]
# @param arr array of floats
# @param minrange float of minimum
# @parma maxrange float of maximum


def mapping_beta_slope(arr):
    beta = np.zeros(len(arr))
    beta[0] = 2.*arr[0]-1. # starting value uniformly distributed between -1 and 1
    for i in np.arange(1,len(arr)):
        # TODO bounds and radius scaling
        beta[i] = beta[i-1] + \
          (arr[i]-0.5)*np.sqrt((gp.xipol[i]-gp.xipol[i-1])/(gp.xipol[-1]))
    return beta
## \fn mapping_beta_slope(arr):
# map [0, 1] to [-infty, 1] for each
# @param arr float array with values in [0,1]


def mapping_beta_lognorm(arr):
    from scipy.stats import lognorm
    return 1.-lognorm.ppf(1.-np.array(arr),0.7)
## \fn mapping_beta_lognorm(arr)
# map [0, 1] to [-infty, 1] for each
# @param arr float array with values in [0,1]


def mapping_beta_log(arr):
    return np.log(np.array(arr)**0.5) + 1.
## \fn mapping_beta_log(arr)
# map [0, 1] to [-infty, 1] for each
# @param arr float array with values in [0,1]


class Cube:
    def __init__ (self, pops):
        self.pops = pops
        # for density and (nu, beta)_i
        self.ndim = gp.nipol + pops*gp.nipol + pops*gp.nbeta
        self.cube = np.zeros(self.ndim)
        return
    ## \fn __init__ (self, pops)
    # constructor, with modes depending on locpop
    
    def copy(self, cub):
        self.cube = cub
        return self
    ## \fn copy(self, cub)
    # copy constructor

    def get_nu(self, pop):
        off = nu_offset(pop)
        return phys.rho(gp.xipol, self.cube[off:off+gp.nipol])
    ## \fn get_nu(self,pop)
    # return tracer densities
    
    def set_nu(self, pop, newnu):
        off = nu_offset(pop)
        for i in range(nipol):
            self.cube[off+i] = newnu[i]
        return self.cube
    ## \fn set_nu(self, pop, newnu)
    # set all tracer densities


    def get_beta(self, pop):
        off = beta_offset(pop)
        return phys.beta(gp.xipol, self.cube[off:off+gp.nbeta])
    ## \fn get_beta(self, pop)
    # return back all tracer densities

    
    def set_beta(self, pop, newbeta):
        off = beta_offset(pop)
        for i in range(gp.nbeta):
            self.cube[off+i] = newbeta[i]
        return self.cube
    ## \fn set_beta(self,pop,newbeta)
    # set all tracer densities

    def convert_to_parameter_space(self):
        # if we want any priors, here they have to enter:
        pc = self.cube
        # rho: density parameters -> adjust halflight density
        # use [0,1]**3 to increase probability of sampling close to 0
        pc[0] = 1./(1-pc[0]**3)-1. # [0,1[ => [0, infty[

        # rho slope for approaching r=0 asymptotically should be smaller than -3
        # to exclude infinite enclosed mass
        pc[1] = (pc[1]**3)*2.999
        for i in range(2, gp.nepol-1):
            # all   -dlog(rho)/dlog(r) at data points and 2,4,8rmax can lie in between 0 and gp.maxrhoslope
            pc[i] *= gp.maxrhoslope
        # rho slope for asymptotically reaching r = \infty must lie below -3
        # to ensure we have a finite mass at all radii 0<r<\infty
        pc[gp.nepol-1] = pc[gp.nepol-1] * gp.maxrhoslope + 3.
        off = gp.nepol
        for pop in range(gp.pops):
            # nu parameters -> adjust density at halflight radius
            # nu(r_half)
            pc[off] = 1.-pc[off]**10           # [0,1[, skewed towards 1
            pc[off] = 1./pc[off]-1.          # [0,1[ => [0, infty[
            # nu inner asymptote
            pc[off+1] = (pc[off+1]**3)*2.999 # see above
            for i in range(2, gp.nepol-1):
                pc[off+i] *= gp.maxnuslope # see above
            # pc[off+gp.nepol-1] = (pc[off+gp.nepol-1]**3) * 3. + 3. # see above
            pc[off+gp.nepol-1] = pc[off+gp.nepol-1] * gp.maxnuslope + 3.
            off += gp.nepol

            # beta* parameters : [0,1] ->  some range, e.g. [-1,1]
            # starting offset in range [-1,1]
            # pc[off] = 2.*np.array(pc[off])-1.
            # betatmp = np.array(pc[off:off+gp.nbeta]) # 1:1, only take [0, 1] (<=> rising beta prior)

            tmp = 2*(pc[off]-0.5)
            # cluster around 0, go symmetrically into both directions,
            # out to maxbetaslope
            # here we allow |beta_star,0| > 1, so that any models with
            # beta(<r_i) = 1, beta(>r_i) < 1
            # are searched as well
            pc[off] = np.sign(tmp)*tmp**2 # between -1 and 1 for first parameter
            off += 1

            tmp = pc[off] # [-1,1]
            pc[off] = tmp*gp.maxbetaslope # rising beta prior
            off += 1
            
        return pc
    ## \fn convert_to_parameter_space(self)
    # convert [0,1] to parameter space
    # such that values in cube are the parameters we need for phys.nu/rho/beta
    
## \class Cube
# Common base class for all parameter sets
