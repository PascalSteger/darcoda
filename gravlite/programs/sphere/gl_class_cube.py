#!/usr/bin/env ipython3

##
# @file
# parameters for MultiNest.
# cube = [0,1]^ndim
# Holds representations for overall density,
# [tracer density, anisotropy]_{for each population}
# Gives access to density, population profiles,
# and calculate physical values from [0,1]
# populations are counted from 0 = first component,
# 1 = first additional component

# (c) 2013 ETHZ Pascal S.P. Steger

import numpy as np
import pdb
import gl_helper as gh


def map_nr(pa, prof, pop, gp):
    # get offset and n(r) profiles, calculate rho
    if prof=='rho':
        scale = gp.rhohalf
        width = gp.rhospread
        iscale = gp.iscale
        maxrhoslope = gp.maxrhoslope
        nrscale = gp.nrtol
        monotonic = gp.monotonic
    elif prof=='nu':
        scale = gp.dat.nuhalf[pop]
        width = gp.nuspread
        iscale = gp.iscale_nu
        maxrhoslope = gp.maxrhoslope_nu
        nrscale = gp.nrtol_nu
        monotonic = gp.monotonic_nu
    else:
        raise Exception('bad profile in gl_class_cube.map_nr')

    # first parameter gives half-light radius value of rho directly
    # use [0,1]**3 to increase probability of sampling close to 0
    # fix value with tracer densities,
    # sample a flat distribution over log(rho_half)
    pa[0] = 10**((pa[0]*2.*width)-width+np.log10(scale))

    # nr(r=0) is = rho slope for approaching r=0 asymptotically, given directly
    # should be smaller than -3 to exclude infinite enclosed mass
    if 1 <= iscale + 1:
        pa[1] = (pa[1]**1)*2.0
    else:
        pa[1] = (pa[1]**1)*2.999

    # offset for the integration of dn(r)/dlog(r) at smallest radius
    if 2 <= iscale + 1:
        pa[2] = (pa[2]**1)*2.0
    else:
        pa[2] = (pa[2]**1)*maxrhoslope

    rdef = gp.xepol # [pc]
    for i in range(3, gp.nrho-1):
        # all -dlog(rho)/dlog(r) at data points and 2,4,8rmax can
        # lie in between 0 and gp.maxrhoslope
        if monotonic:
            # only increase n(r), use pa[i]>=0 directly
            pa[i] = pa[i-1]+pa[i]*nrscale*(np.log(rdef[i-2])-np.log(rdef[i-3]))
        else:
            # use pa => [-1, 1] for full interval
            pa[i] = pa[i-1]+(pa[i]-0.5)*2. * nrscale * (np.log(rdef[i-2])-np.log(rdef[i-3]))

        pa[i] = max(0., pa[i])
        if i <= iscale+1: # iscale: no. bins with xipol<Rscale[all]
            pa[i] = min(2.0, pa[i])
        else:
            pa[i] = min(maxrhoslope, pa[i])
    # rho slope for asymptotically reaching r = \infty is given directly
    # must lie below -3
    # to ensure we have a finite mass at all radii 0<r<=\infty
    pa[gp.nrho-1] = pa[gp.nrho-1] * maxrhoslope
    if monotonic:
        pa[gp.nrho-1] += pa[gp.nrho-2]
    # finite mass prior: to bound between 3 and gp.maxrhoslope, favoring 3:
    # pa[gp.nrho-1] = max(pa[gp.nrho-1], 3.)
    return pa
## \fn map_nr(pa, prof, pop, gp)
# mapping nr and rho parameters from [0,1] to full parameter space
# setting all n(r<r_{iscale})<=2.0
# and possibly a monotonically increasing function
# first parameter is offset for rho_half
# second parameter is asymptotic n(r to 0) value
# @param pa cube [0,1]^ndim
# @param prof string nu, rho, rhostar
# @param pop population int, 0 for rho*, 1,2,... for tracer densities
# @param gp global parameters


def map_nu_directly(pa, gp):
    for i in range(gp.nepol):
        pa[i] = 10**(pa[i]*(gp.maxlog10nu-gp.minlog10nu)+gp.minlog10nu)
    return pa
## \fn map_nu_directly(pa, gp)
# map tracer densities, directly
# @param pa cube [0,1]^n
# @param gp global parameters


def map_betastar_old(pa, gp):
    off_beta = 0
    # beta* parameters : [0,1] ->  some range, e.g. [-1,1]
    # starting offset in range [-1,1]
    # cluster around 0, go symmetrically in both directions,
    pa[0] = 1.98*(pa[0]-0.5)
    # out to maxbetaslope
    # here we allow |beta_star,0| > 1, so that any models with
    # beta(<r_i) = 1, beta(>r_i) < 1
    # are searched as well
    # pa[0] = np.sign(tmp)*tmp**2 # between -1 and 1 for first parameter
    off_beta += 1
    for i in range(gp.nbeta-1):
        pa[off_beta] = (2*(pa[off_beta]-0.5))*gp.maxbetaslope
        # rising beta prior would remove -0.5
        off_beta += 1

    return pa
## \fn map_betastar_old(pa, gp)
# mapping beta parameters from [0,1] to full parameter space,
# using consecutive polynomials
# NOT USED ANYMORE
# @param pa parameter array
# @param gp global parameters


def map_betastar_sigmoid(pa, gp):
    gh.sanitize_vector(pa, 5, 0, 1)
    # s0 = np.log(r0/r0turn)
    # kappa = (a0-a1)/(betastar(r_s) - a1)-1
    # beta = (a0-a1)/(1+kappa*exp(alpha*s0))
    bmin = gp.minbetastar
    bmax = gp.maxbetastar
    bdiff = bmax-bmin
    pa[0] = pa[0]*bdiff + bmin  # a0
    if gp.beta00prior:
        pa[0] = 0
    pa[1] = pa[1]*bdiff + bmin  # a1
    pa[2] = pa[2]*5             # alpha
    pa[3] = pa[3]*max(gp.xepol) # r_s
    return pa
## \fn map_betastar(pa, gp)
# mapping beta parameters from [0,1] to full param space
# @param pa parameter vector, size 4
# @param gp global parameters


def map_betastar_j(pa, gp):
    gh.sanitize_vector(pa, 4, 0, 1)
    # betastar = exp(-(r/r0)^n)*(a0-a1)+a1
    pa[0] = pa[0]*1.98-1 # a_0, betastar(r=0), in between -0.99 and +0.99
    pa[1] = pa[1]*1.98-1 # a_1, betastar(r->infty), same range
    pa[2] = pa[2]*max(gp.xipol) # r_0, scale radius for transition from a0->a1
    pa[3] = pa[3]*3 # n, rate of transition
    return pa
## \fn map_betastar(pa, gp)
# mapping beta parameters from [0,1] to full param space
# @param pa parameter vector, size 4
# @param gp global parameters


def map_MtoL(pa, gp):
    gh.sanitize_scalar(pa, 0, 1)
    scale = gp.MtoLmax - gp.MtoLmin
    pa = pa*scale+gp.MtoLmin
    return pa
## \fn map_MtoL(pa, gp)
# map [0,1] to MtoL flat prior
# @param pa scalar
# @param gp global parameters holding MtoL{min,max}


class Cube:
    def __init__ (self, gp):
        self.pops = gp.pops
        # for density and (nu, beta)_i
        self.cube = np.zeros(gp.ndim)
        return
    ## \fn __init__ (self, gp)
    # constructor, with modes depending on locpop
    # @param gp


    def convert_to_parameter_space(self, gp):
        # if we want any priors, here they have to enter:
        off = 0
        pc = self.cube
        # DM density rho, set in parametrization of n(r)
        offstep = gp.nrho
        tmp_nr = map_nr(pc[0:offstep], 'rho', 0, gp)
        for i in range(offstep):
            pc[off+i] = tmp_nr[i]
        off += offstep

        # rho* only for observations
        if gp.investigate == 'obs':
            offstep = gp.nrho
            tmp_rhostar = map_nr(pc[off:off+offstep], 'nu', 0, gp)
            for i in range(offstep):
                pc[off+i] = tmp_rhostar[i]
            off += offstep

            offstep = 1
            pc[off] = map_MtoL(pc[off], gp)
            off += offstep

        for pop in range(1,gp.pops+1): # nu1, nu2, ...
            offstep = gp.nrho
            tmp_nu = map_nr(pc[off:off+offstep], 'nu', pop, gp)
            for i in range(offstep):
                pc[off+i] = tmp_nu[i]
            off += offstep

            offstep = gp.nbeta
            tmp_betastar = map_betastar_sigmoid(pc[off:off+offstep], gp)
            for i in range(offstep):
                pc[off+i] = tmp_betastar[i]
            off += offstep

        if off != gp.ndim:
            gh.LOG(1,'wrong subscripts in gl_class_cube')
            raise Exception('wrong subscripts in gl_class_cube')

        return pc
    ## \fn convert_to_parameter_space(self, gp)
    # convert [0,1]^ndim to parameter space
    # such that values in cube are the parameters we need for rho, nu_i, beta_i
    # @param gp


    def __repr__(self):
        return "Cube: "+self.pops+" populations"
    ## \fn __repr__(self)
    # string representation for ipython


    def copy(self, cub):
        self.cube = cub
        return self
    ## \fn copy(self, cub)
    # copy constructor

## \class Cube
# Common base class for all parameter sets
