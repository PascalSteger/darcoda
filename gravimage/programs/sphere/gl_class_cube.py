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

def ginv(x, mu, sig):
    # TODO
    return x
## \fn ginv(x, mu, sig)
# take range [0,1] and return [-inf, inf]
# sampled with a normal distribution


def map_nr(params, prof, pop, gp):
    gh.sanitize_vector(params, gp.nrho, 0, 1, gp.debug)
    nr = np.zeros(gp.nepol) # to hold the n(r) = dlog(rho)/dlog(r) values

    # get offset and n(r) profiles, calculate rho
    if prof=='rho':
        rhoscale = gp.rhohalf
        Rscale = gp.Xscale[0]
        width = gp.log10rhospread
        rlimnr = gp.rlimnr
        innerslope = gp.innerslope
        nrscale = gp.nrtol/(max(np.log(gp.xipol))-min(np.log(gp.xipol)))
        monotonic = gp.monotonic
    elif prof=='nu':
        rhoscale = gp.dat.nuhalf[pop]
        Rscale = gp.Xscale[pop]
        width = gp.log10nuspread
        rlimnr = gp.rlimnr_nu
        innerslope = gp.innerslope
        nrscale = gp.nrtol_nu/(max(np.log(gp.xipol))-min(np.log(gp.xipol)))
        monotonic = gp.monotonic_nu
    else:
        raise Exception('wrong profile in gl_class_cube.map_nr')

    # first parameter gives half-light radius value of rho directly
    # use [0,1]**3 to increase probability of sampling close to 0
    # fix value with tracer densities,
    # sample a flat distribution over log(rho_half)
    rhohalf = 10**((params[0]-0.5)*2.*width+np.log10(rhoscale))

    # nr(r=0) is = rho slope for approaching r=0 asymptotically, given directly
    # should be smaller than -3 to exclude infinite enclosed mass
    nrasym0 = params[1]*innerslope

    # work directly with the dn(r)/dlog(r) parameters here
    dnrdlrparams = params[2:]
    # offset for the integration of dn(r)/dlog(r) at smallest radius
    nr[0] = dnrdlrparams[0]*innerslope

    for k in range(1, gp.nepol):
        deltalogr = (np.log(gp.xepol[k-1])-np.log(gp.xepol[k-2]))
        # construct n(r_k+1) from n(r_k)+dn/dlogr*Delta log r, integrated
        if monotonic:
            # only increase n(r), use pa[i]>=0 directly
            nr[k] = nr[k-1] + dnrdlrparams[k] * nrscale/2. * deltalogr
        else:
            # use pa => [-1, 1] for full interval
            nr[k] = nr[k-1] + (dnrdlrparams[k]-0.5)*2. * nrscale * deltalogr

        # cut at zero: we do not want to have density rising outwards
        nr[k] = max(0., nr[k])

        # restrict n(r)
        if gp.xepol[k] <= rlimnr*Rscale:
            nr[k] = min(maxrhoslope/2, nr[k])
        else:
            nr[k] = min(maxrhoslope, nr[k])
    # rho slope for asymptotically reaching r = \infty is given directly
    # must lie below -3, thus n(r)>3
    deltalogrlast = (np.log(gp.xepol[-1])-np.log(gp.xepol[-2]))
    # to ensure we have a finite mass at all radii 0<r<=\infty
    if monotonic:
        nrasyminfty = nr[-1]+params[-1] * nrscale/2. * deltalogrlast
    else:
        nrasyminfty = nr[-1]+(params[-1]-0.5)*2 * nrscale * deltalogrlast

    # finite mass prior: to bound between 3 and gp.maxrhoslope, favoring 3:
    nrasyminfty = max(nrasyminfty, 3.001)

    params = np.hstack([rhohalf, nrasym0, nr, nrasyminfty])
    return params
## \fn map_nr(params, prof, pop, gp)
# mapping rho parameters from [0,1] to full parameter space
# setting all n(r<r_nrlim)<=2.0
# and possibly a monotonically increasing function
# first parameter is offset for rho_half
# second parameter is asymptotic n(r to 0) value
# @param params cube [0,1]^ndim
# @param prof string nu, rho, rhostar
# @param pop population int, 0 for rho*, 1,2,... for tracer densities
# @param gp global parameters


def map_betastar_poly(params, gp):
    gh.sanitize_vector(params, gp.nbeta, 0, 1, gp.debug)
    off_beta = 0
    # beta* parameters : [0,1] ->  some range, e.g. [-1,1]
    # starting offset in range [-1,1]
    # cluster around 0, go symmetrically in both directions,
    params[0] = 1.98*(params[0]-0.5)
    # out to maxbetaslope
    # here we allow |beta_star,0| > 1, so that any models with
    # beta(<r_i) = 1, beta(>r_i) < 1
    # are searched as well
    # pa[0] = np.sign(tmp)*tmp**2 # between -1 and 1 for first parameter
    off_beta += 1
    for i in range(gp.nbeta-1):
        params[off_beta] = (2*(params[off_beta]-0.5))*gp.maxbetaslope
        # rising beta prior would remove -0.5
        off_beta += 1

    return params
## \fn map_betastar_old(params, gp)
# mapping beta parameters from [0,1] to full parameter space,
# using consecutive polynomials
# NOT USED ANYMORE
# @param params parameter array
# @param gp global parameters


def map_betastar_sigmoid(params, gp):
    gh.sanitize_vector(params, gp.nbeta, 0, 1, gp.debug)
    # s0 = np.log(r0/r0turn)
    # kappa = (a0-a1)/(betastar(r_s) - a1)-1
    # beta = (a0-a1)/(1+kappa*exp(alpha*s0))
    bmin = gp.minbetastar
    bmax = gp.maxbetastar
    bdiff = bmax-bmin
    a0 = params[0]*bdiff + bmin  # a0
    if gp.beta00prior:
        params[0] = 0
    a1 = params[1]*bdiff + bmin  # a1
    alpha = params[2]*5             # alpha
    # r_s, sampled in log space over all radii,
    # as we want flat prior in log space
    logrs = params[3]*(np.log(max(gp.xepol))-np.log(min(gp.xepol)))+np.log(min(gp.xepol))
    return np.hstack([a0, a1, alpha, logrs])
## \fn map_betastar(pa, gp)
# mapping beta parameters from [0,1] to full param space
# @param pa parameter vector, size 4
# @param gp global parameters


def map_betastar_j(params, gp):
    gh.sanitize_vector(params, 4, 0, 1, gp.debug)
    # betastar = exp(-(r/r0)^n)*(a0-a1)+a1
    a0 = params[0]*1.98-1 # a_0, betastar(r=0), in between -0.99 and +0.99
    a1 = params[1]*1.98-1 # a_1, betastar(r->infty), same range
    # r_0, scale radius for transition from a0->a1:
    logrs = params[2]*(np.log(max(gp.xepol))-np.log(min(gp.xepol)))+np.log(min(gp.xepol))
    n0 = params[3]*3 # n, rate of transition
    return np.hstack([a0, a1, logrs, n0])
## \fn map_betastar(params, gp)
# mapping beta parameters from [0,1] to full param space
# @param params parameter vector, size 4
# @param gp global parameters


def map_MtoL(param, gp):
    gh.sanitize_scalar(param, 0, 1, gp.debug)
    scale = gp.MtoLmax - gp.MtoLmin
    MtoL = param*scale+gp.MtoLmin
    return MtoL
## \fn map_MtoL(param, gp)
# map [0,1] to MtoL flat prior
# @param param scalar
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
