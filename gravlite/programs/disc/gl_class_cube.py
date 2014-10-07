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
# disc version, done

# (c) 2013 ETHZ Pascal S.P. Steger

# SS testing (7 october 2014)

import ipdb
import numpy as np
import numpy.random as npr
import gl_helper as gh

def map_tilt_slope(vec, gp):
    T = np.zeros(len(vec))
    T[0] = vec[0] # uniformly between 0 and 1
    for i in np.arange(1,len(vec)):
        # TODO bounds and radius scaling
        T[i] = T[i-1] + \
          (vec[i]-0.5)*np.sqrt((gp.xipol[i]-gp.xipol[i-1])/(gp.xipol[-1]))
    return T
## \fn map_tilt_slope(vec, gp):
# map [0, 1] to [-infty, 1] for each
# @param vec float array with values in [0,1]
# @param gp


def map_tiltstar(pa, gp):
    off = 0
    # tilt parameters : [0,1] ->  some range, e.g. [-1,1]
    # starting offset in range [-1,1]
    # cluster around 0, go symmetrically in both directions,
    tmp = 2*(pa[0]-0.5)
    pa[0] = np.sign(tmp)*tmp**2 # between -1 and 1 for first parameter
    off += 1
    for i in range(gp.nbeta-1):
        pa[off] = ((pa[off]-0.5)*2.)*gp.maxbetaslope
        off += 1
    return pa
## \fn map_tiltstar(pa)
# mapping tilt parameters from [0,1] to full parameter space
# @param pa parameter array


def map_nr(params, prof, pop, gp):
    gh.sanitize_vector(params, gp.nrho, 0, 1)
    nr = np.zeros(gp.nepol) # to hold the n(r) = dlog(rho)/dlog(r) values

    # get offset and n(r) profiles, calculate rho
    if prof=='rho':
        rhoscale = gp.rhohalf
        Rscale = gp.Xscale[0]
        width = gp.log10rhospread
        rlimnr = gp.rlimnr
        maxrhoslope = gp.maxrhoslope
        nrscale = gp.nrtol/(max(np.log(gp.xipol))-min(np.log(gp.xipol)))
        monotonic = gp.monotonic
    elif prof=='nu':
        rhoscale = gp.dat.nuhalf[pop]
        Rscale = gp.Xscale[pop]
        width = gp.nuspread
        rlimnr = gp.rlimnr_nu
        maxrhoslope = gp.maxrhoslope_nu
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
    if gp.xepol[0] <= rlimnr*Rscale:
        nrasym0 = (params[1]**1)*min(maxrhoslope/2, 2.99)
    else:
        nrasym0 = (params[1]**1)*2.99

    # work directly with the dn(r)/dlog(r) parameters here
    dnrdlrparams = params[2:-1]
    # offset for the integration of dn(r)/dlog(r) at smallest radius
    if gp.xepol[1] <= rlimnr*Rscale:
        nr[0] = (dnrdlrparams[0]**1)*min(maxrhoslope/2, 2.99)
    else:
        nr[0] = (dnrdlrparams[0]**1)*maxrhoslope

    for k in range(1, gp.nepol):
        # all -dlog(rho)/dlog(r) at data points and 2,4,8rmax can
        # lie in between 0 and gp.maxrhoslope
        deltalogr = (np.log(gp.xepol[k-1])-np.log(gp.xepol[k-2]))
        # construct n(r_k+1) from n(r_k)+dn/dlogr*Delta log r, integrated
        if monotonic:
            # only increase n(r), use pa[i]>=0 directly
            nr[k] = nr[k-1] + dnrdlrparams[k] * nrscale * deltalogr
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
        nrasyminfty = nr[-1]+params[-1] * nrscale * deltalogrlast
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


def map_nu(pa, gp):
    # TODO: assertion len(pa)=gp.nepol
    for i in range(len(pa)):
        pa[i] = 10**(pa[i]*(gp.maxlog10nu-gp.minlog10nu)+gp.minlog10nu)
    return pa
## \fn map_nu(pa, gp)
# map tracer densities, directly
# @param pa cube [0,1]^n
# @param gp global parameters


def map_MtoL(pa, gp):
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
        # for density and (nu, tilt)_i
        self.cube = np.zeros(gp.ndim)
        return
    ## \fn __init__ (self, pops)
    # constructor, with modes depending on locpop


    def convert_to_parameter_space(self, gp):
        # if we want any priors, here they have to enter:
        pc = self.cube
        off = 0
        offstep = 1
        pc[off] = pc[off]*200-100+17**2 # for normalization C
        off += offstep

        offstep = gp.nrho
        tmp = map_nr(pc[off:off+offstep], 'rho', 0, gp)
        for i in range(offstep):
            pc[off+i] = tmp[i]
        off += offstep

        # rho_baryons
        offstep = gp.nrho
        tmp_rho_baryons = map_nr(pc[off:off+offstep], 'nu', 0, gp)
        for i in range(offstep):
            pc[off+i] = tmp_rho_baryons[i]
        off += offstep

        offstep = 1
        pc[off] = map_MtoL(pc[off], gp)
        off += offstep

        for pop in range(1, gp.pops+1):
            offstep = gp.nrho
            tmp = map_nr(pc[off:off+offstep], 'nu', pop, gp)
            for i in range(offstep):
                pc[off+i] = tmp[i]
            off += offstep

            offstep = gp.nbeta
            tmp = map_tiltstar(pc[off:off+offstep], gp)
            for i in range(offstep):
                pc[off+i] = tmp[i]
            off += offstep

        if off != gp.ndim:
            gh.LOG(1,'wrong subscripts in gl_class_cube')
            raise Exception('wrong subscripts in gl_class_cube')

        return pc
    ## \fn convert_to_parameter_space(self, gp)
    # convert [0,1]^ndim to parameter space
    # such that values in cube are the parameters we need for rho, nu_i, beta_i
    # @param gp global parameters


    def __repr__(self):
        return "Cube (disc) with "+str(gp.pops)+" pops "

    def copy(self, cub):
        self.cube = cub
        return self
    ## \fn copy(self, cub)
    # copy constructor


## \class Cube
# Common base class for all parameter sets
