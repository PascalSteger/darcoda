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

import pdb
import numpy as np
import numpy.random as npr

def map_tilt_slope(arr, gp):
    T = np.zeros(len(arr))
    T[0] = arr[0] # uniformly between 0 and 1
    for i in np.arange(1,len(arr)):
        # TODO bounds and radius scaling
        T[i] = T[i-1] + \
          (arr[i]-0.5)*np.sqrt((gp.xipol[i]-gp.xipol[i-1])/(gp.xipol[-1]))
    return T
## \fn map_tilt_slope(arr, gp):
# map [0, 1] to [-infty, 1] for each
# @param arr float array with values in [0,1]
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


def map_nr(pa, prof, pop, gp):
    # first parameter gives half-light radius value of rho directly
    pa[0] = 10**((pa[0]*2.*width)-width+np.log10(scale))
    if prof=='rho':
        maxrhoslope = gp.maxrhoslope
        monotonic = gp.monotonic
        iscale = gp.iscale
        nrtol = gp.nrtol
        scale = gp.rhohalf
        width = gp.rhospread
    elif prof == 'nu':
        maxrhoslope = gp.maxrhoslope_nu
        monotonic = gp.monotonic_nu
        iscale = gp.iscale_nu
        nrtol = gp.nrtol_nu
        scale = gp.dat.nuhalf[pop]
        width = gp.nuspread[pop]
    else:
        raise Exception('wrong prof in gl_class_cube.map_nr')

    # nr(r=0) is = rho slope for approaching r=0 asymptotically, given directly
    # should be smaller than -3 to exclude infinite enclosed mass
    if 1 <= iscale + 1:
        pa[1] = (pa[1]**1)*maxrhoslope/2.
    else:
        pa[1] = (pa[1]**1)*min(maxrhoslope, 2.99)

    # offset for the integration of dn(r)/dlog(r) at smallest radius
    if 2 <= iscale + 1:
        pa[2] = (pa[2]**1)*maxrhoslope/2.
    else:
        pa[2] = (pa[2]**1)*maxrhoslope

    rdef = gp.xepol
    for i in range(3, gp.nrho-1):
        # all -dlog(rho)/dlog(r) at data points and 2,4,8rmax can
        # lie in between 0 and maxrhoslope
        scale = nrtol
        if monotonic:
            # only increase n(r), use pa[i]>=0 directly
            pa[i] = pa[i-1]+pa[i]        * scale * (np.log(rdef[i-2])-np.log(rdef[i-3]))
        else:
            # use pa => [-1, 1] for full interval
            pa[i] = pa[i-1]+(pa[i]-0.5)*2. * scale * (np.log(rdef[i-2])-np.log(rdef[i-3]))

        pa[i] = max(0., pa[i])
        if i <= iscale+1: # iscale: no. bins with xipol<Xscale[all]
            pa[i] = min(maxrhoslope/2., pa[i])
        else:
            pa[i] = min(maxrhoslope, pa[i])
    # rho slope for asymptotically reaching r = \infty is given directly
    # must lie below -1
    # to ensure we have a finite mass at all radii 0<r<=\infty
    pa[gp.nrho-1] = pa[gp.nrho-1]*maxrhoslope# *(maxrhoslope-1)+1
    if monotonic:
        pa[gp.nrho-1] += pa[gp.nrho-2]

    # finite mass prior: to bound between 3 and maxrhoslope, favoring 3:
    # pa[gp.nrho-1] = max(pa[gp.nrho-1], 3.)
    return pa
## \fn map_nr(pa, prof, pop, gp)
# mapping nr and rho parameters from [0,1] to full parameter space
# setting all n(r<r_{iscale})<=2.0
# and possibly a monotonically increasing function
# first parameter is offset for rho_half
# second parameter is asymptotic n(r to 0) value
# @param pa [0,1]^Nipol
# @param prof string for which profile to map, nu or rho
# @param pop integer for population number, 0 for rho*, 1,2 for tracers
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

        # rho*
        offstep = gp.nrho
        tmp_rhostar = map_nr(pc[off:off+offstep], 'nu', 0, gp)
        for i in range(offstep):
            pc[off+i] = tmp_rhostar[i]
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
