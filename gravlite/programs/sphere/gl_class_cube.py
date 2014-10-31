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

def dierfc(y):
    z = y
    if (y > 1):
        z = 2 - y
    w = 0.916461398268964 - np.log(z)
    u = np.sqrt(w)
    s = (np.log(u) + 0.488826640273108) / w
    t = 1 / (u + 0.231729200323405)
    x = u * (1 - s * (s * 0.124610454613712 + 0.5)) -\
        ((((-0.0728846765585675 * t + 0.269999308670029) * t +\
        0.150689047360223) * t + 0.116065025341614) * t +\
        0.499999303439796) * t
    t = 3.97886080735226 / (x + 3.97886080735226)
    u = t - 0.5
    s = (((((((((0.00112648096188977922 * u +\
        1.05739299623423047e-4) * u - 0.00351287146129100025) * u -\
        7.71708358954120939e-4) * u + 0.00685649426074558612) * u +\
        0.00339721910367775861) * u - 0.011274916933250487) * u -\
        0.0118598117047771104) * u + 0.0142961988697898018) * u +\
        0.0346494207789099922) * u + 0.00220995927012179067
    s = ((((((((((((s * u - 0.0743424357241784861) * u -\
        0.105872177941595488) * u + 0.0147297938331485121) * u +\
        0.316847638520135944) * u + 0.713657635868730364) * u +\
        1.05375024970847138) * u + 1.21448730779995237) * u +\
        1.16374581931560831) * u + 0.956464974744799006) * u +\
        0.686265948274097816) * u + 0.434397492331430115) * u +
        0.244044510593190935) * t -\
        z * np.exp(x * x - 0.120782237635245222)
    x += s * (x * s + 1)
    if (y > 1):
        x = -x
    return x
## \fn dierfc(y)
# inverse error function
# @param y

def ginv(x, mu, sigma):
    # TODO Uniform[0:1]  ->  Gaussian[mean=mu,variance=sigma**2]
    # double precision r,mu,sigma,GaussianPrior
    # double precision SqrtTwo
    # parameter(SqrtTwo=1.414213562d0)
    if x <= 1.0e-16 or (1.-x) <= 0:
        GaussianPrior = -1.0e32
    else:
        GaussianPrior = mu+sigma*np.sqrt(2)*dierfc(2.*(1.-x))
    return GaussianPrior
## \fn ginv(x, mu, sigma)
# take range [0,1] and return [-inf, inf]
# sampled with a normal distribution
# @param x random variable sampled by MultiNest, uniformly in range [0,1]
# @param mu central value
# @param sigma width of normal distribution

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
