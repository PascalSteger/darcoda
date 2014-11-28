#!/usr/bin/env python3

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

# (c) GPL v3 2014 ETHZ Pascal S.P. Steger

import numpy as np
import pdb
import gl_helper as gh

def dierfc(y):
    y = np.array(y)
    qa=9.16461398268964e-01
    qb=2.31729200323405e-01
    qc=4.88826640273108e-01
    qd=1.24610454613712e-01
    q0=4.99999303439796e-01
    q1=1.16065025341614e-01
    q2=1.50689047360223e-01
    q3=2.69999308670029e-01
    q4=-7.28846765585675e-02
    pa=3.97886080735226000e+00
    pb=1.20782237635245222e-01
    p0=2.44044510593190935e-01
    p1=4.34397492331430115e-01
    p2=6.86265948274097816e-01
    p3=9.56464974744799006e-01
    p4=1.16374581931560831e+00
    p5=1.21448730779995237e+00
    p6=1.05375024970847138e+00
    p7=7.13657635868730364e-01
    p8=3.16847638520135944e-01
    p9=1.47297938331485121e-02
    p10=-1.05872177941595488e-01
    p11=-7.43424357241784861e-02
    p12=2.20995927012179067e-03
    p13=3.46494207789099922e-02
    p14=1.42961988697898018e-02
    p15=-1.18598117047771104e-02
    p16=-1.12749169332504870e-02
    p17=3.39721910367775861e-03
    p18=6.85649426074558612e-03
    p19=-7.71708358954120939e-04
    p20=-3.51287146129100025e-03
    p21=1.05739299623423047e-04
    p22=1.12648096188977922e-03

    # remove insensible ranges
    y[y==0] = 1e-15
    y[y==1] = 1.-1e-15

    z=1.*y
    w=qa-np.log(z)
    u=np.sqrt(w)
    s=(qc+np.log(u))/w
    t=1/(u+qb)
    x=u*(1-s*(0.5+s*qd))-((((q4*t+q3)*t+q2)*t+q1)*t+q0)*t
    t=pa/(pa+x)
    u=t-0.5
    s=(((((((((p22*u+p21)*u+p20)*u+p19)*u+p18)*u+p17)*u+p16)*u+p15)*u+p14)*u+p13)*u+p12
    s=((((((((((((s*u+p11)*u+p10)*u+p9)*u+p8)*u+p7)*u+p6)*u+p5)*u+p4)*u+p3)*u+p2)*u+p1)*u+p0)*t-z*np.exp(x*x-pb)
    x=x+s*(1+x*s)
    return x
## \fn dierfc(y)
# inverse of complimentary error function from MultiNest implementation
# @param y

def ginv(x, mu, sigma):
    # parameter(SqrtTwo=1.414213562d0)
    return mu+sigma*np.sqrt(2)*dierfc(2.*(1.-x))
## \fn ginv(x, mu, sigma)
# take uniform range [0,1] and return normal pdf centered around mu, with width sigma
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
        width = gp.log10rhospread
        innerslope = gp.innerslope
        nrscale = gp.nrtol/(max(np.log(gp.xipol))-min(np.log(gp.xipol)))
        monotonic = gp.monotonic
    elif prof=='nu':
        rhoscale = gp.dat.nuhalf[pop]
        width = gp.log10nuspread
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
    dnrdlrparams = params[1:]

    for k in range(0, gp.nepol):
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

    # rho slope for asymptotically reaching r = \infty is given directly
    # must lie below -3, thus n(r)>3
    deltalogrlast = (np.log(gp.xepol[-1])-np.log(gp.xepol[-2]))

    # finite mass prior: to bound between 3 and gp.maxrhoslope, favoring 3:
    nrasyminfty = max(nr[-1], 3.001)

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
