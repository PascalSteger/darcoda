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

import pdb
import numpy as np
import gl_helper as gh
import scipy.special

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
    #pdb.set_trace()
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

def map_kr(params, prof, pop, gp):
    # params[0] for central value of rho or nu
    # params[1] for central value of kz (for rho or nu)
    #pdb.set_trace()
    if prof == 'rho':
        rhonu_C_max = gp.rho_C_max
        rhonu_C_min = gp.rho_C_min
        rhonu_C_prior_type = gp.rho_C_prior_type
        kz_C_max = gp.kz_rho_C_max
        kz_C_min = gp.kz_rho_C_min
        monotonic = gp.monotonic_rho
        kz_rhonu_selection = gp.kz_rho_selection
    elif prof == 'nu':
        rhonu_C_max = gp.nu_C_max
        rhonu_C_min = gp.nu_C_min
        rhonu_C_prior_type = gp.nu_C_prior_type
        kz_C_max = gp.kz_nu_C_max
        kz_C_min = gp.kz_nu_C_min
        monotonic = gp.monotonic_nu
        kz_rhonu_selection = gp.kz_nu_selection

    # rho or nu central value
    if rhonu_C_prior_type == 'linear':
        rhonu_C = rhonu_C_min + (rhonu_C_max-rhonu_C_min)*params[0]
    elif rhonu_C_prior_type == 'log':
        rhonu_C = np.log10(rhonu_C_min) + (np.log10(rhonu_C_max)-np.log10(rhonu_C_min))*params[0]
        rhonu_C = 10**rhonu_C

    # kz for rho or nu, central value
    kz_C = kz_C_min + (kz_C_max - kz_C_min)*params[1]

    # Starting from the central value walk the kz value, limited
    # by max_kz_slope (=dk/dz)
    kz_vector = []
    kz_i_m1 = kz_C #kz_(i-1)

    for jter in range(0, gp.nbins):
        z_diff = gp.z_all_pts[jter+1] - gp.z_all_pts[jter]
        if kz_rhonu_selection == 'tophat':
            kz_max = kz_i_m1 + gp.max_kz_slope*z_diff
            kz_min = kz_i_m1 - gp.max_kz_slope*z_diff
            kz_i = kz_min + (kz_max - kz_min)*params[2+jter]

        elif kz_rhonu_selection == 'gaussian':
            sigma_kz = gp.max_kz_slope*z_diff/2. # TODO: Think on this some more
            #if monotonic:
            #    kz_i_m1 = max(kz_i_m1, 0)
            kz_i = kz_i_m1 + np.sqrt(2) * sigma_kz * scipy.special.erfinv(2*(params[jter]-0.5))

        #if monotonic:
        #    kz_i = max(kz_i, 0.) #Sledgehammer method for monotonic implementation
            #if monotonic:
            #    kz_min=max(kz_min, 0.)

        kz_i_m1 = kz_i
        kz_vector.append(kz_i)

    ##kz Last Star
    #kz_LS_max = kz_i + gp.max_kz_slope*(gp.z_all_pts[-1] - gp.z_all_pts[-2])
    #kz_LS_max = kz_i - gp.max_kz_slope*(gp.z_all_pts[-1] - gp.z_all_pts[-2])
    #kz_LS = kz_LS_max*params[3+jter] + kz_min*(1-params[3+jter])

    return np.hstack([rhonu_C, kz_C, kz_vector])


def map_simplenu_baryon(params, gp):
    # Input: two multinest cube params, [0,1]
    # Output: K and D parameters for simplenu baryon model
    # H Silverwood 3/02/15

    K_range = gp.simplenu_baryon_K_max - gp.simplenu_baryon_K_min
    K = gp.simplenu_baryon_K_min + K_range*params[0]
    D_range = gp.simplenu_baryon_D_max - gp.simplenu_baryon_D_min
    D = gp.simplenu_baryon_D_min + D_range*params[1]

    return np.array([K, D])


def map_nr(params, prof, pop, gp):
    gh.sanitize_vector(params, gp.nrho, 0, 1, gp.debug)
    nr = np.zeros(gp.nepol) # to hold the n(r) = dlog(rho)/dlog(r) values

    # get offset and n(r) profiles, calculate rho
    if prof=='rho':
        rhoscale = gp.rhohalf
        Rscale = gp.Xscale[0]
        width = gp.log10rhospread
        nrscale = gp.nztol/(max(np.log(gp.xipol))-min(np.log(gp.xipol)))
        monotonic = gp.monotonic
    elif prof=='nu':
        rhoscale = gp.dat.nuhalf[pop]
        Rscale = gp.Xscale[pop]
        width = gp.log10nuspread
        nrscale = gp.nztol_nu/(max(np.log(gp.xipol))-min(np.log(gp.xipol)))
        monotonic = gp.monotonic_nu
    else:
        raise Exception('wrong profile in gl_class_cube.map_nr')

    # zeroth parameter gives half-light radius value of rho directly
    # use [0,1]**3 to increase probability of sampling close to 0
    # fix value with tracer densities,
    # sample a flat distribution over log(rho_half)
    rhohalf = 10**((params[0]-0.5)*2.*width+np.log10(rhoscale))

    # nr(r=0) is = rho slope for approaching r=0 asymptotically, given directly
    # should be smaller than -1 to exclude infinite enclosed mass
    nrasym0 = params[1]*0.99

    # work directly with the dn(r)/dlog(r) parameters here
    dnrdlrparams = params[2:]
    for k in range(0, gp.nepol):
        deltalogr = (np.log(gp.xepol[k-1])-np.log(gp.xepol[k-2]))
        # construct n(r_k+1) from n(r_k)+dn/dlogr*Delta log r, integrated
        if monotonic:
            # only increase n(r), use pa[i]>=0 directly
            nr[k] = nr[k-1] + dnrdlrparams[k-1] * nrscale * deltalogr
        else:
            # use pa => [-1, 1] for full interval
            nr[k] = nr[k-1] + (dnrdlrparams[k-1]-0.5)*2. * nrscale * deltalogr
        # cut at zero: we do not want to have density rising outwards
        nr[k] = max(0., nr[k])
    # rho slope for asymptotically reaching r = \infty is given directly
    # must lie below -1, thus n(r)>1, to ensure finite total mass
    deltalogrlast = (np.log(gp.xepol[-1])-np.log(gp.xepol[-2]))
    if monotonic:
        nrasyminfty = nr[-1]+dnrdlrparams[-1] * nrscale * deltalogrlast
    else:
        nrasyminfty = nr[-1]+(dnrdlrparams[-1]-0.5)*2 * nrscale * deltalogrlast

    # finite mass prior: nrasyminfty must be > 1:
    if nrasyminfty < 1.001:
        if nr[-1]+1*nrscale*deltalogrlast > 1.001:
            nrasyminfty = (nr[-1]+1*nrscale*deltalogrlast-1.001)*dnrdlrrarams[-1]+1.001
        else:
            nrasyminfty = -9999
            raise Exception('too low asymptotic n(r)')

    params = np.hstack([rhohalf, nrasym0, nr, nrasyminfty])
    return params
## \fn map_nr(params, prof, pop, gp)
# mapping rho parameters from [0,1] to full parameter space
# prior on dn/dlnr: < nrscale
# and possibly a monotonically increasing function
# zeroth parameter is offset for rho_half
# first parameter is asymptotic n(r to 0) value
# @param params cube [0,1]^ndim
# @param prof string nu, rho, rhostar
# @param pop population int, 0 for rho*, 1,2,... for tracer densities
# @param gp global parameters


#def map_nu(pa, gp):
#    # TODO: assertion len(pa)=gp.nepol
#    for i in range(len(pa)):
#        pa[i] = 10**(pa[i]*(gp.maxlog10nu-gp.minlog10nu)+gp.minlog10nu)
#    return pa
### \fn map_nu(pa, gp)
## map tracer densities, directly
## if used, put
##         self.maxlog10nu = 4.     # direct sampling of nu: min value
##        self.minlog10nu = 0.     # direct sampling of nu: max value
## in gl_params
## @param pa cube [0,1]^n
## @param gp global parameters


def map_MtoL(pa, gp):
    scale = gp.MtoLmax - gp.MtoLmin
    pa = pa*scale+gp.MtoLmin
    return pa
## \fn map_MtoL(pa, gp)
# map [0,1] to MtoL flat prior
# @param pa scalar
# @param gp global parameters holding MtoL{min,max}

def map_rhonu(params, prof, pop, gp):
    if prof == 'rho':
        rhonu_C_max = gp.rho_C_max
        rhonu_C_min = gp.rho_C_min
        monotonic = gp.monotonic_rho
        prior_type = gp.prior_type_rho
    elif prof == 'nu':
        rhonu_C_max = gp.nu_C_max
        rhonu_C_min = gp.nu_C_min
        monotonic = gp.monotonic_nu
        prior_type = gp.prior_type_nu

    rhonu_C = rhonu_C_min + (rhonu_C_max-rhonu_C_min)*params[0]

    # Pick rho/nu values between rhonu_C_max and rhonu_C_min, or if monotonic
    # pick between rhonu_i_m1 and rhonu_C_min
    rhonu_vector = []
    rhonu_i_m1 = rhonu_C
    for jter in range(1, gp.nbins+1):
        if monotonic:
            rhonu_C_max = rhonu_i_m1

        if prior_type == 'linear':
            rhonu_value = rhonu_C_min + (rhonu_C_max-rhonu_C_min)*params[jter]
        elif prior_type == 'log':
            rhonu_value = np.log10(rhonu_C_min) + (np.log10(rhonu_C_max)-np.log10(rhonu_C_min))*params[jter]
            rhonu_value = 10**rhonu_value
        elif prior_type == 'gaussian':
            sigma_rhonu = (rhonu_C_max - rhonu_C_min)/4.
            rhonu_value = rhonu_i_m1 + np.sqrt(2) * sigma_rhonu * scipy.special.erfinv(2*(params[jter]-0.5))
            rhonu_value = max(rhonu_C_min, rhonu_value)
            rhonu_value = min(rhonu_C_max, rhonu_value)

        rhonu_vector.append(rhonu_value)
        rhonu_i_m1 = rhonu_value

    return np.hstack([rhonu_C, rhonu_vector])

def map_constdm(params, prof, pop, gp):
    rhoC = gp.rho_C_min + (gp.rho_C_max-gp.rho_C_min)*params[0]
    return [rhoC]



class Cube:
    def __init__ (self, gp):
        self.pops = gp.ntracer_pops
        # for density and (nu, tilt)_i
        self.cube = np.zeros(gp.ndim)
        return
    ## \fn __init__ (self, pops)
    # constructor, with modes depending on locpop


    def convert_to_parameter_space(self, gp):
        # if we want any priors, here they have to enter:

        #Normalisation constant C, for the calculation of sigma_z via Jeans Eq.
        pc = self.cube
        off = 0
        offstep = 1
        IntC_max=(gp.sigz_C_max**2)*gp.nu_C_max
        IntC_min=(gp.sigz_C_min**2)*gp.nu_C_min
        pc[off] = IntC_min+(IntC_max-IntC_min)*pc[off]
        off += offstep

        #Dark Matter mass profile parameters: rho_C, kz_C, kz_vector
        if gp.darkmattermodel == 'const_dm':
            offstep = 1
            tmp_DM = map_constdm(pc[off:off+offstep], 'rho', 0, gp)
        elif gp.darkmattermodel == 'kz_dm':
            offstep = gp.nrhonu + 1
            if gp.scan_rhonu_space:
                tmp_DM = map_rhonu(pc[off:off+offstep], 'rho', 0, gp)
            else:
                tmp_DM = map_kr(pc[off:off+offstep], 'rho', 0, gp)

        for i in range(offstep):
            pc[off+i] = tmp_DM[i]
        off += offstep

        #Baryon mass profile parameters
        #Redo this when we introduce baryons
        for baryon_pop in range(0, gp.nbaryon_pops):
            offstep = gp.nbaryon_params
            if gp.baryonmodel == 'simplenu_baryon':
                tmp_baryon = map_simplenu_baryon(pc[off:off+offstep], gp)
            elif gp.baryonmodel == 'kz_baryon':
                tmp_baryon = map_kr(pc[off:off+offstep], 'rho', baryon_pop, gp)
            for i in range(offstep):
                pc[off+i] = tmp_baryon[i]
            off += offstep

        #Tracer profile parameters: nu_C, kz_nu_C, kz_nu_vector # kz_nu_LS
        for tracer_pop in range(0, gp.ntracer_pops):
            offstep = gp.nrhonu + 1
            if gp.scan_rhonu_space:
                tmp_tracer = map_rhonu(pc[off:off+offstep], 'nu', tracer_pop, gp)
            else:
                tmp_tracer = map_kr(pc[off:off+offstep], 'nu', tracer_pop, gp)

            for i in range(offstep):
                pc[off+i] = tmp_tracer[i]
            off += offstep

        #print('pc = ', pc[0:gp.ndim])
        #pdb.set_trace()

        if off != gp.ndim:
            gh.LOG(1,'wrong subscripts in gl_class_cube')
            raise Exception('wrong subscripts in gl_class_cube')

        return pc
    ## \fn convert_to_parameter_space(self, gp)
    # convert [0,1]^ndim to parameter space
    # such that values in cube are the parameters we need for rho, nu_i, beta_i
    # @param gp global parameters


    def __repr__(self):
        return "Cube (disc) with "+str(gp.ntracer_pops)+" pops "

    def copy(self, cub):
        self.cube = cub
        return self
    ## \fn copy(self, cub)
    # copy constructor


## \class Cube
# Common base class for all parameter sets
