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
import gl_physics as phys
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

def map_kr(params, prof, pop, gp):   # Not used # not supported
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


    z_points_tmp = np.append(0., gp.z_bincenter_vecs[pop])
    for jter in range(0, gp.nbins[pop]):
        z_diff = z_points_tmp[jter+1] - z_points_tmp[jter]
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
    D_range = gp.simplenu_baryon_D_max - gp.simplenu_baryon_D_min

    if gp.prior_type_simplenu_baryon == 'linear':
        K = gp.simplenu_baryon_K_min + K_range*params[0]
        D = gp.simplenu_baryon_D_min + D_range*params[1]
    elif gp.prior_type_simplenu_baryon == 'gaussian':
        K_mid = gp.simplenu_baryon_K_mid
        K_sd  = gp.simplenu_baryon_K_sd
        D_mid = gp.simplenu_baryon_D_mid
        D_sd = gp.simplenu_baryon_D_sd
        K = K_mid + np.sqrt(2) * K_sd * scipy.special.erfinv(2*(params[0]-0.5))
        D = D_mid + np.sqrt(2) * D_sd * scipy.special.erfinv(2*(params[1]-0.5))

    return np.array([K, D])

def map_trivial_baryon(params, gp):
    Sig_bary_min = gp.trivial_baryon_Sig_min
    #print ('before')
    #pdb.set_trace()
    Sig_bary_range = gp.trivial_baryon_Sig_max - gp.trivial_baryon_Sig_min
    Sig_bary = Sig_bary_min + Sig_bary_range*params[0]
    return np.array([Sig_bary])

def map_obs_baryon(params, gp):
    Sig_tot = gp.obs_baryon_Sig_tot*(1.+ gp.obs_baryon_Sig_tot_err*(2.*params[0]-1.))
    Sig_dwarf = gp.obs_baryon_Sig_dwarf*(1.+ gp.obs_baryon_Sig_dwarf_err*(2.*params[1]-1.))

    rho0_MS_thick = gp.obs_baryon_rho0_MS * gp.obs_baryon_MS_beta_max*params[2]

    h_bary  = gp.obs_baryon_h *(1. + gp.obs_baryon_h_err*(2.*params[3]-1.))
    h1_bary = gp.obs_baryon_h1*(1. + gp.obs_baryon_h1_err*(2.*params[4]-1.))
    h2_bary = gp.obs_baryon_h2*(1.+ gp.obs_baryon_h2_err*(2.*params[5]-1.))
    h3_bary = gp.obs_baryon_h3*(1.+ gp.obs_baryon_h3_err*(2.*params[6]-1.))
    x_bary  = gp.obs_baryon_x* (1.+ gp.obs_baryon_x_err*(2.*params[7]-1.))

    Sig_HII = gp.obs_baryon_Sig_HII*(1.+ gp.obs_baryon_Sig_HII_err*(2.*params[8]-1.))
    h_HII = gp.obs_baryon_h_HII*(1.+ gp.obs_baryon_h_HII_err*(2.*params[9]-1.))


    if gp.wide_bary_range:  # rescale to allow for smaller bary dens
        Sig_dwarf = Sig_dwarf*Sig_tot/gp.obs_baryon_Sig_tot
        Sig_HII = Sig_HII*Sig_tot/gp.obs_baryon_Sig_tot
        rho0_MS_thick = rho0_MS_thick*Sig_tot/gp.obs_baryon_Sig_tot


    return np.array([Sig_tot, Sig_dwarf, rho0_MS_thick, h_bary, h1_bary, h2_bary, h3_bary, x_bary, Sig_HII, h_HII])


def map_simplenu_baryon_gaussian(params, gp): # not used  # not supported
    rho_mid_vector = gp.gaussian_rho_baryon_mid_vector
    rho_SD_vector = gp.gaussian_rho_baryon_SD_vector
    rho_baryon_out = []

    for jter in range(0,sum(gp.nbins)+1):
        rho_baryon_jter = rho_mid_vector[jter] + np.sqrt(2) * rho_SD_vector[jter] * scipy.special.erfinv(2*(params[jter]-0.5))
        #rho_baryon_jter = rho_mid_vector[jter] + rho_SD_vector[jter]*(2*params[jter] - 1) #flat prior test
        rho_baryon_out.append(rho_baryon_jter)

    return np.array(rho_baryon_out)


def map_gaussian_per_bin_dm(params, gp): # not used  # not supported
    med = gp.rho_DM_gaussian_med = 10.0E6
    sd = gp.rho_DM_gaussian_SD = 5.0E6

    rho_DM_out = []

    for jter in range(0,sum(gp.nbins)+1):
        rho_DM_jter = med + np.sqrt(2) * sd * scipy.special.erfinv(2*(params[jter]-0.5))
        rho_DM_out.append(rho_DM_jter)

    return np.array(rho_DM_out)



def map_hypererr(param, prof, pop, gp):
    if prof == 'nu':
        meanerr = gp.dat.meannuerr
    elif prof == 'sig':
        meanerr = gp.dat.meansigz2err
    minhyper = 0.1*meanerr   # Hard coded, put in gl_params
    maxhyper = 10.*meanerr   # Hard coded, put in gl_params
    lmax = 1./minhyper
    lmin = 1./maxhyper
    lam = 1/(param[0]*(lmax-lmin)+lmin)  # lam : maxhyper -> minhyper
    return lam
## \fn map_hypererr
# return hyperparameter

def map_sigRz2_model(params,tracer_pop,gp):
    A_range = gp.tilt_A_max - gp.tilt_A_min
    A = gp.tilt_A_min + A_range*params[0]

    #print('min A = ', gp.tilt_A_min, ', A = ', A, ', max A = ', gp.tilt_A_max)
    n_range = gp.tilt_n_max - gp.tilt_n_min
    n = gp.tilt_n_min + n_range*params[1]
    if n < 0:
        print('Negative n found')
        print('gp.tilt_n_max = ', gp.tilt_n_max)
        print('gp.tilt_n_min = ', gp.tilt_n_min)
        print('n = ', n)


    if gp.prior_type_tilt_R == 'linear':
        #R_range = gp.tilt_R_max - gp.tilt_R_min
        #R = gp.tilt_R_min + R_range*params[2]
        k_range = gp.tilt_k_max[tracer_pop] - gp.tilt_k_min[tracer_pop]
        k = gp.tilt_k_min[tracer_pop] + k_range*params[2]
    elif gp.prior_type_tilt_R == 'gaussian':
        R = gp.tilt_R_med + np.sqrt(2) * gp.tilt_R_sd * scipy.special.erfinv(2*(params[2]-0.5))


    #return np.array([A,n,R])
    return np.array([A,n,k])
# Input: 3 multinest cube params,
# Output: tilt parameters A, n and R = 2*R_0*R_1/(R_0 + R_1)
# SS: 19 May 2015

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

def map_rhonu(params, prof, pop, gp):   # Not used
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

def map_nu_data(params, pop, gp): # not used  # not supported
    nudat    = gp.dat.nu[pop]
    nuerr    = gp.dat.nuerr[pop]
    nu_out_vector = []

    nu_C = nudat[0] + np.sqrt(2)*5*nuerr[0]* scipy.special.erfinv(2*(params[0]-0.5))

    for jter in range(0, gp.nbins[pop]):
        nu_i = nudat[jter] + np.sqrt(2) * nuerr[jter] * scipy.special.erfinv(2*(params[jter+1]-0.5))
        nu_out_vector.append(nu_i)

    return np.hstack([nu_C, nu_out_vector])




def map_constdm(params, prof, pop, gp):
    if gp.rho_C_prior_type == 'linear':   # Used
        rho_C = gp.rho_C_min + (gp.rho_C_max-gp.rho_C_min)*params[0]
    elif gp.rho_C_prior_type == 'log':
        rho_C = np.log10(gp.rho_C_min) + (np.log10(gp.rho_C_max)-np.log10(gp.rho_C_min))*params[0]
        rho_C = 10**rho_C
    return [rho_C]


def map_simplenu_dm(params, prof, pop, gp):
    if gp.rho_C_prior_type == 'linear':
        rho_C = gp.rho_C_min + (gp.rho_C_max-gp.rho_C_min)*params[0]
    elif gp.rho_C_prior_type == 'log':
        rho_C = np.log10(gp.rho_C_min) + (np.log10(gp.rho_C_max)-np.log10(gp.rho_C_min))*params[0]
        rho_C = 10**rho_C
    #K_range = gp.simplenu_dm_K_max - gp.simplenu_dm_K_min
    #K = gp.simplenu_dm_K_min + K_range*params[1]
    Sig_DD_inf_range = gp.simplenu_DD_Sig_inf_max - gp.simplenu_DD_Sig_inf_min
    Sig_DD_inf = gp.simplenu_DD_Sig_inf_min + Sig_DD_inf_range*params[1]
    D_range = gp.simplenu_dm_D_max - gp.simplenu_dm_D_min
    D = gp.simplenu_dm_D_min + D_range*params[2]
    return np.array([rho_C, Sig_DD_inf, D])


def map_nu_exponential_sum(params, pop, gp):
    #[nuC_median, nuC_SD, nu_h, nu_h_SD]
    #exp_ii = which exponential in the sum
    output = np.array([])
    nu_C_median  = gp.nu_exp_sum_priors[pop][0][0]
    nu_C_SD      = gp.nu_exp_sum_priors[pop][0][1]
    nu_h_median = gp.nu_exp_sum_priors[pop][0][2]
    nu_h_SD     = gp.nu_exp_sum_priors[pop][0][3]
    
    # SS: Picking x from flat distr for nu=10**x:
    nu_C_exp = nu_C_median-nu_C_SD + 2.*nu_C_SD*params[0] 
    nu_C = 10**nu_C_exp
    #nuC = nuC_median + np.sqrt(2) * nuC_SD * scipy.special.erfinv(2*(params[exp_ii*2]-0.5)) #Gaussian
    nu_h_min = nu_h_median-nu_h_SD
    nu_h = nu_h_min + 2*nu_h_SD*params[1]
   #nu_h = nu_h_median + np.sqrt(2) * nu_h_SD * scipy.special.erfinv(2*(params[exp_ii*2 + 1] - 0.5))


    output = np.append(output, (nu_C, nu_h))

    if gp.N_nu_model_exps == 2:  # nu = A*exp(-x/h) - B*exp(-x/k) (>0 så A>B och h>k)
        nu_C_median  = gp.nu_exp_sum_priors[pop][1][0]
        nu_C_SD      = gp.nu_exp_sum_priors[pop][1][1]
        nu_h_median = gp.nu_exp_sum_priors[pop][1][2]
        nu_h_SD     = gp.nu_exp_sum_priors[pop][1][3]

        C_exp_low = nu_C_median - nu_C_SD
        C_exp_high = min(nu_C_median + nu_C_SD , nu_C_exp) 
        h_low = nu_h_median - nu_h_SD
        h_high = min(nu_h_median + nu_h_SD , nu_h)

        nu_C_exp = C_exp_low + (C_exp_high-C_exp_low)*params[2]
        nu_C = 10**nu_C_exp
        nu_h = h_low + (h_high-h_low)*params[3]
        output = np.append(output, (-1.*nu_C, nu_h)) # The second term for nu-fit is negative

        #nuC = params[exp_ii*2] * (nuC_median + nuC_SD) # Linear RE DO THIS
        #nuC = nuC_median-nuC_SD + 2.*nuC_SD*params[exp_ii*2]  # SS
        #print ('gl_class_cube: nu_C:',nu_C,' nu_h:',nu_h)
        #if output[2] > 0.5*output[0]: print ('output:',output)

    return output


class Cube:
    def __init__ (self, gp):
        self.pops = gp.ntracer_pops
        # for density and (nu, tilt)_i
        self.cube = np.zeros(gp.ndim)
        return
    ## \fn __init__ (self, pops)
    # constructor, with modes depending on locpop


    def convert_to_parameter_space(self, gp):
        #print ('In convert_to_parameter_space  !!!!!!!!!!!!!!!')
        # if we want any priors, here they have to enter:

        pc = self.cube
        off = 0

        #Dark Matter mass profile parameters: rho_C, kz_C, kz_vector
        if gp.darkmattermodel == 'const_dm':
            offstep = 1
            tmp_DM = map_constdm(pc[off:off+offstep], 'rho', 0, gp)
        elif gp.darkmattermodel == 'ConstPlusDD':
            offstep = 3
            tmp_DM = map_simplenu_dm(pc[off:off+offstep], 'rho', 0, gp)
        elif gp.darkmattermodel == 'gaussian_per_bin': # not used
            offstep = gp.nrhonu
            tmp_DM = map_gaussian_per_bin_dm(params, gp)
        elif gp.darkmattermodel == 'kz_dm': # not used
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
            elif gp.baryonmodel == 'trivial_baryon':
                tmp_baryon = map_trivial_baryon(pc[off:off+offstep], gp)
            elif gp.baryonmodel == 'obs_baryon':
                tmp_baryon = map_obs_baryon(pc[off:off+offstep], gp)
            elif gp.baryonmodel == 'kz_baryon':
                tmp_baryon = map_kr(pc[off:off+offstep], 'rho', baryon_pop, gp)
            elif gp.baryonmodel == 'simplenu_baryon_gaussian':
                tmp_baryon = map_simplenu_baryon_gaussian(pc[off:off+offstep], gp)

            for i in range(offstep):
                pc[off+i] = tmp_baryon[i]
            off += offstep

        #Tracer profile parameters: nu_C, kz_nu_C, kz_nu_vector # kz_nu_LS
        tmp_tracer_params = [None] * gp.ntracer_pops
        for tracer_pop in range(0, gp.ntracer_pops):
            if gp.scan_rhonu_space:
                tmp_tracer = map_rhonu(pc[off:off+offstep], 'nu', tracer_pop, gp)
            elif gp.nu_model=='kz_nu': # not supported
                offstep = gp.nbins[tracer_pop] + 1 + 1 #kz on bincenters, and zC=0, and nu_C
                tmp_tracer = map_kr(pc[off:off+offstep], 'nu', tracer_pop, gp)
            elif gp.nu_model=='gaussian_data': # not supported
                offstep = gp.nbins[tracer_pop] + 1 #nu on bincenters, and nu_C
                tmp_tracer = map_nu_data(pc[off:off+offstep], tracer_pop, gp)
            elif gp.nu_model=='exponential_sum':
                offstep = gp.N_nu_model_exps*2
                tmp_tracer = map_nu_exponential_sum(pc[off:off+offstep],tracer_pop, gp)
                tmp_tracer_params[tracer_pop] = tmp_tracer
                
            for i in range(offstep):
                pc[off+i] = tmp_tracer[i]
            off += offstep

        #print('pc = ', pc[0:gp.ndim])

        # Introducing tilt term:
        if gp.tilt:
            for tracer_pop in range(0, gp.ntracer_pops):
                offstep = gp.ntilt_params
                tmp_tilt = map_sigRz2_model(pc[off:off+offstep],tracer_pop,gp)
                for i in range(offstep):
                    pc[off+i] = tmp_tilt[i]
                off += offstep

        # Possibility to use hyperparameters:
        if gp.hyperparams:
            offstep = 1
            tmp_hypernu = map_hypererr(pc[off:off+offstep], 'nu', 0, gp)
            pc[off] = tmp_hypernu
            off += offstep

            offstep = 1
            tmp_hypersig = map_hypererr(pc[off:off+offstep], 'sig', 0, gp)
            pc[off] = tmp_hypersig
            off += offstep

        if gp.analytic_sigz2 == False:
        #Normalisation constant C, for the calculation of sigma_z via Jeans Eq.
        #one for each tracer population
            offstep = gp.ntracer_pops
            for t_pop in range(0, gp.ntracer_pops):
                if gp.nu_model == 'exponential_sum':
                    tmp_nu0 = 0.
                    for jter in range(0,2*gp.N_nu_model_exps,2):
                    #tmp_nu0 += tmp_tracer_params[t_pop][jter]*np.exp(-1.*gp.z_bincenter_vecs[t_pop][0]/tmp_tracer_params[t_pop][jter+1]) 
                    #z_temp = gp.z_bincenter_vecs[t_pop][0] - gp.nu_z_bincenter_vecs[t_pop][0] # Redefined z=0 to be at first nu bin
                        z_temp = gp.z_bincenter_vecs[t_pop][-1] - gp.nu_z_bincenter_vecs[t_pop][0]
                        # Above: C defined at gp.z_bincenter_vecs[pop][-1]
                        # For nu: z=0 defined at gp.nu_z_bincenter_vecs[pop][0]
                        tmp_nu0 += tmp_tracer_params[t_pop][jter]*np.exp(-z_temp/tmp_tracer_params[t_pop][jter+1]) 
                #print ('tmp_nu0:',tmp_nu0,' z0:',gp.z_bincenter_vecs[t_pop][0],' pop:',t_pop)
                    IntC_max=(gp.sigz_C_max[t_pop]**2)*tmp_nu0  # SS
                    IntC_min=(gp.sigz_C_min[t_pop]**2)*tmp_nu0  # SS

                #IntC_max=(gp.sigz_C_max**2)*tmp_tracer_params[t_pop][0]  # SS
                #IntC_min=(gp.sigz_C_min**2)*tmp_tracer_params[t_pop][0]  # SS
                    if gp.print_flag:
                        print ('sigz_C_max:',gp.sigz_C_max,'****************')
                        print ('sigz_C_min:',gp.sigz_C_min)
                        print ('tmp_nu0:',tmp_nu0,' A:',tmp_tracer_params[t_pop][0])
                        print ('pop:',t_pop,' IntC_max:',IntC_max,' IntC_min:',IntC_min)
                else:  # False
                    IntC_max=(gp.sigz_C_max**2)*gp.nu_C_max
                    IntC_min=(gp.sigz_C_min**2)*gp.nu_C_min

                IntC = IntC_min+(IntC_max-IntC_min)*pc[off+t_pop]
                pc[off+t_pop] = IntC  # C at z_min now (not at 0)

            off += offstep

        if off != gp.ndim:
            print ('in gl_class_cube, off:',off,'gp.ndim:',gp.ndim)
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
