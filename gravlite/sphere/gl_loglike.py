#!/usr/bin/env python3

## @file
# define log likelihood function to be called by pyMultinest
# spherical version

import numpy as np
import pdb
from gl_physics import rho, beta, sig_kap_los
from gl_class_profiles import Profiles
from gl_priors import check_rho, check_nr, check_bprior, check_beta
from gl_chi import calc_chi2
import gl_helper as gh
    
def geom_loglike(cube, ndim, nparams, gp):
    tmp_profs = Profiles(gp.pops, gp.nipol)
    off = 0

    rho_param = np.array(cube[off:off+gp.nepol])
    # if check_nr(rho_param[2:-1]):
    #     print('dn/dr too big!')
    #     return gh.err(0.7, gp)

    tmp_rho = rho(gp.xepol, rho_param, gp)
    if(check_rho(tmp_rho, gp.rprior, gp.rhotol)):
        print('rho slope error') # should be ensured by [0,1] priors and maxslope
        return gh.err(1.)
    tmp_profs.set_rho(tmp_rho[:gp.nipol])
    off += gp.nepol

    nuparstore = []
    for pop in np.arange(gp.pops)+1:
        nu_param = cube[off:off+gp.nepol]
        nuparstore.append(nu_param)
        
        tmp_nu = rho(gp.xepol, nu_param, gp) #  [1], [pc]
        # if check_rho(tmp_nu, gp):
        #     print('nu error')
        #     return err/2.
        # if gp.bprior and check_bprior(tmp_rho, tmp_nu):
        #      print('bprior error')
        #      return gh.err(1.5, gp)
        tmp_profs.set_nu(pop, tmp_nu[:gp.nipol]) # [munit/pc^3]
        off += gp.nepol

        beta_param = np.array(cube[off:off+gp.nbeta])
        tmp_beta = beta(gp.xipol, beta_param, gp)
        if check_beta(tmp_beta, gp):
            print('beta error')
            return gh.err(2., gp)
        tmp_profs.set_beta(pop, tmp_beta)
        off += gp.nbeta

        try:
            sig, kap = sig_kap_los(gp.xepol, pop, rho_param, nu_param, beta_param, gp)
            # sig and kap already are on data radii only, so no extension by 3 bins here
        except Exception as detail:
            return gh.err(3., gp)
        tmp_profs.set_sig_kap(pop, sig, kap)
    
    # determine log likelihood (*not* reduced chi2)
    chi2 = calc_chi2(tmp_profs, nuparstore, gp)
    print('found log likelihood = ', -chi2/2.)
    return -chi2/2.   # from   likelihood L = exp(-\chi^2/2), want log of that
