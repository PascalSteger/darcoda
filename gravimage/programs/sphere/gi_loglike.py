#!/usr/bin/env ipython3

## @file
# define log likelihood function to be called by multinest
# spherical version

# (c) GPL v3 2014 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import pdb
import time
from scipy.interpolate import splev, splrep
from multiprocessing import Pool

# import matplotlib.pyplot as plt
# fig, ax=plt.subplots()
# import matplotlib.animation as animation

import gi_physics as phys
from gi_class_profiles import Profiles
from gi_priors import check_beta
from gi_chi import calc_chi2
import gi_helper as gh
import gi_project as glp


def geom_loglike(cube, ndim, nparams, gp):
    tmp_profs = Profiles(gp.pops, gp.nepol)
    off = 0
    offstep = gp.nrho
    if gp.chi2_Sig_converged <= 0:
        rhodmpar = np.array(cube[off:off+offstep])
        tmp_rho0 = phys.rho(gp.xepol, rhodmpar, 0, gp)
        # rhodmpar hold [rho(rhalf), nr to be used for integration
        # from halflight radius, defined on gp.xepol]
        # (only calculate) M, check
        tmp_M0 = glp.rho_SUM_Mr(gp.xepol, tmp_rho0)
         # store profiles
        tmp_profs.set_prof('nr', 1.*rhodmpar[1+1:-1], 0, gp)
        tmp_profs.set_prof('rho', tmp_rho0, 0, gp)
        tmp_profs.set_prof('M', tmp_M0, 0, gp)

    off += offstep # anyhow, even if Sig not yet converged

    # get profile for rho*
    if gp.investigate == 'obs':
        offstep = gp.nrho
        lbaryonpar = np.array(cube[off:off+offstep])
        rhostar = phys.rho(gp.xepol, lbaryonpar, 0, gp)


        off += offstep
        Signu = glp.rho_param_INT_Sig(gp.xepol, lbaryonpar, 0, gp) # [Munit/pc^2]
        MtoL = cube[off]
        off += 1

        # store these profiles every time
        tmp_profs.set_prof('nu', rhostar, 0, gp)
        tmp_profs.set_prof('Sig', Signu, 0, gp)
        tmp_profs.set_MtoL(MtoL)

    else:
        lbaryonpar = np.zeros(gp.nrho)
        MtoL = 0.

    for pop in np.arange(1, gp.pops+1):  # [1, 2, ..., gp.pops]
        offstep = gp.nrho
        nupar = np.array(cube[off:off+offstep])
        tmp_nrnu = 1.*nupar[1+1:-1]


        tmp_nu = phys.rho(gp.xepol, nupar, pop, gp)
        tmp_Signu = glp.rho_param_INT_Sig(gp.xepol, nupar, pop, gp)
        #tmp_nu = pool.apply_async(phys.rho, [gp.xepol, nupar, pop, gp])
        #tmp_Signu = pool.apply_async(glp.rho_param_INT_Sig, [gp.xepol, nupar, pop, gp])


        off += offstep
        offstep = gp.nbeta
        if gp.chi2_Sig_converged <= 0:
            betapar = np.array(cube[off:off+offstep])
            tmp_beta, tmp_betastar = phys.beta(gp.xepol, betapar, gp)
            if check_beta(tmp_beta, gp):
                gh.LOG(2, 'beta error')
                tmp_profs.chi2 = gh.err(1., gp)
                return tmp_profs
            try:
                if gp.checksig:
                    import gi_analytic as ga
                    anrho = ga.rho(gp.xepol, gp)[0]
                    rhodmpar_half = np.exp(splev(gp.dat.rhalf[0], splrep(gp.xepol, np.log(anrho))))
                    nr = -gh.derivipol(np.log(anrho), np.log(gp.xepol))
                    dlr = np.hstack([nr[0], nr, nr[-1]])
                    if gp.investigate =='gaia':
                        dlr[-1] = 4
                        rhodmpar = np.hstack([rhodmpar_half, dlr])
                    lbaryonpar = 0.0*rhodmpar
                    MtoL = 0.0
                    if gp.investigiiate == 'gaia':
                        if gp.case == 5:
                            betapar = np.array([  4.24378376e-14,   1, 2, 1.41421356e+02])
                        elif gp.case == 6:
                            betapar = np.array([  5.82811566e-13,   1, 2, 3.53553391e+02])
                        elif gp.case == 7:
                            betapar = np.array([  8.31604545e-25,  -1.69605095e-29,   1, 1])
                    elif gp.investigate == 'hern':
                        betapar = np.array([0, 0, 2, max(gp.xipol)/2]) # for hern
                    annu = ga.rho(gp.xepol, gp)[1]
                    nupar_half = np.exp(splev(gp.dat.rhalf[1], splrep(gp.xepol, np.log(annu))))
                    nrnu = -gh.derivipol(np.log(annu), np.log(gp.xepol))
                    dlrnu = np.hstack([nrnu[0], nrnu, nrnu[-1]])
                    if gp.investigate == 'gaia':
                        dlrnu[-1] = 6
                    nupar = np.hstack([nupar_half, dlrnu])
                sig,kap,zetaa,zetab=phys.sig_kap_zet(gp.xepol, rhodmpar, lbaryonpar, MtoL, nupar, betapar, pop, gp)

            except Exception:
                gh.LOG(1, 'sigma error')
                tmp_profs.chi2 = gh.err(2., gp)
                return tmp_profs

            # now store the profiles
            gh.sanitize_vector(tmp_beta, len(tmp_profs.x0), -200, 1, gp.debug)
            tmp_profs.set_prof('beta', tmp_beta, pop, gp)
            gh.sanitize_vector(tmp_betastar, len(tmp_profs.x0), -1, 1, gp.debug)
            tmp_profs.set_prof('betastar', tmp_betastar, pop, gp)
            tmp_profs.set_prof('sig', sig, pop, gp)
            tmp_profs.set_prof('kap', kap, pop, gp)
            tmp_profs.set_zeta(zetaa, zetab, pop)

        tmp_profs.set_prof('nrnu', tmp_nrnu, pop, gp)
        tmp_profs.set_prof('nu', tmp_nu, pop, gp) # pool: tmp_nu.get()

        # following profile needs to be stored at all times, to calculate chi
        tmp_profs.set_prof('Sig', tmp_Signu, pop, gp)

        off += offstep # still do this even if gp.chi2_Sig_converged is False
    if off != gp.ndim:
        gh.LOG(1, 'wrong subscripts in gi_loglike')
        pdb.set_trace()

    # determine log likelihood

    chi2 = calc_chi2(tmp_profs, gp)
    gh.LOG(-1, '  ln L = ', gh.pretty(-chi2/2.))
    # x=gp.dat.rbin
    # linedat,=ax.loglog(x, gp.dat.Sig[1], 'b')
    # line,=ax.loglog(x, tmp_profs.get_prof("Sig", 1), 'r', alpha=0.1)
    # plt.draw()
    # plt.show()
    tmp_profs.chi2 = chi2

    # after some predefined wallclock time and Sig convergence, plot all profiles
    if time.time() - gp.last_plot >= gp.plot_after and gp.chi2_Sig_converged <= 0:
        gp.last_plot = time.time()
        try:
            import plotting.plot_profiles
            plotting.plot_profiles.run(gp.files.timestamp, gp.files.outdir, gp)
        except:
            print('plotting error in gi_loglike!')

    # close pool automatically after with clause
    return tmp_profs
## \fn geom_loglike(cube, ndim, nparams, gp)
# define log likelihood function to be called by Multinest and plot_profiles
# spherical version
# @param cube parameter cube as defined by gi_class_cube, in physical space already (not [0,1] cube anymore)
# @param ndim number of dimensions, needed as argument by multinest
# @param nparams number of parameters, needed as argument by multinest
# @param gp global parameters
