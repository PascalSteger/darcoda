#!/usr/bin/env ipython3

## @file
# define log likelihood function to be called by multinest
# spherical version

# (c) GPL v3 2015 ETHZ Pascal Steger, pascal@steger.aero

import numpy as np
import pdb
from scipy.interpolate import splev, splrep
from pylab import *
#from multiprocessing import Pool

# import matplotlib.pyplot as plt
# fig, ax=plt.subplots()
# import matplotlib.animation as animation

import gi_physics as phys
from gi_class_profiles import Profiles
from gi_priors import check_beta
from gi_chi import calc_chi2
import gi_helper as gh
import gi_project as gip

def geom_loglike(cube, ndim, nparams, gp):
    tmp_profs = Profiles(gp.pops, gp.nepol)
    off = 0
    offstep = gp.nrho
    if gp.chi2_Sig_converged <= 0:
        rhodmpar = np.array(cube[off:off+offstep])
        tmp_rho0 = phys.rho(gp.xepol, rhodmpar, 0, gp)
        # for J factor calculation (has been deferred to output routine)
        #tmp_rhofine = phys.rho(gp.xfine, rhodmpar, 0, gp)
        #tmp_Jfine = gip.Jpar(gp.xfine, tmp_rhofine, gp) #tmp_rhofine
        #tck = splrep(gp.xfine[:-3], tmp_Jfine)
        #tmp_J = splev(gp.xepol, tck)
        # rhodmpar hold [rho(rhalf), nr to be used for integration
        # from halflight radius, defined on gp.xepol]
        # (only calculate) M, check
        tmp_M0 = gip.rho_SUM_Mr(gp.xepol, tmp_rho0)
         # store profiles
        tmp_profs.set_prof('nr', 1.*rhodmpar[1+1:-1], 0, gp)
        tmp_profs.set_prof('rho', tmp_rho0, 0, gp)
        #tmp_profs.set_prof('J', tmp_J, 0, gp)
        tmp_profs.set_prof('M', tmp_M0, 0, gp)
    off += offstep # anyhow, even if Sig not yet converged

    # get profile for rho*
    if gp.investigate == 'obs':
        offstep = gp.nrho
        lbaryonpar = np.array(cube[off:off+offstep])
        rhostar = phys.rho(gp.xepol, lbaryonpar, 0, gp)
        off += offstep
        Signu = gip.rho_param_INT_Sig(gp.xepol, lbaryonpar, 0, gp) # [Munit/pc^2]
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
        tmp_Signu = gip.rho_param_INT_Sig(gp.xepol, nupar, pop, gp)
        #tmp_nu = pool.apply_async(phys.rho, [gp.xepol, nupar, pop, gp])
        #tmp_Signu = pool.apply_async(gip.rho_param_INT_Sig, [gp.xepol, nupar, pop, gp])
        off += offstep

        offstep = 1
        tmp_hyperSig = cube[off:off+offstep]
        off += offstep

        offstep = 1
        tmp_hypersig = cube[off:off+offstep]
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
            #if True:
                if gp.checksig and gp.investigate == 'hern':
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
                    betapar = np.array([0, 0, 2, max(gp.xipol)/2]) # for hern
                    annu = ga.rho(gp.xepol, gp)[1]
                    nupar_half = np.exp(splev(gp.dat.rhalf[1], splrep(gp.xepol, np.log(annu))))
                    nrnu = -gh.derivipol(np.log(annu), np.log(gp.xepol))
                    dlrnu = np.hstack([nrnu[0], nrnu, nrnu[-1]])
                    if gp.investigate == 'gaia':
                        dlrnu[-1] = 6
                    nupar = np.hstack([nupar_half, dlrnu])
                elif gp.checkbeta and gp.investigate == 'gaia':
#                    rhodmpar = np.array([ 0.41586608, 0.38655515, 0.60898657, 0.50936769, 0.52601378, 0.54526758,  0.5755599, 0.57900806, 0.60252357, 0.60668445, 0.62252721, 0.63173754, 0.64555439, 0.65777175, 0.67083556, 0.68506606, 0.69139872, 0.66304763, 0.61462276, 0.70916575, 0.53287872])
                    rhodmpar = np.array([ 0.18235821,  0.4719348,   0.,          0.,          0.10029569,  0.11309553,  0.25637863,  0.31815175,  0.40621336,  0.46247927,  0.53545415,  0.60874961,  0.68978141,  0.79781574,  0.91218048,  1.08482356,  1.36074895,  1.88041885,  2.31792908,  2.62089078,  3.001     ])

                    betapar = np.array([  1.23555034e-03,   9.89999994e-01,   2.03722518e+00,   5.85640906e+00])
                    nupar =  np.array([ 0.15649498,  6.65618254,  0.10293663,  0.1087109,   0.13849277,  0.24371261, 0.62633345,  1.05913181,  1.43774113,  1.82346043,  2.20091446,  2.60007997,  2.98745825,  3.423104,    3.80766658,  4.2089698,   4.62950843,  4.91166037,  4.97380638,  4.99718073,  5.2277589 ])
                    gp.dat.nrnu = [np.array([ 0.15476906,  0.85086798,  0.9342867 ,  0.88161169,  0.83254241, 0.85086798,  0.99930431,  1.22211638,  1.47184763,  1.78910057, 2.1987677 ,  2.51961046,  2.80345393,  3.10336133,  3.88504346,  4.52442727,  4.88817769,  5.07880404,  4.83455511,  6.32165657,  4.88817769]),
               np.array([ 0.15476906,  0.85086798,  0.9342867 ,  0.88161169,  0.83254241,  0.85086798,  0.99930431,  1.22211638,  1.47184763,  1.78910057, 2.1987677 ,  2.51961046,  2.80345393,  3.10336133,  3.88504346,  4.52442727,  4.88817769,  5.07880404,  4.83455511,  6.32165657,  4.88817769]),
               np.array([ 0.15476906,  0.85086798,  0.9342867 ,  0.88161169,  0.83254241, 0.85086798,  0.99930431,  1.22211638,  1.47184763,  1.78910057,  2.1987677 ,  2.51961046,  2.80345393,  3.10336133,  3.88504346,  4.52442727,  4.88817769,  5.07880404,  4.83455511,  6.32165657,    4.88817769]),
               np.array([ 0.15476906,  0.85086798,  0.9342867 ,  0.88161169,  0.83254241, 0.85086798,  0.99930431,  1.22211638,  1.47184763,  1.78910057,  2.1987677 ,  2.51961046,  2.80345393,  3.10336133,  3.88504346,  4.52442727,  4.88817769,  5.07880404,  4.83455511,  6.32165657,  4.88817769])]
                    gp.dat.nrnuerr = [np.array([  0.05158969,  12.22044422,   2.44408884,   2.44408884, 2.44408884,   2.44408884,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777, 0.48881777,   2.44408884,   2.44408884,   2.44408884,   2.44408884]),
                  np.array([  0.05158969,  12.22044422,   2.44408884,   2.44408884, 2.44408884,   2.44408884,   0.48881777,   0.48881777,    0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777, 0.48881777,   2.44408884,   2.44408884,   2.44408884,   2.44408884]),
                  np.array([  0.05158969,  12.22044422,   2.44408884,   2.44408884, 2.44408884,   2.44408884,   0.48881777,   0.48881777, 0.48881777,   0.48881777,   0.48881777,   0.48881777,  0.48881777,   0.48881777,   0.48881777,   0.48881777, 0.48881777,   2.44408884,   2.44408884,   2.44408884,   2.44408884]),
                  np.array([  0.05158969,  12.22044422,   2.44408884,   2.44408884, 2.44408884,   2.44408884,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777, 0.48881777,   2.44408884,   2.44408884,   2.44408884,   2.44408884])]

                    lbaryonpar = 0.0*rhodmpar
                    MtoL = 0.0

                sig,kap,zetaa,zetab = phys.sig_kap_zet(gp.xepol, rhodmpar, lbaryonpar, MtoL, nupar, betapar, pop, gp)
                #fill_between(gp.xipol, gp.dat.sig[1]-gp.dat.sigerr[1], gp.dat.sig[1]+gp.dat.sigerr[1])
                #plot(gp.xepol, sig, 'r')
                #xscale('log')
                #ylim([0, 30])
                #xlabel('$r$ [pc]')
                #ylabel('$\sigma_{LOS}$ [km/s]')
                #savefig('siglos_gaia_2.pdf')
                #pdb.set_trace()
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
            tmp_profs.hypersig = tmp_hypersig
            tmp_profs.set_prof('kap', kap, pop, gp)
            tmp_profs.set_zeta(zetaa, zetab, pop)

        tmp_profs.set_prof('nrnu', tmp_nrnu, pop, gp)
        tmp_profs.set_prof('nu', tmp_nu, pop, gp) # pool: tmp_nu.get()

        # following profile needs to be stored at all times, to calculate chi
        tmp_profs.set_prof('Sig', tmp_Signu, pop, gp)
        tmp_profs.hyperSig = tmp_hyperSig

        off += offstep # still do this even if gp.chi2_Sig_converged is False
    if off != gp.ndim:
        gh.LOG(1, 'wrong subscripts in gi_loglike')
        pdb.set_trace()

    # determine log likelihood
    chi2 = calc_chi2(tmp_profs, gp)
    gh.LOG(-1, gp.investigate+'/'+str(gp.case)+'/'+gp.files.timestamp+':  ln L = ', gh.pretty(-chi2/2.))
    # x=gp.dat.rbin
    # linedat,=ax.loglog(x, gp.dat.Sig[1], 'b')
    # line,=ax.loglog(x, tmp_profs.get_prof("Sig", 1), 'r', alpha=0.1)
    # plt.draw()
    # plt.show()
    tmp_profs.chi2 = chi2

    # after some predefined wallclock time and Sig convergence, plot all profiles
    #if time.time() - gp.last_plot >= gp.plot_after and gp.chi2_Sig_converged <= 0:
    #    gp.last_plot = time.time()
    #    try:
    #        import plotting.plot_profiles
    #        plotting.plot_profiles.run(gp.files.timestamp, gp.files.outdir, gp)
    #    except:
    #        print('plotting error in gi_loglike!')
    # close pool automatically after with clause
    return tmp_profs
## \fn geom_loglike(cube, ndim, nparams, gp)
# define log likelihood function to be called by Multinest and plot_profiles
# spherical version
# @param cube parameter cube as defined by gi_class_cube, in physical space already (not [0,1] cube anymore)
# @param ndim number of dimensions, needed as argument by multinest
# @param nparams number of parameters, needed as argument by multinest
# @param gp global parameters
