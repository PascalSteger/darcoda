#!/usr/bin/python2.7
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
print ('GravLite: A Non-Parametric Mass Modelling Routine for discs and spheres')

import numpy as np
import gl_params as gp
import gl_funs   as gfun
import gl_file   as gfile
import gl_plot   as gpl
import gl_init   as ginit
import gl_priors as gprio
import pdb
from gl_class_params import Params

if gp.getnewdata: gfile.bin_data()

gfile.get_data()
gfile.ipol_data()

ginit.mcmc_init()                       # TODO: thread

gpl.prepare_plots()                     # TODO: thread
gpl.plot_data()


n = 0
plotfirst = True
while ( n < gp.niter-1):
    if not gp.initphase:
        gfile.store_old_params(gp.pars,gp.chi2)
        
    if not gp.checksigma:
        gfun.get_new_parameters() # or: gp.parst.assign(
        if gprio.check_density():         continue
        if gprio.check_mass():            continue
        if gprio.check_delta():           continue


    # only calculate sigma and other variables if previous test succeeded
    if gp.geom == 'sphere':
        gfun.calc_M_nu_sig_sphere()
    elif gp.geom == 'disc':
        gfun.calc_M_nu_sig_disc()

    if not gp.checksigma:
        if gprio.check_sigma(): continue

    gfun.calc_chi2()

    if plotfirst:
        gpl.plot_first_guess()
        plotfirst = False

    gfun.accept_reject(n)       # includes adapt_stepsize
    if gp.checksigma: break

    if gp.metalpop and np.random.rand()<0.001: # get new data for another population split
        # every 1000th model passing through priors (exclude 100 burn-in)
        gfile.bin_data();        gfile.get_data();        gfile.ipol_data()
        gpl.plot_data();         gpl.plot_first_guess()
        # TODO: increase stepsize again?
        
    # gfun.wait()
    if not gp.initphase:
        gfile.write_outfile()
        n = n + 1


gp.LOG.warning( 'Finished!' )
gpl.show_plots()
