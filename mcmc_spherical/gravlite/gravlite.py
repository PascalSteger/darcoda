#!/usr/bin/python2.7
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
print ('GravLite: A Non-Parametric Mass Modelling Routine')
#print ('A non-parametric method to determine the total enclosed mass of a dwarf spheroidal implemented for two independent tracer populations')

import numpy as np
import gl_params as gp
import gl_funs   as gfun
import gl_file   as gfile
import gl_plot   as gpl
import gl_init   as ginit
import gl_priors as gprio
import pdb
from gl_class_params import Params

if gp.getnewdata:
    import grw_com
    import grw_dens
    import grw_siglos

gpl.prepare_plots()
gfile.get_data()

gpl.plot_data()

gfile.ipol_data()
ginit.mcmc_init()
gfile.adump()

gp.parst = Params(0)
gp.parst.set(gp.pars)

if gp.logprior: gp.lparst.set(gp.lpars)



n = 0
plotfirst = True
while ( n < gp.niter-1):
    # print('n = ',n)
    if not gp.initphase:
        gfile.store_old_params(gp.pars,gp.chi2)

    if not gp.checksigma:
        gfun.get_new_parameters() # or: gp.parst.set(

    gfun.calc_M_nu_sig()

    if not gp.checksigma:
        if gprio.check_priors(): continue

    gfun.calc_chi2()

    if plotfirst:
        gpl.plot_first_guess()
        plotfirst = False

    gfun.accept_reject(n)
    gfun.adapt_stepsize()

    if gp.checksigma: break
    # gfun.wait()
    gfile.write_outfile()
    if not gp.initphase: n = n + 1

gp.LOG.warning( 'Finished!' )
gpl.show_plots()
