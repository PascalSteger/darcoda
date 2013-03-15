#!/usr/bin/python
print ('GravLite: A Non-Parametric Mass Modelling Routine')
#print ('A non-parametric method to determine the total enclosed mass of a dwarf spheroidal implemented for two independent tracer populations')

import numpy as np
import gl_params as gp
import gl_funs   as gfun
import gl_file   as gfile
import gl_plot   as gplot
import gl_init   as ginit
import gl_priors as gprio
import pdb
from gl_class_params import Params



gplot.prepare_plots()
gfile.get_data()

gplot.plot_data()

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
    n = n + 1
    if gp.initphase == 'over':
        gfile.store_old_params(gp.pars,gp.chisq)

    if not gp.checksigma:
        gfun.get_new_parameters() # or: gp.parst.set(

    gfun.calc_nu_sig()

    if not gp.checksigma:
        if gprio.check_priors(): continue

    gfun.calc_chisqt()

    if plotfirst:
        gplot.plot_first_guess()
        plotfirst = False

    gfun.accept_reject(n)
    gfun.adapt_stepsize()

    if gp.checksigma: break
    # gfun.wait()
    gfile.write_outfile()

gp.LOG.warning( 'Finished!' )
gplot.show_plots()
