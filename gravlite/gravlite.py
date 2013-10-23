#!/usr/bin/env python3

##
# @file
# main file
# @ingroup gravlite
# @defgroup gravlite all functions called from gravlite
#

##
# @package gravlite Non-Parametric Mass Modelling Routine for Discs and Spheres
# (c) ETHZ 2013, Pascal Steger, psteger@phys.ethz.ch
#

# or for local snowball usage, from ipython in emacs: !/usr/bin/env ipython-python3.2

print('GravLite: Non-Parametric Mass Modelling Routine for discs and spheres')
print('(c) 2013 Pascal S.P. Steger, psteger@phys.ethz.ch')



import numpy as np
import gl_params as gp
import gl_funs   as gfun
import gl_file   as gfile
import gl_plot   as gpl
import gl_init   as ginit
import gl_priors as gprio
import pdb
from gl_class_params import Params


## run MCMC, called from __init__ after setting up all variables and plotting area
def run_MCMC():
    ## counter variable for the MCMC
    n = 0 
    ## bool to determine whether the first plot has to be drawn
    plotfirst = True 

    while ( n < gp.niter-1):
        if not gp.initphase:
            gfile.store_old_params(gp.pars,gp.chi2)

        if not gp.checkint:
            gfun.get_new_parameters() # or: gp.parst.assign() for special wishes
            if gprio.check_density():   continue
            if gprio.check_mass():      continue
            if gprio.check_delta():     continue

        # only calculate sigma and other variables if previous tests succeeded
        try:
            if gp.geom == 'sphere':
                    gfun.calc_M_nu_sig_kap_sphere()
            elif gp.geom == 'disc':
                    gfun.calc_M_nu_sig_disc()
        except Exception as detail:
            print('handling error in calc_M_nu_sig_kap_sphere:', detail)
            # gave back old, working values, have to get new parameters now, so jump to end of loop
            continue

        if not gp.checkint:
            if gprio.check_sigma(): continue

        try:
            gfun.calc_chi2() # determine likelihood function (*not* reduced chi2)
        except Exception as detail:
            print('handling error in calc_chi2')
            continue

        if plotfirst:
            gpl.plot_first_guess()                  # plot the first model
            plotfirst = False

        gfun.accept_reject(n) # accept/reject new pars, adapt_stepsize handling
        if gp.checkint: break # only 1 iteration if checkint set

        if gp.initphase=='over' and gp.metalpop and np.random.rand()<0.0001:
            # get new data for another population split
            # every 10000th model passing through priors (exclude 100 burn-in)
            gfile.bin_data();        gfile.get_data();        gfile.ipol_data()
            gpl.plot_data();         gpl.plot_first_guess()
            # TODO: increase stepsize again?

        # gfun.wait()
        if not gp.initphase:
            gfile.write_outfile()
            n = n + 1                       # n is increased even if model is rejected!

    return



if __name__=="__main__":
    if gp.getnewdata:    gfile.bin_data()
    gpl.prepare_plots()                     # TODO: thread
    gfile.get_data()
    gfile.ipol_data() 
    # interpolate to regular (lin or log) array, or keep stuff if gp.consttr
    
    ginit.mcmc_init()                       # TODO: thread
    gpl.plot_data()
    
    run_MCMC()

    print('Finished!')
    gpl.show_plots()


## @mainpage Gravlite - Non-Parametric Mass Modelling Routine for Discs and Spheres
#
# @section intro_sec Introduction
#
# Gravlite is a tool to determine the mass distribution in disc-like and spherical systems.
# It takes as input a tracer density distribution, a line-of-sight velocity dispersion, and the velocity's fourth moment as a function of radius.
# It then generates a highly-dimensional parameter space for tracer density, overall density distribution, and possibly a velocity anisotropy profile ands traces it with a simple Markov Chain Monte Carlo method.
#
# @section install_sec Installation
# Following packets need to be installed on your system:
# * python3
# * matplotlib/pylab
# * scipy
#
# Unzip the file gravlite.zip, and run
#
# python3 gravlite.py
# 
# @subsection pars Parameter files: Main configuration
# Sample parameter files for several possible scenarios are stored in the subfolder ./params.
# The file gl_params.py is a soft link to one of them.
# Following mass modelling methods have been implemented so far:
# * gl_params_walker.py: spherical Walker mock data from the Gaia challenge catalogue, 2 populations
# * gl_params_gaia.py:   spherical mock data from the Gaia challenge catalogue, 1 population
# * gl_params_hern.py:   spherical mock data taken from a Hernquist profile
# * gl_params_simple.py: disk-like mock data, generated on the fly
# * gl_params_sim.py:    disk-like mock data, from a simulation by S. Garbari
# 
# @subsection init Init values
# Initial values for all profiles and the respective stepsizes are set in gl_init.py. We use following representations:
# * overall density rho: can either be set as rho_i=rho(r_i), where r_i are the radii of bins, or the coefficients of the Legendre polynomial;
# * tracer densities nu_i = nu(r_i) in linear space or nu_i = log(nu(r_i)) in logarithmic space;
# * velocity anisotropy beta_i are incremental changes from beta(r=0) = 0.
# 
# Have fun!
