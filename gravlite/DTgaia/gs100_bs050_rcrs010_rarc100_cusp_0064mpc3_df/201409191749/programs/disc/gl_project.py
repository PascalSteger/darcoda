#!/usr/bin/env ipython3

##
# @file
# Functions related to projection and deprojection of density in spherical models.
# Conventions:
# rho, r, Mr     denote 3D density, 3D radius, M(<3D radius)
# Sig, R, MR     denote 2D density, 2D radius, M(<2D radius)
# *SUM*          denotes main method = summing
# *INT*          denotes main method = integrating
# *NORM*         denotes main method = renormalization

# (c) 2013 Pascal Steger, ETH Zurich, psteger@phys.ethz.ch

import numpy as np
import pdb
from scipy.integrate import simps
from scipy.integrate import quad, fixed_quad, quadrature, romberg, cumtrapz
from scipy.interpolate import splrep, splev, splint
import gl_helper as gh
import gl_plot as gpl
import gl_physics as phys


def rho_param_INT_Sig_disc(z0, rhopar, pop, gp):
    # TODO: check integration for z direction only

    # use splines on variable transformed integral
    # \Sigma(R) = \int_{r=0}^{R} \rho(r) dr
    xmin = z0[0]/30. # tweaked. z0[0]/1e4 gives error in quad()
    z0left = np.array([xmin, z0[0]*0.25, z0[0]*0.50, z0[0]*0.75])
    z0nu = np.hstack([z0left, z0])

    rhonu = phys.rho(z0nu, rhopar, pop, gp) # rho takes rho(rhalf) and n(r) parameters
    Sig = np.zeros(len(z0nu)-gp.nexp)
    for i in range(len(z0nu)-gp.nexp):
        # TODO speed up using the same function for integrating all parts
        Sig[i] = gh.quadinflog(z0nu, rhonu, xmin, z0nu[i])

    gh.checkpositive(Sig, 'Sig in rho_param_INT_Sig_disc')
    return Sig[len(z0left):] # @z0 (z0nu without z0left, and without 3 extension bins)
## \fn rho_param_INT_Sig_disc(z0, rhopar, pop, gp)
# take 3D density parameters, calculate projected surface density
# @param z0 radii of bins, (nrho-nexp entries) [pc]
# @param rhopar 3D density, (nrho entries) [Munit/pc^3]
# @param pop int for population (0 both, 1, 2, ..)
# @param gp global parameters


def nu_param_INT_Sig_disc(z0, nupar, pop, gp):
    # TODO: check integration for z direction only
    # use splines on variable transformed integral
    # \Sigma(R) = \int_{r=0}^{R} \rho(r) dr
    z0nu = z0

    nunu = phys.nu_decrease(z0nu, nupar, gp)
    Sig = np.zeros(len(z0nu))
    for i in range(len(z0nu)):
        Sig[i] = gh.quadinflog(z0nu, nunu, z0nu[0], z0nu[i])

    gh.checkpositive(Sig, 'Sig in nu_param_INT_Sig_disc')
    return Sig
## \fn nu_param_INT_Sig_disc(z0, nupar, pop, gp)
# take 3D density parameters, calculate projected surface density
# @param z0 radii of bins, (nrho-nexp entries) [pc]
# @param nupar 3D density, (nrho entries) [Munit/pc^3]
# @param pop int for population (0 both, 1, 2, ..)
# @param gp global parameters

    


