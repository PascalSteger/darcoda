#!/bin/env ipython3

##
# @file
# global params for data analysis step 1 used to generate input for spherical MCMC

# (c) 2013 ETHZ Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import ipdb
from pylab import *

import gl_params
gp = gl_params.Params()
from gl_class_files import Files
import gl_helper as gh

showplots = False
n = 100 # number of iterations in MCMC_*
# nbins = gp.nipol

Rerr  = 0. #0.01      # distance error in [Xscale]
vrerr = 0. #2.0      # [km/s] 0.01 # velocity error. only raises sig_los

lograd  = False    # log steps for radial bin in readout,
                        # show x-axis in log scale

case = gp.case

if gp.investigate == 'hern':
    repr  = 1     # choose simulation representation
    Rcut = 1.e10  # [Rvir]
    Rmin = 0. * gp.Xscale[0] # [pc]
    Rmax = 3. * gp.Xscale[0] # [pc]

    simname = gp.files.get_sim_name(gp) # dir+'simulation/'+prename+'unit_hern_%i' %(repr)
    simpos = gp.files.dir+'simulation/'+simname+'stars_pos.txt'
    simvel = gp.files.dir+'simulation/'+simname+'stars_vel.txt'


elif gp.investigate == 'walk': # or just want to try some other generic pymc stuff:
    r_DM  = 1000.

    def rhodm(r):
        exp = ((1.*gamma_DM-beta_DM)/alpha_DM)
        rho = rho0*(r/r_DM)**(-gamma_DM)*(1.+(r/r_DM)**alpha_DM)**exp
        return rho

    fi = Files(gp)
    fi.set_walk(gp)
    dir = fi.dir
    fil = dir+'mem2'

    pmsplit = 0.9 # minimum probability of membership required for analysis
    # use 0 if grw_* should be called from within gravlite
    fileposcartesian = dir+'simulation/pos.txt'
    filevelcartesian = dir+'simulation/vel_my.txt'

elif gp.investigate == 'gaia':
    fi  = Files(gp)
    fi.set_gaia(gp)
    dir = fi.dir
    fil = dir + 'dat'
    r_DM = 1000.

elif gp.investigate == 'triax':
    fi  = Files(gp)
    fi.set_triax(gp)
    dir = fi.dir
    fil = dir + 'dat'
    r_DM = 1500.                  # [pc]

elif gp.investigate == 'obs':
    fi = Files(gp)
    fi.set_obs(gp)
    dir = fi.dir
    fil = dir+'mem2'
    pmsplit = 0.9

def volume_circular_ring(Binmin, Binmax, gp):
    Vol = np.zeros(gp.nipol)
    for k in range(gp.nipol):
        Vol[k] = np.pi*(Binmax[k]**2-Binmin[k]**2) # [Rscale^2] or [pc^2]
    return Vol
## \fn volume_circular_ring(Binmin, Binmax, gp)
# volume of a circular ring from binmin to binmax, in [Rscale^2] or [pc^2]
# @param Binmin in [Rscale] or [pc]
# @param Binmax in [Rscale] or [pc]
# @param gp global parameters


def get_com_file(n):
    gh.sanitize_scalar(n, 0, 2)
    return gp.files.dir+'centeredpos_' + str(n) + '.txt'
## \fn get_com_file(n)
# get filename of COM file
# @param n population


def determine_radius(R, Rmin, Rmax, gp):
    if lograd:
        print(gp.nipol, ' bins in log spacings')
        return gh.bin_r_log(Rmax/gp.nipol, Rmax, gp.nipol)
    elif gp.consttr:
        print(len(R)/gp.nipol,' particles per bin')
        return gh.bin_r_const_tracers(R, gp.nipol)
    else:
        print(gp.nipol, ' bins in linear spacings')
        return gh.bin_r_linear(Rmin, Rmax, gp.nipol)
## \fn determine_radius(R, Rmin, Rmax, gp)
# determine radius once and for all. this must not be changed between
# readout and gravlite run. if you wish to change: set gp.getnewdata =
# True in gl_params.py
# @param R
# @param Rmin float
# @param Rmax float
# @param gp global parameters


def show_part_pos(x, y, pmn, Xscale, comp):
    res = (abs(x)<3)*(abs(y)<3)
    x = x[res]; y = y[res]           # [Xscale]
    en = len(x)
    if en == 0:
        return

    scatter(x[:en], y[:en], c=pmn[:en], s=35, \
            vmin=0.95, vmax=1.0, lw=0.0, alpha=0.2)
    # xscale('log'); yscale('log')
    colorbar()
    circ_HL=Circle((0,0), radius=Xscale/Xscale, fc='None', ec='g', lw=1)
    gca().add_patch(circ_HL)
    if gp.investigate != 'obs':
        circ_DM=Circle((0,0), radius=r_DM/Xscale, fc='None', ec='r', lw=1)
        gca().add_patch(circ_DM)

    # visible region
    maxr = max(np.abs(x));  mayr = max(np.abs(y)) # [rscale]
    width2 = max([maxr,mayr]) # [rscale]
    xlim([-width2,width2])
    ylim([-width2,width2])
    axes().set_aspect('equal')

    xlabel(r'$x [R_s]$')
    ylabel(r'$y [R_s]$')
    # title(fil)
    # savefig(gp.files.dir+'centeredpos_' + str(n) + '.png')
    return
## \fn show_part_pos(x, y, pmn, Xscale, comp) show 2D scatter plot of particle
# positions
# @param x coordinate
# @param y coordinate
# @param pmn probability of membership
# @param Xscale scale radius in 2D
# @param comp population number: 0, 1, 2


def show_plots_dens_3D(comp, rbin, p_dens, p_edens, gp):
    ion(); subplot(111)
    plot(rbin, p_dens, 'b', lw=1)
    lbound = p_dens-p_edens; lbound[lbound<1e-6] = 1e-6
    ubound = p_dens+p_edens;
    fill_between(rbin, lbound, ubound, alpha=0.5, color='r')
    yscale('log')
    xlim([0, gp.maxR])
    ylim([np.min(lbound), np.max(ubound)])
    xlabel(r'$r [r_c]$')
    ylabel(r'$\nu(r)/\nu(0)$')
    savefig( gp.files.dir+'siglos/siglos_' + str(n) + '.png')
    ioff(); show(); clf()
## \fn show_plots_dens_3D(Rbin, p_dens, p_edens, gp)
# show density
# @param Rbin
# @param p_dens
# @param p_edens
# @param gp global parameters


def show_plots_dens_2D(comp, Rbin, P_dens, P_edens, Dens0pc):
    ion(); subplot(111)

    # plot density
    plot(Rbin, P_dens*Dens0pc, 'b', lw=1)
    lbound = (P_dens-P_edens)*Dens0pc; lbound[lbound<1e-6] = 1e-6
    ubound = (P_dens+P_edens)*Dens0pc
    plot(Rbin, lbound, 'k')
    plot(Rbin, ubound, 'k')
    fill_between(Rbin, lbound, ubound, alpha=0.5, color='r')
    yscale('log')
    # xlim([0, gp.maxR]); ylim([np.min(lbound),np.max(ubound)])
    xlabel(r'$R [R_c]$')
    ylabel(r'$\nu_{2D}(R) [\mathrm{Munit/pc/pc}]$')
    savefig(gp.files.dir+'Sigma/Sig_' + str(n) + '.png')
    ioff(); show(); clf()
## \fn show_plots_dens_2D(comp, Rbin, P_dens, P_edens, Dens0pc)
# show density
# @param comp int component
# @param Rbin bin radii, array, [pc]
# @param P_dens density
# @param P_edens error on density
# @param Dens0pc central density in Munit/pc^3


def show_plots_sigma(comp, Rbin, p_dvlos, p_edvlos):
    ion(); subplot(111)
    plot(Rbin, p_dvlos, 'b', lw=1)
    fill_between(Rbin, p_dvlos-p_edvlos, p_dvlos+p_edvlos, alpha=0.5, color='r')
    # [rscale],2*[km/s]

    xlabel(r'$R [\mathrm{Xscale}]$')
    ylabel(r'$\langle\sigma_{\mathrm{LOS}}\rangle [\mathrm{km/s}]$')
    ylim([-1, 30])
    # xlim([0, gp.maxR])
    savefig(get_siglos_png(comp))
    ioff(); show(); clf()
    return
## \fn show_plots_sigma(Rbin, p_dvlos, p_edvlos)
# show sigma
# @param comp int
# @param Rbin [pc]
# @param p_dvlos
# @param p_edvlos


def show_plots_vlos(rbin, p_dvlos, p_edvlos):
    ion(); subplot(111)
    print('rbin = ', rbin,' rscale')
    print('p_dvlos = ', p_dvlos,' km/s')
    print('p_edvlos = ', p_edvlos, 'km/s')
    plot(rbin, p_dvlos, 'b', linewidth=3)
    fill_between(rbin, p_dvlos-p_edvlos, p_dvlos+p_edvlos, alpha=0.5, color='r')
    # [rscale], [km/s], [km/s]

    xlabel(r'$r [rscale]$')
    ylabel(r'$\langle\sigma_{LOS}\rangle [km/s]$')
    ylim([-5,30])
    # xscale('log')
    xlim([np.min(rbin),np.max(rbin)])
    savefig(gpr.get_siglos_png(0))
    ioff();show();clf()
## \fn show_plots_vlos(rbin, p_dvlos, p_edvlos)
# show line-of-sight velocity profile with error bars
# @param rbin [pc]
# @param p_dvlos profile, array of floats, defined in the bins given by rbin
# @param p_edvlos error profile



def show_plots_kappa(comp, Rbin, p_kappa, p_ekappa):
    ion(); subplot(111)
    plot(Rbin, p_kappa, 'b', lw=1)
    fill_between(Rbin, p_kappa-p_ekappa, p_kappa+p_ekappa, alpha=0.5, color='r')
    # [rscale], 2*[1]
    xlabel(r'$R [\mathrm{Xscale}]$')
    ylabel(r'$\langle\kappa_{\mathrm{LOS}}\rangle [1]$')
    ylim([0, 5.])
    # xlim([0, gp.maxR])
    savefig( gp.files.dir+'kappalos/kappalos_' + str(n) + '.png')
    ioff(); show(); clf()
## \fn show_plots_kappa(comp, Rbin, P_dens, P_edens, p_dvlos, p_edvlos, p_kappa, p_ekappa, Dens0pc)
# show kappa profile with errors
# @param comp int
# @param Rbin [pc]
# @param p_kappa
# @param p_ekappa
