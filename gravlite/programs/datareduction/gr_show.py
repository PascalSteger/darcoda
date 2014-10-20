#!/usr/bin/python3

##
# @file
# draw preprocessed data
from pylab import *
ion()

def show_part_pos(x, y, pmn, Xscale):
    clf()
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
    if gp.investigate == 'walk' or gp.investigate == 'gaia':
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
    savefig(gp.files.dir+'centeredpos_' + str(n) + '.png')
    pdb.set_trace()
    return
## \fn show_part_pos(x, y, pmn, Xscale)
# show 2D scatter plot of particle positions
# @param x coordinate
# @param y coordinate
# @param pmn probability of membership
# @param Xscale scale radius in 2D


def show_plots_dens_3D(rbin, p_dens, p_edens, gp):
    clf()
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
    pdb.set_trace()
## \fn show_plots_dens_3D(Rbin, p_dens, p_edens, gp)
# show density
# @param Rbin
# @param p_dens
# @param p_edens
# @param gp global parameters


def show_plots_dens_2D(Rbin, P_dens, P_edens, Dens0pc):
    clf()
    # plot density
    plot(Rbin, P_dens*Dens0pc, 'b', lw=1)
    lbound = (P_dens-P_edens)*Dens0pc; lbound[lbound<1e-6] = 1e-6
    ubound = (P_dens+P_edens)*Dens0pc
    plot(Rbin, lbound, 'k')
    plot(Rbin, ubound, 'k')
    fill_between(Rbin, lbound, ubound, alpha=0.5, color='r')
    yscale('log')
    xlabel(r'$R [R_c]$')
    ylabel(r'$\nu_{2D}(R) [\mathrm{Munit/pc/pc}]$')
    savefig(gp.files.dir+'Sigma/Sig_' + str(n) + '.png')
    pdb.set_trace()
## \fn show_plots_dens_2D(Rbin, P_dens, P_edens, Dens0pc)
# show density
# @param Rbin bin radii, array, [pc]
# @param P_dens density
# @param P_edens error on density
# @param Dens0pc central density in Munit/pc^3


def show_plots_sigma(Rbin, p_dvlos, p_edvlos):
    clf()
    plot(Rbin, p_dvlos, 'b', lw=1, label='data')
    fill_between(Rbin, p_dvlos-p_edvlos, p_dvlos+p_edvlos, alpha=0.5, color='r')
    # [rscale],2*[km/s]
    xlabel(r'$R [\mathrm{Xscale}]$')
    ylabel(r'$\langle\sigma_{\mathrm{LOS}}\rangle [\mathrm{km/s}]$')
    savefig( gp.files.dir+'siglos/siglos_' + str(n) + '.png')
    pdb.set_trace()
    return
## \fn show_plots_sigma(Rbin, p_dvlos, p_edvlos)
# show sigma
# @param Rbin [pc]
# @param p_dvlos
# @param p_edvlos


def show_plots_vlos(rbin, p_dvlos, p_edvlos):
    clf()
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
    savefig( gp.files.dir+'siglos/siglos_' + str(n) + '.png')
    pdb.set_trace()
## \fn show_plots_vlos(rbin, p_dvlos, p_edvlos)
# show line-of-sight velocity profile with error bars
# @param rbin [pc]
# @param p_dvlos profile, array of floats, defined in the bins given by rbin
# @param p_edvlos error profile



def show_plots_kappa(Rbin, p_kappa, p_ekappa):
    plot(Rbin, p_kappa, 'b', lw=1)
    fill_between(Rbin, p_kappa-p_ekappa, p_kappa+p_ekappa, alpha=0.5, color='r')
    # [rscale], 2*[1]
    xlabel(r'$R [\mathrm{Xscale}]$')
    ylabel(r'$\langle\kappa_{\mathrm{LOS}}\rangle [1]$')
    ylim([0, 5.])
    # xlim([0, gp.maxR])
    savefig( gp.files.dir+'kappalos/kappalos_' + str(n) + '.png')
    pdb.set_trace()
    return
## \fn show_plots_kappa(Rbin, P_dens, P_edens, p_dvlos, p_edvlos, p_kappa, p_ekappa, Dens0pc)
# show kappa profile with errors
# @param Rbin [pc]
# @param p_kappa
# @param p_ekappa
