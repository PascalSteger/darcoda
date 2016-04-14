#!/usr/bin/env python3

## @file
# collect profiles and perform actions on them

# (c) 2014 ETHZ Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import numpy.random as npr
import math
import pdb
import scipy

import matplotlib.pyplot as plt
plt.ioff()

#from gl_class_profiles import Profiles
#import plotting.gl_output as go
#import gl_helper as gh
#import gl_analytic as ga
#import gl_project as glp
#import gl_physics as phys

import barrett.posterior as posterior
import barrett.data as data


def plot_barrett_histograms_all(h5_filename, plot_filename, nbins, gp):
    plt.ioff()
    fig, axes = plt.subplots(nrows=gp.ndim,
                         ncols=1,
                         sharex=False,
                         sharey=False)

    plt.subplots_adjust(wspace=0, hspace=1)

    for jter in range(0, gp.ndim):
        ax = axes[jter]
        limits = None

        P = posterior.oneD(h5_filename, gp.param_headers[jter], limits=limits, nbins=nbins)
        P.marginalise()
        P.plot(ax)

        ax.set_yticklabels([])

    fig.set_size_inches(4,gp.ndim*2)
    fig.savefig(plot_filename,  bbox_inches='tight')

    plt.close(fig)

def plot_barrett_histograms_split(h5_filename, plot_filename_root, nbins, gp):
    plt.ioff()
    for jter in range(0, gp.ndim):
        fig = plt.figure(figsize=(4, 1))
        ax  = fig.add_subplot(111)

        limits = None

        P = posterior.oneD(h5_filename, gp.param_headers[jter], limits=limits, nbins=nbins)
        P.marginalise()
        P.plot(ax)

        ax.set_yticklabels([])
        #ax.tick_params(axis='x', labelsize=8)
        #pdb.set_trace()

        if gp.nu_model=='gaussian_data':
            #Parse the param_header
            p_head = gp.param_headers[jter]
            #pdb.set_trace()
            if 'nu' in p_head and 'k_' not in p_head and 'C' not in p_head:
                bin_num = int(''.join([s for s in p_head if s.isdigit()]))
                nu_val = gp.dat.nu[0][bin_num]
                nu_plus_err = nu_val + gp.dat.nuerr[0][bin_num]
                nu_minus_err = nu_val - gp.dat.nuerr[0][bin_num]
                ax.axvline(nu_val, color='k', linestyle='-')
                ax.axvline(nu_plus_err, color='k', linestyle='-', alpha=0.6)
                ax.axvline(nu_minus_err, color='k', linestyle='-', alpha=0.6)
            #elif 'nu' in p_head and 'C' in p_head:
            #    ax.axvline(gp.nu_C_max, color='k', linestyle='-')
            #    ax.axvline(gp.nu_C_min, color='k', linestyle='-')
            #    ax.axvline((gp.nu_C_max + gp.nu_C_min)/2, color='k', linestyle='-')
            #    pdb.set_trace()
            #elif 'C' in p_head and '0' in p_head:
            #    IntC_max=(gp.sigz_C_max**2)*gp.nu_C_max
            #    IntC_min=(gp.sigz_C_min**2)*gp.nu_C_min
            #    ax.axvline(IntC_max, color='k', linestyle='-')
            #    ax.axvline(IntC_min, color='k', linestyle='-')
            #    ax.axvline((IntC_max + IntC_min)/2, color='k', linestyle='-')


            if 'rho' in p_head and 'baryon' in p_head and 'C' in p_head:
                rho_model_mid = gp.gaussian_rho_baryon_mid_vector[0]
                rho_plus_err = rho_model_mid + gp.gaussian_rho_baryon_SD_vector[0]
                rho_minus_err = rho_model_mid - gp.gaussian_rho_baryon_SD_vector[0]
                ax.axvline(rho_model_mid, color='k', linestyle='-')
                ax.axvline(rho_plus_err, color='k', linestyle='-', alpha=0.6)
                ax.axvline(rho_minus_err, color='k', linestyle='-', alpha=0.6)
            elif 'rho' in p_head and 'baryon' in p_head:
                bin_num = int(''.join([s for s in p_head if s.isdigit()]))
                rho_model_mid = gp.gaussian_rho_baryon_mid_vector[bin_num+1]
                rho_plus_err = rho_model_mid + gp.gaussian_rho_baryon_SD_vector[bin_num+1]
                rho_minus_err = rho_model_mid - gp.gaussian_rho_baryon_SD_vector[bin_num+1]
                ax.axvline(rho_model_mid, color='k', linestyle='-')
                ax.axvline(rho_plus_err, color='k', linestyle='-', alpha=0.6)
                ax.axvline(rho_minus_err, color='k', linestyle='-', alpha=0.6)


        fig.savefig(plot_filename_root + str(jter) + '.pdf', bbox_inches='tight')

        plt.close(fig)


def plot_barrett_2D_hist_sep(h5_filename, plot_filename_root, nbins, gp):

    import matplotlib.ticker as ticker
    var1 = '$R_{0,{\\rm tilt}}$'
    var2 = '$\\rho_{\\rm DM, const}$'

    fig = plt.figure(figsize=(5, 5))
    ax  = fig.add_subplot(111)

    limits = None

    P = posterior.twoD(h5_filename, var1, var2, ylimits=limits, nbins=nbins)


    P.marginalise()
    P.plot(ax)
    pdb.set_trace()
    ax.set_xlabel(r'${\rm Tilt}\,\,R_{0} \quad[{\rm kpc}]$', fontsize=14)
    ax.set_ylabel(r'$\rho_{\rm DM, const} \quad[10^{-3}\rm{M}_\odot/\rm{pc}^3]$', fontsize=14)


    kpc = 3.0857E19   # kpc in m
    Msun = 1.9891E30  # Sun's mass in kg
    GeV = 1.78266E-27 # GeV in kg

    ax2 = ax.twinx()

    ax_ylim = np.array([0,20])*1E6    # SAMPLING
    ax_ylim2 = 1E6*Msun/(GeV*(100.*kpc)**3)*ax_ylim*(1E-6)
    ax.set_ylim(ax_ylim)
    ax2.set_ylim(ax_ylim2)

    ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/(1E6)))
    ax.yaxis.set_major_formatter(ticks_y)

    #ticks_y_tex=[r'$0$', r'$5$', r'$10$', r'$15$', r'$20']
    #ax.set_yticks(ticks_y_tex)


    ax2.set_ylabel(r'$\rho_{\rm{DM}}\quad[\rm{GeV}/\rm{cm}^3]$', fontsize=14)


    fig.savefig(plot_filename_root +'_' + var1 + '_' + var2 + '_68_95_99_smoothed.pdf', bbox_inches='tight')

    plt.close(fig)

#def plot_vs_mass(dataset, vars, filename, nbins=60):
#
#    if not os.path.isdir(dataset):
#         os.mkdir(dataset)
#
#    n = len(vars)
#
#    fig, axes = plt.subplots(nrows=n,
#                             ncols=1,
#                             sharex='col',
#                             sharey=False)
#
#    plt.subplots_adjust(wspace=0, hspace=0)
#
#    for i, y in enumerate(vars):
#        ax = axes[i]
#
#        if y[0:6] == 'log(C_':
#            ylimits = (-30, 30)
#        else:
#            ylimits = None
#
#        P = posterior.twoD(dataset+'.h5', 'log(m_{\chi})', y, ylimits=ylimits, nbins=nbins)
#        P.marginalise()
#        P.plot(ax)
#
#    fig.set_size_inches(8,n*6)
#    fig.savefig('%s/%s.png' % (dataset, filename), dpi=100)
#
#    plt.close(fig)
