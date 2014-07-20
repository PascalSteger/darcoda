#!/usr/bin/env ipython3

##
# @file
# all functions to work with plots

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import pdb
import matplotlib.pyplot as plt


def startlog():
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.set_yscale('log')
    return fig, ax
## \fn startlog
# start new figure in semilog. to be used in debug mode


def startloglog():
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')
    return fig, ax
## \fn startloglog
# start new figure in loglog. to be used in debug mode


def start():
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    # fig.clf()
    # fig.set_yscale('linear')
    return fig, ax
## \fn start()
# start new figure. to be used in debug mode


def setlims(ax, xlim, ylim):
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    return
## \fn setlims(ax, xlim, ylim)
# set limits on axis ax
# @param ax axis object
# @param xlim [xmin, xmax]
# @param ylim [ymin, ymax]


def save_plot():
    plt.savefig(gp.files.get_outpng())
    return
## \fn save_plot()
# save plot to png file
# NOT USED ANYMORE
