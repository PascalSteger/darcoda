#!/usr/bin/env python3

##
# @file
# all functions to work with plots

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import pdb
from pylab import *
from matplotlib.ticker import MaxNLocator

import gl_physics as phys
import gl_chi as gc
from gl_analytic import *
from gl_int import *

global f, axs


def startlog():
    fig = figure(figsize=(1,1))
    fig.clf()
    fig.yscale('log')
    return fig
## \fn startlog
# start new figure in semilog. to be used in debug mode


def startloglog():
    clf()
    xscale('log')
    yscale('log')
    return
## \fn startloglog
# start new figure in loglog. to be used in debug mode


def start():
    fig = figure(figsize=(2,2))
    fig.clf()
    # fig.set_yscale('linear')
    return fig
## \fn start()
# start new figure. to be used in debug mode


def setlims(ax,xlim,ylim):
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
