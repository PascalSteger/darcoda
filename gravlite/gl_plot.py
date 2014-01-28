#!/usr/bin/env python3

##
# @file
# all functions to work with plots

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import pdb
from pylab import *
from matplotlib.ticker import MaxNLocator

import gl_params as gp
import gl_physics as phys
import gl_chi as gc
from gl_analytic import *
from gl_int import *
import gl_units as units

global f, axs

# TODO: needs function to plot any profile in semilogarithmic plot,
#       with custom labels


# TODO: function for plotting median, 1sigma, 2sigma profiles in semilog plot
def startlog():
    clf()
    yscale('log')
    return
## \fn start_log
# start new figure in semilog. to be used in debug mode


def start():
    clf()
    yscale('linear')
    return
## \fn start()
# start new figure. to be used in debug mode
    

def setlims(ax,xlim,ylim):
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    return
## \fn setlims(ax, xlim, ylim)
# set limits on axis ax
# @param xlim [xmin, xmax]
# @param ylim [ymin, ymax]

def save_plot():
    plt.savefig(gp.files.get_outpng())
    return
## \fn save_plot()
# save plot to png file
