#!/usr/bin/env python3

##
# @file
# parameters for MultiNest.
# cube = [0,1]^ndim,  ndim = nepol + 2*pops*(nepol+nbeta), 
# Holds representations for overall density,
# [tracer density, anisotropy]_{for each population}
# Gives access to density, population profiles,
# and calculate physical values from [0,1]
# populations are counted from 0 = first component,
# 1 = first additional component

# (c) 2013 ETHZ Pascal S.P. Steger

import pdb
import csv
import numpy as np
import numpy.random as npr
import random
import gl_params
gp = gl_params.Params()
import gl_physics as phys
import gl_helper as gh

class Output:
    def __init__ (self):
        self.container = []
        return
    ## \fn __init__ (self)
    # constructor, with descriptor and data in array format
    
    def add(self, d, arr):
        self.container.append(np.hstack([d, arr]))
    # \fn add(d, arr)
    # add a column to ASCII data file
    # @param d string of description
    # @param arr np.array of floats

    def write(self, filename):
        with open(filename, 'w') as csvfile:
            csvw = csv.writer(csvfile, delimiter=',')
            for column in zip(*[s for s in self.container]):
                csvw.writerow(column)

        return 0
        
## \class Output
# class for generating ASCII output files from plot_multinest
