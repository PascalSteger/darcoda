#!/usr/bin/env ipython

##
# @file
# store profiles to files

# (c) 2014 ETHZ Pascal S.P. Steger


import pdb
import csv
import numpy as np

class Output:
    def __init__ (self):
        self.descriptions = []
        self.arrays = []
        return
    ## \fn __init__ (self)
    # constructor, with descriptor and data in array format


    def add(self, d, vec):
        self.descriptions.append(d)
        self.arrays.append(vec)
    # \fn add(d, vec)
    # add a column to ASCII data file
    # @param d string of description
    # @param vec np.array of floats


    def write_old(self, filename):
        with open(filename, 'w') as csvfile:
            csvw = csv.writer(csvfile, delimiter=',')
            for column in zip(*[s for s in self.container]):
                csvw.writerow(column)
        return 0
    ## \fn write_old(self, filename)
    # old routine for output to file
    # @param filename string


    def write(self, filename):
        headers = ",".join(self.descriptions)
        csvfile = open(filename, 'wb')
        np.savetxt(csvfile, np.transpose(np.array(self.arrays)), header=headers, delimiter=",")
        csvfile.close()
    ## \fn write(self, filename)
    # write file to csv file
    # @param filename string


    def __repr__(self):
        return "Output: "+self.descriptions+" "+self.arrays
    ## \fn __repr__(self)
    # string representation for ipython


## \class Output
# class for generating ASCII output files from plot_multinest
