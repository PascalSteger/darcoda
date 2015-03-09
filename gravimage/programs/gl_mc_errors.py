#!/usr/bin/env python3

## Estimates errors using Monte Carlo sampling

# Hamish Silverwood, GRAPPA, UvA, 23 February 2015

import numpy as np
import gl_helper as gh
import pdb
import pickle
import sys
import numpy.random as rand
import matplotlib.pyplot as plt


#TEST this will eventually go outside
def ErSamp_gauss_linear_w_z():
    fraction_err = 0.05
    datafile = '/home/hsilverw/LoDaM/darcoda/Data_Sets/simplenu/simplenu_sigz_raw_sdz_p05_sdvz_5.dat'
    data = np.loadtxt(datafile)
    z_data = data[:, 0]

    z_sampled = []
    for z_val in z_data:
        z_sampled.append(rand.normal(loc = z_val, scale= z_val*fraction_err))

    return z_sampled


z_data_flat_distro = rand.random(2000000)
def ErSamp_flat_distro_test():
    fraction_err = 0.001
    z_data = z_data_flat_distro

    z_sampled = []
    for z_val in z_data:
        z_sampled.append(abs(rand.normal(loc = z_val, scale = fraction_err)))

    return z_sampled




def mc_nu_error(sampled_z_func, number_mcs, binmin, binmax, bincenter):
    # sampled_z_func - returns a vector of z points

    nu_vectors=[]
    for jter in range(0, number_mcs):
        jter_z_data = sampled_z_func()
        jter_nu, dummy, dummy, dummy, dummy = gh.nu_sig_from_bins(binmin, binmax, jter_z_data, np.ones(len(jter_z_data)))
        nu_vectors.append(jter_nu)

    #Calculate standard deviations of nu
    nu_vectors = np.array(nu_vectors)
    nu_stdevs  = []
    nu_means   = []
    nu_medians = []
    for pter in range(0, len(binmin)):
        nu_stdevs.append(np.std(nu_vectors[:, pter]))
        nu_means.append(np.mean(nu_vectors[:, pter]))
        nu_medians.append(np.median(nu_vectors[:, pter]))

    #pdb.set_trace()
    #fig = plt.figure()
    #ax = fig.add_subplot(111)

    #no, bins, patches = ax.hist(nu_vectors[:,0], 100)

    return np.array(nu_stdevs)



if __name__=="__main__":

    binmin, binmax, bincenter = gh.bin_r_linear(0.2, 0.8, 12)

    nu_stdevs = mc_nu_error(ErSamp_flat_distro_test, 100, binmin, binmax, bincenter)
    pdb.set_trace()
