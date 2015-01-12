#!/usr/bin/env python3

## @file
# plot all profiles of a given run
# modified by Hamish Silverwood for batch running, 24 October 2014

# (c) 2014 ETHZ Pascal Steger, psteger@phys.ethz.ch

# start off with plotting/ in front of path
#import sys
import pdb
import pickle
import os
import sys
import numpy as np
import numpy.random as npr
import time
npr.seed(int(time.time())) # 1989 for random events that are reproducible
from optparse import OptionParser
import gl_helper as gh

import matplotlib
matplotlib.use('Agg')

def prepare_output_folder(basename):
    os.system('mkdir -p '+ basename + 'output/data/')
    os.system('mkdir -p '+ basename + 'output/ascii/')
    os.system('mkdir -p '+ basename + 'output/pdf/')
    os.system('mkdir -p '+ basename + 'output/analytic/')
    return 0
## \fn prepare_output_folder(basename)
# create directory structure in output folder
# @param basename string, data folder


#def read_scale(basename, gp):
#    gp.Xscale = []
#    gp.Sig0pc = []
#    gp.maxsiglos = []
#    for pop in range(gp.pops+1):
#        # basename + 'scale_' + str(i) + '.txt'
#        A = np.loadtxt(gp.files.get_scale_file(pop),\
#                       unpack=False, skiprows=1)
#        gp.Xscale.append(A[0])  # [pc]
#        gp.Sig0pc.append(A[1])  # [Munit/pc^2]
#        # totmass_tracers is A[2] # [Munit]
#        gp.nu0pc.append(A[3])   # [Munit/pc^3]
#        gp.maxsiglos.append(A[4]) # [km/s]
### \fn read_scale(basename, gp)
## read scale file, store into gp.*scale
## @param basename string
## @param gp


def correct_E_error(filename):
    import os
    os.system("sed -i 's/\\([0-9]\\)-\\([0-9]\\)/\\1E-\\2/g' " + filename)
    os.system("sed -i 's/\\([0-9]\\)+\\([0-9]\\)/\\1E+\\2/g' " + filename)
    return
## \fn correct_E_error(filename)
# replace 3.4230210-301 and ...+301 by 3.4230210E-301 and ...+301 the sed way
# @param filename string


def read_models(basename):
    # read in all accepted models
    print(basename+'/{ev.dat, phys_live.points}')
    correct_E_error(basename + '/ev.dat')
    correct_E_error(basename + '/phys_live.points')
    REJECTED = np.loadtxt(basename+'/ev.dat', skiprows=0, unpack=False)
    LIVE = np.loadtxt(basename+'/phys_live.points', skiprows=0, unpack=False)
    ALL = np.vstack([REJECTED[:,:-3], LIVE[:,:-2]])
    npr.shuffle(ALL)
    # for debugging, random 10 models
    # ALL = ALL[:10]
    # for debugging, based on live points only
    # ALL = np.vstack([LIVE[:12,:-2]])
    # for debugging, fixed 100 models after N=1000 iterations
    #ALL = np.vstack([REJECTED[:6,:-2]])
    return ALL
## \fn read_models(basename)
# read in all models, concatenate them
# @param basename string


def pcload_single_entries(basename, gp):
    import gl_collection as glc
    pc = glc.ProfileCollection(gp.pops, gp.nipol)
    with open(basename+'pc2.save', 'rb') as fi:
        dum = pickle.load(fi) # dummy variable, was used to create file
        while 1:
            try:
                MODEL = pickle.load(fi)
                if npr.random() < 0.001: #HS bodge to only pull in a fraction of models
                    pc.add(MODEL)
            except EOFError:
                break
    return pc
## \fn pcload_single_entries(basename, gp)
# load all data into [chi^2, profiles] pairs, and add to a profile collection
# @param basename string
# @param gp global parameters


def run(timestamp, basename, gp):
    prepare_output_folder(basename)
    import gl_file as glf
    gp.dat = glf.get_binned_data_noscale(gp)
    import gl_helper as gh
    bincenters, binmins, binmaxs, nudat, nuerr = gh.readcol5(gp.files.nufiles[0])
    max_z = max(bincenters) #?
    min_z = min(bincenters) #?

    pc = pcload_single_entries(basename, gp)
    if len(pc.chis) == 0:
        gh.LOG(1, 'no profiles found for plotting')
        return
    # first plot all chi^2 values in histogram
    pc.plot_profile(basename, 'chi2', 0, gp)

    # then select only the best models for plotting the profiles
    pc.cut_subset()
    pc.set_x0(gp.z_bincenters) # [kpc]
    print('here')

    pc.sort_profiles_disc(gp)
    pc.write_all_disc(basename, gp)
    pc.plot_profile(basename, 'nu_vec', 0, gp)
    pc.plot_profile(basename, 'sig_vec', 0, gp)
    pc.plot_profile(basename, 'rho_DM_vec', 0, gp)
    pc.plot_profile(basename, 'kz_rho_DM_vec', 0, gp)
    pc.plot_profile(basename, 'kz_nu_vec', 0, gp)
    pc.plot_profile(basename, 'Sig_DM_vec', 0, gp)
## \fn run(timestamp, basename, gp)
# call all model read-in, and profile-plotting routines
# @param timestamp string
# @param basename string
# @param gp global parameters


if __name__ == '__main__':
    print(sys.argv)
    basename = sys.argv[1]
    timestamp = sys.argv[1][-13:-1]
    print('Timestamp = ', timestamp)

    import select_run as sr

    # include runtime gl_params, but other files all from current directory
    import import_path as ip
    # load stored parameters
    ip.insert_sys_path(basename+'programs/')
    import gl_params as glp
    ip.remove_first()
    # uncomment following to use stored collection, loglike, all other modules
    #ip.insert_sys_path(basename+'sphere')
    #import gl_collection as glc
    ##ip.remove_first(); ip.remove_first() # uncomment to include most recent
    gp = glp.Params(timestamp)
    gp.pops = sr.get_pops(basename)
    print('working with ', gp.pops, ' populations')
    run(timestamp, basename, gp)
