#!/usr/bin/env python3

## @file
# plot all profiles of a given run

# (c) 2013 ETHZ Pascal Steger, psteger@phys.ethz.ch

# start off with plotting/ in front of path
#import sys
import pdb, pickle, os
import numpy as np
import numpy.random as npr
npr.seed(1989) # for random events that are reproducible

import matplotlib
import matplotlib.pyplot as plt
plt.ioff()

calculate_anew = False

def prepare_output_folder(basename):
    os.system('mkdir -p '+ basename + 'output/data/')
    os.system('mkdir -p '+ basename + 'output/ascii/')
    os.system('mkdir -p '+ basename + 'output/png/')
    os.system('mkdir -p '+ basename + 'output/pdf/')
    os.system('mkdir -p '+ basename + 'output/analytic/')
    return 0
## \fn prepare_output_folder(basename)
# create directory structure in output folder
# @param basename string, data folder


def read_scale(basename, gp):
    gp.Xscale = []
    gp.Sig0pc = []
    gp.maxsiglos = []
    for pop in range(gp.pops+1):
        # basename + 'scale_' + str(i) + '.txt'
        A = np.loadtxt(gp.files.get_scale_file(pop),\
                       unpack=False, skiprows=1)
        gp.Xscale.append(A[0])  # [pc]
        gp.Sig0pc.append(A[1])  # [Munit/pc^2]
        # totmass is A[2] # [Munit]
        gp.nu0pc.append(A[3])   # [Munit/pc^3]
        gp.maxsiglos.append(A[4]) # [km/s]
## \fn read_scale(basename, gp)
# read scale file, store into gp.*scale
# @param basename string
# @param gp


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


def pcload_single_entries(basename):
    import gl_collection as glc
    pc = glc.ProfileCollection(gp.pops, gp.nipol)
    with open(basename+'pc2.save', 'rb') as fi:
        dum = pickle.load(fi) # dummy variable, was used to create file
        while 1:
            try:
                MODEL = pickle.load(fi)
                pc.add(MODEL)
            except EOFError:
                break
    return pc
## \fn pcload_single_entries(basename)
# load all data into [chi^2, profiles] pairs, and add to a profile collection
# @param basename string


if __name__ == '__main__':
    import select_run as sr
    timestamp, basename = sr.run()
    prepare_output_folder(basename)

    # include runtime gl_params, but other files all from current directory
    import import_path as ip

    # load stored parameters
    ip.insert_sys_path(basename+'programs/')
    import gl_params as glp
    ip.remove_first()

    # uncomment following to use stored collection, loglike, and all other modules as well
    #ip.insert_sys_path(basename+'sphere')
    #import gl_collection as glc
    ##ip.remove_first(); ip.remove_first() # uncomment to include most recent representation code from class_cube

    gp = glp.Params(timestamp)
    import gl_file as glf
    gp.dat = glf.get_data(gp)

    gp.pops = sr.get_pops(basename)
    print('working with ', gp.pops, ' populations')

    read_scale(basename, gp) # store half-light radii in  gp.Xscale
    import gl_helper as gh
    Radii, Binmin, Binmax, Sigdat1, Sigerr1 = gh.readcol5(gp.files.Sigfiles[0]) # [Xscale0], [Munit/Xscale0^2]
    gp.xipol = Radii * gp.Xscale[0]       # [pc]
    maxR = max(Radii)                     # [pc]
    Radii = np.hstack([Radii, 2*maxR, 4*maxR, 8*maxR]) # [pc]z
    gp.xepol = Radii * gp.Xscale[0]       # [pc]

    #if calculate_anew:
    #    pc = calculate_profiles(gp)
    #    pcsave(basename, pc)
    #else:
    pc = pcload_single_entries(basename)
    # pc = pcload(basename)
    pc.cut_subset()
    pc.set_x0(gp.xipol) # [pc]

    if gp.investigate =='walk' or gp.investigate=='gaia':
        r0analytic = np.logspace(np.log10(1.),\
                                 np.log10(max(gp.xepol)), 100)
        pc.set_analytic(r0analytic, gp)

    pc.sort_profiles(gp)

    pc.write_all(basename, gp)

    pc.plot_profile(basename, 'chi2', 0, gp)
    pc.plot_profile(basename, 'rho', 0, gp)

    if gp.investigate == 'obs':
        pc.plot_profile(basename, 'Sig', 0, gp)
        pc.plot_profile(basename, 'nu', 0, gp)
        pc.plot_profile(basename, 'nrnu', 0, gp)

    pc.plot_profile(basename, 'nr', 0, gp)
    pc.plot_profile(basename, 'M', 0, gp)

    for pop in np.arange(1, gp.pops+1):
        pc.plot_profile(basename, 'betastar', pop, gp)
        pc.plot_profile(basename, 'beta', pop, gp)
        pc.plot_profile(basename, 'Sig', pop, gp)
        pc.plot_profile(basename, 'nu', pop, gp)
        pc.plot_profile(basename, 'nrnu', pop, gp)
        pc.plot_profile(basename, 'sig', pop, gp)

    # pc.plot_overview(basename, gp)
