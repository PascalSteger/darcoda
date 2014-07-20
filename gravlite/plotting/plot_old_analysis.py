#!/usr/bin/env python3

## @file
# plot all profiles of a given run

# (c) 2013 ETHZ Pascal Steger, psteger@phys.ethz.ch

# start off with plotting/ in front of path
#import sys
import pdb, pickle
import numpy as np
import numpy.random as npr
npr.seed(1989) # for random events that are reproducible

import matplotlib
import matplotlib.pyplot as plt
plt.ioff()

import gl_collection as glc
import gl_loglike as gll

calculate_anew = False



def calculate_profiles(gp):
    models = read_models(basename)
    pc = glc.ProfileCollection(gp.pops, gp.nipol)
    i = 0
    for m in models:
        i += 1
        print('model '+str(i)+' of '+str(len(models)))
        tmp_profs = gll.geom_loglike(m, 0, 0, gp)
        pc.add(tmp_profs)
    pc.cut_subset()
    # or: use [minchi, maxchi] for whole range
    # or: use [maxchi*0.99, maxchi*1.01] for only highest chi2
    return pc
## \fn calculate_profiles(gp)
# calculate all profiles once, then store them in ascii output files
# NOT USED ANYMORE
# @param gp global profiles


def pcsave(basename, pc):
    with open(basename+'pc.save', 'wb') as f:
        pickle.dump(pc, f)
## \fn pcsave(basename, pc)
# save profile collection
# NOT USED ANYMORE
# @param basename string for directory
# @param pc profile collection


def pcload(basename):
    with open(basename+'pc.save', 'rb') as f:
        return pickle.load(f)
## \fn pcload(basename)
# NOT USED ANYMORE
# load saved profile collection
# @param basename string directory



def read_scale(basename, gp):
    gp.Rscale = []
    gp.Sig0pc = []
    gp.maxvlos = []
    for i in range(gp.pops+1):
        # basename + 'scale_' + str(i) + '.txt'
        A = np.loadtxt(gp.files.get_scale_file(i),\
                       unpack=False, skiprows=1)
        gp.Rscale.append(A[0])  # [pc]
        gp.Sig0pc.append(A[2])  # [Munit/pc^2]
        # totmass is A[2] # [Munit]
        gp.maxvlos.append(A[4]) # [km/s]
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
    # ALL = np.vstack([REJECTED[1000:1100,:-2]])
    return ALL
## \fn read_models(basename)
# read in all models, concatenate them
# @param basename string


def calculate_profiles(gp):
    models = read_models(basename)
    pc = glc.ProfileCollection(gp.pops, gp.nipol)
    i = 0
    for m in models:
        i += 1
        print('model '+str(i)+' of '+str(len(models)))
        tmp_profs = gll.geom_loglike(m, 0, 0, gp)
        pc.add(tmp_profs)
    pc.cut_subset()
    # or: use [minchi, maxchi] for whole range
    # or: use [maxchi*0.99, maxchi*1.01] for only highest chi2
    return pc
## \fn calculate_profiles(gp)
# calculate all profiles once, then store them in ascii output files
# @param gp global profiles


def pcsave(basename, pc):
    with open(basename+'pc.save', 'wb') as f:
        pickle.dump(pc, f)
## \fn pcsave(basename, pc)
# save profile collection
# @param basename string for directory
# @param pc profile collection


def pcload(basename):
    with open(basename+'pc.save', 'rb') as f:
        return pickle.load(f)
## \fn pcload(basename)
# load saved profile collection
# @param basename string directory


if __name__ == '__main__':
    import select_run as sr
    timestamp, basename = sr.run()

    # include runtime gl_params, but other files all from current directory
    import import_path as ip

    # load stored parameters
    ip.insert_sys_path(basename)
    import gl_params as glp
    ip.remove_first()

    # uncomment following to use stored collection, loglike, and all other modules as well
    #ip.insert_sys_path(basename+'sphere')
    #import gl_collection as glc
    #import gl_loglike as gll
    ##ip.remove_first(); ip.remove_first() # uncomment to include most recent representation code from class_cube

    gp = glp.Params(timestamp)
    import gl_file as glf
    gp.dat = glf.get_data(gp)

    gp.pops = sr.get_pops(basename)
    print('working with ', gp.pops, ' populations')

    read_scale(basename, gp) # store half-light radii in  gp.Rscale
    import gl_helper as gh
    Radii, Binmin, Binmax, Sigdat1, Sigerr1 = gh.readcol5(gp.files.Sigfiles[0]) # [Rscale0], [Munit/Rscale0^2]
    gp.xipol = Radii * gp.Rscale[0]       # [pc] TODO check
    maxR = max(Radii) # [pc]
    Radii = np.hstack([Radii, 2*maxR, 4*maxR, 8*maxR])
    gp.xepol = Radii * gp.Rscale[0]       # [pc] TODO check

    if calculate_anew:
        pc = calculate_profiles(gp)
        pcsave(basename, pc)
    else:
        pc = pcload(basename)
    pc.cut_subset()
    pc.set_x0(gp.xipol) # [pc]

    if gp.investigate =='walk' or gp.investigate=='gaia':
        r0analytic = np.logspace(np.log10(1.), np.log10(max(gp.xepol)), 100)
        pc.set_analytic(r0analytic, gp)

    pc.calculate_M(gp)
    pc.sort_profiles(gp)

    pc.write_all(basename, gp)

    pc.plot_profile(basename, 'chi2', 0, gp)
    pc.plot_profile(basename, 'rho', 0, gp)
    pc.plot_profile(basename, 'M', 0, gp)
    pc.plot_profile(basename, 'nr', 0, gp)
    pc.plot_profile(basename, 'betastar', 1, gp)
    pc.plot_profile(basename, 'beta', 1, gp)
    pc.plot_profile(basename, 'Sig', 1, gp)
    pc.plot_profile(basename, 'sig', 1, gp)
    if gp.pops == 2:
        pc.plot_profile(basename, 'betastar', 2, gp)
        pc.plot_profile(basename, 'beta', 2, gp)
        pc.plot_profile(basename, 'Sig', 2, gp)
        pc.plot_profile(basename, 'sig', 2, gp)
        
    # pc.plot_overview(basename, gp)
    # plt.draw()
