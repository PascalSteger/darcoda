#!/usr/bin/env python3

## @file
# plot all profiles of a given run

# (c) 2015 ETHZ Pascal Steger, pascal@steger.aero

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

import gi_helper as gh
import select_run as sr

def prepare_output_folder(basename):
    os.system('mkdir -p '+ basename + 'output/data/')
    os.system('mkdir -p '+ basename + 'output/ascii/')
    os.system('mkdir -p '+ basename + 'output/pdf/')
    os.system('mkdir -p '+ basename + 'output/analytic/')
    return 0
## \fn prepare_output_folder(basename)
# create directory structure in output folder
# @param basename string, data folder

def read_scale(gp):
    gp.Xscale = []
    gp.Sig0pc = []
    gp.maxsiglos = []
    for pop in range(gp.pops+1):
        A = np.loadtxt(gp.files.get_scale_file(pop), unpack=False, skiprows=1)
        gp.Xscale.append(A[0])  # [pc]
        gp.Sig0pc.append(A[1])  # [Munit/pc^2]
        # totmass_tracers is A[2] # [Munit]
        gp.nu0pc.append(A[3])   # [Munit/pc^3]
        gp.maxsiglos.append(A[4]) # [km/s]
## \fn read_scale(gp)
# read scale file, store into gp.*scale
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
    print(basename+'{ev.dat, phys_live.points}')
    correct_E_error(basename + 'ev.dat')
    correct_E_error(basename + 'phys_live.points')
    REJECTED = np.loadtxt(basename+'ev.dat', skiprows=0, unpack=False)
    LIVE = np.loadtxt(basename+'phys_live.points', skiprows=0, unpack=False)
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
    import gi_collection as glc
    pc = glc.ProfileCollection(gp.pops, gp.nepol)
    import re
    tmp = re.split('/DT', basename)[-1]
    path = str.join('/', re.split('/', tmp)[:-1])
    # get number of iterations from run_info file
    import gi_base as gb
    bp = gb.get_basepath()
    fil = open(bp+"/run_info", "r")
    for line in fil:
        if re.search(path, line):
            line2 = re.sub(r'\n', '', line)
            if not re.search("File ", line2):
                runparams = line2
    fil.close()
    numofmodels = int(re.split('\t', runparams)[2])
    current = 0
    with open(basename+'pc2.save', 'rb') as fi:
        dum = pickle.load(fi) # dummy variable, was used to create file
        try:
            while True:
                # and current<22360: # if needing to cut to max no. iterations
                current += 1
                if current%100 == 0:
                    gh.progressbar((1.0*current)/numofmodels)
                MODEL = pickle.load(fi)
                pc.add(MODEL)
        except EOFError:
            pass
    print("")
    return pc
## \fn pcload_single_entries(basename, gp)
# load all data into [chi^2, profiles] pairs, and add to a profile collection
# @param basename string
# @param gp global parameters

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", "--investigation", dest="investigate", default="", help="investigation to plot")
    parser.add_option("-c", "--case",          dest="case",        default=-1, help="case to plot")
    parser.add_option("-l", "--latest", help='plot latest one', dest = 'latest', default = False, action = 'store_true')
    parser.add_option("-t", "--timestamps", help='take all models from following timestamps', dest='timestamps', default=[], action='append')
    # parser.add_option("-t", "--timestamp", dest="timestamp", default=-1, help="timestamp of run to plot. Overrides -l")
    # parser.add_option("-a", "--action", dest="action", default="p", help="action to take: p: print, k: kill")
    (options, args) = parser.parse_args()

    # gh.LOG(1, 'plot_profiles '+str(options.investigate)+' '+str(options.case)+' '+str(options.latest))
    if len(options.timestamps) == 0:
        import select_run as sr
        timestamp, basename = sr.run(options.investigate, options.case, options.latest)
        options.timestamps=[timestamp]

    else:
        import gi_base as gb
        basepath = gb.get_basepath()
        basename = os.path.abspath(basepath+'DT'+options.investigate+'/'+str(options.case))+'/'


    # use last timestamp for sample gp
    import import_path as ip
    tt = options.timestamps[-1]
    ip.insert_sys_path(basename+tt+'/programs/')
    ip.insert_sys_path(basename+tt+'/programs/sphere')
    import gi_params as gip
    gp = gip.Params(tt)
    gp.pops = sr.get_pops(basename+tt+'/')

    # find minimum number of iterations for all runs
    minlinelen = 1e99
    #import gi_collection as gc
    for tt in options.timestamps:
        linelen = 0
        with open(basename+tt+'/pc2.save', 'rb') as fi:
            dum = pickle.load(fi) # dummy variable
            try:
                while True:
                    linelen += 1
                    MODEL = pickle.load(fi)
            except EOFError:
                pass
        print('linelen = ',linelen)
        minlinelen = min(minlinelen, linelen)
    ip.remove_first(); ip.remove_first()
    ip.remove_first(); ip.remove_first()


    # set up overall profile collection
    import import_path as ip
    ip.insert_sys_path(basename+tt+'/programs/')
    ip.insert_sys_path(basename+tt+'/programs/sphere')
    import gi_params as gip
    gp = gip.Params(tt)
    gp.pops = sr.get_pops(basename+tt+'/')
    print('working with ', gp.pops, ' populations')
    ip.remove_first(); ip.remove_first() # programs/sphere and programs/reducedata from gip.Params() call
    import gi_collection
    pcall = gi_collection.ProfileCollection(gp.pops, gp.nepol)
    ip.remove_first(); ip.remove_first() # tt/programs/sphere and tt/programs

    # load all profiles of all timestamps into pcall
    for tt in options.timestamps:
        # call external program to work with independent gp
        os.system('external_pc.py '+basename+' '+tt+' '+str(minlinelen))
        # TODO check: read in pickled values
        with open(basename+tt+'/pc', 'rb') as fn:
            pc = pickle.load(fn)

        pcall.merge(pc)

    # use last timestamp for output
    import import_path as ip
    ip.insert_sys_path(basename+tt+'/programs/')
    ip.insert_sys_path(basename+tt+'/programs/sphere')
    import gi_params as gip
    gp = gip.Params(tt)
    gp.pops = sr.get_pops(basename+tt+'/')
    print('working with ', gp.pops, ' populations')

    if len(pcall.chis) == 0:
        gh.LOG(1, 'no profiles found for plotting')
        sys.exit(1)
    # first plot all chi^2 values in histogram
    pcall.plot_profile(basename+tt+'/', 'chi2', 0, gp)
    # then select only the best models for plotting the profiles
    pcall.cut_subset()
    import gi_helper as gh
    read_scale(gp) # store half-light radii in  gp.Xscale
    Radii, Binmin, Binmax, Sigdat1, Sigerr1 = gh.readcol5(gp.files.Sigfiles[0])
    # [Xscale0], [Munit/Xscale0^2]
    # verified that indeed the stored files in the run directory are used
    gp.xipol = Radii * gp.Xscale[0]       # [pc]
    maxR = max(Radii)                     # [pc]
    minR = min(Radii)                     # [pc]
    Radii = np.hstack([minR/8, minR/4, minR/2, Radii, 2*maxR, 4*maxR, 8*maxR])
    gp.xepol = Radii * gp.Xscale[0]       # [pc]

    pcall.set_x0(gp.xepol, Binmin*gp.Xscale[0], Binmax*gp.Xscale[0]) # [pc]
    #if gp.investigate == 'obs':
    #pc.calculate_J(gp)

    if gp.investigate =='walk' or gp.investigate=='gaia' or gp.investigate =='triax':
        r0analytic = np.logspace(np.log10(1.), np.log10(max(gp.xepol)), 100)
        pcall.set_analytic(r0analytic, gp)
    pcall.sort_profiles(gp)
    pcall.write_all(basename+tt+'/', gp)

    # start following plotting routines as threads, in parallel
    #from multiprocessing import Pool
    #with Pool(processes=8) as pool:
    #    pr=[]
    #    pr.append(pool.apply_async(pc.plot_profile, [basename+tt+'/', 'rho', 0, gp]))
    #for lp in pr:
    #    lp.get()

    pcall.plot_profile(basename+tt+'/', 'rho', 0, gp)
    pcall.plot_profile(basename+tt+'/', 'nr', 0, gp)
    #pc.plot_profile(basename+tt+'/', 'J', 0, gp)
    if gp.investigate == 'obs':
        pcall.plot_profile(basename+tt+'/', 'Sig', 0, gp)
        pcall.plot_profile(basename+tt+'/', 'nu', 0, gp)
        pcall.plot_profile(basename+tt+'/', 'nrnu', 0, gp)
    if gp.geom == 'sphere':
        pcall.plot_profile(basename+tt+'/', 'M', 0, gp)
    for pop in np.arange(1, gp.pops+1):
        pcall.plot_profile(basename+tt+'/', 'betastar', pop, gp)
        pcall.plot_profile(basename+tt+'/', 'beta', pop, gp)
        pcall.plot_profile(basename+tt+'/', 'Sig', pop, gp)
        pcall.plot_profile(basename+tt+'/', 'nu', pop, gp)
        pcall.plot_profile(basename+tt+'/', 'nrnu', pop, gp)
        pcall.plot_profile(basename+tt+'/', 'sig', pop, gp)

    cmd = "cd "+basename+tt+"/output/pdf/;"
    if gp.pops == 1:
        cmd += "pdfjam --outfile cat.pdf --nup 3x4 --no-landscape prof_nr_0.pdf prof_rho_0.pdf prof_M_0.pdf prof_nrnu_1.pdf prof_nu_1.pdf prof_Sig_1.pdf prof_betastar_1.pdf prof_beta_1.pdf prof_sig_1.pdf ../prof_chi2_0.pdf;"
    elif gp.pops == 2:
        cmd += "pdfjam --outfile cat.pdf --nup 3x5 --no-landscape prof_nr_0.pdf prof_rho_0.pdf prof_M_0.pdf prof_nrnu_1.pdf prof_nu_1.pdf prof_Sig_1.pdf prof_betastar_1.pdf prof_beta_1.pdf prof_sig_1.pdf prof_nrnu_2.pdf prof_nu_2.pdf prof_Sig_2.pdf prof_betastar_2.pdf prof_beta_2.pdf prof_sig_2.pdf;"

    cmd += "evince cat.pdf &"
    os.system(cmd)
