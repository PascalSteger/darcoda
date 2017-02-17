#!/usr/bin/env python3

## @file
# plot all profiles of a given run

# (c) 2014 ETHZ Pascal Steger, psteger@phys.ethz.ch

# start off with plotting/ in front of path
#import sys
#sys.path.append("/home/sofia/darcoda")
import pdb
import pickle
import os
import numpy as np
import numpy.random as npr
import time
npr.seed(int(time.time())) # 1989 for random events that are reproducible
from optparse import OptionParser
import gl_helper as gh
import gl_multinest_helper as glmh
import h5py
import gl_barrett as glb
#import gl_collection as glc

def prepare_output_folder(basename):
    os.system('mkdir -p '+ basename + 'output/data/')
    os.system('mkdir -p '+ basename + 'output/ascii/')
    os.system('mkdir -p '+ basename + 'output/pdf/')
    os.system('mkdir -p '+ basename + 'output/analytic/')
    os.system('mkdir -p '+ basename + 'output/histograms/')
    return 0
## \fn prepare_output_folder(basename)
# create directory structure in output folder
# @param basename string, data folder

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


def pcload_single_entries(basename, profile_source, gp):
    print ('In plot_prof..., in pcload_simgle_entries')
    import gl_collection as glc
    pc = glc.ProfileCollection(gp.ntracer_pops, gp.nbins, gp.no_Sigrho_bins)

    if profile_source =='standard':
        filename = 'pc2.save'
    elif profile_source == 'livepoints':
        filename = 'phys_live_profiles.save'
    elif profile_source == 'MNoutput':
        filename = 'phys_MNoutput_profiles.save'

    model_load_count = 0
    print ('basename:',basename,' filename:',filename)
    with open(basename+filename, 'rb') as fi:
        dum = pickle.load(fi) # dummy variable, was used to create file
        while 1:
            model_load_count +=1
            try:
                MODEL = pickle.load(fi) # TAG:
                #print ('MODEL:',MODEL)  #out: Profiles (disc): 130.848579 etc
                pc.add(MODEL)  # THIS IS WHERE THE add FUNCTION IS CALLED !!!
            except EOFError:
                break
    print('Loaded ', model_load_count, ' models from ', filename)
    return pc
## \fn pcload_single_entries(basename, gp)
# load all data into [chi^2, profiles] pairs, and add to a profile collection
# @param basename string
# @param gp global parameters


def run(timestamp, basename, profile_source, gp):


    prepare_output_folder(basename)
    import gl_file as glf
    gp.dat = glf.get_binned_data_noscale(gp)
    import gl_helper as gh
    bincenters, binmins, binmaxs, nudat, nuerr = gh.readcol5(gp.files.nufiles[0])

    pc = pcload_single_entries(basename, profile_source, gp)
    if len(pc.chis) == 0:
        gh.LOG(1, 'no profiles found for plotting')
        return

    # then select only the best models for plotting the profiles
    #pc.cut_subset()
    #pc.set_z_values(gp.z_bincenter_vecs, gp.z_all_pts_sorted)
    #pc.set_z_values(gp.z_bincenter_vecs, gp.extend_z_binc_vecs[0])# Bodge TAG
    # Below: if having different binning for nu
    pc.set_z_values(gp.z_bincenter_vecs, gp.nu_z_bincenter_vecs, gp.z_vec_Sigrho)# extend_z_binc_vecs not used for nu (only for rhoDM etc).
    # Think about which consequences above line has !!!   TAG !
    # Last arg of func call needs to have same dim as used DM etc profiles
    # (Maybe skip connection with bin centers and just divide in 100 or so bins)
    #  (then also skip vertical lines in plots for these...)
    # Just below: tried stuff which doesn't work
    #temp = np.sort(np.concatenate((gp.extend_z_binc_vecs[0],gp.z_bincenter_vecs[1])))
    #print ('temp:',temp)
    #pc.set_z_values(gp.z_bincenter_vecs, temp)# Bodge TAG


    #Plot Histograms WORK IN PROGRESS  # SS commented out below
    #glb.plot_barrett_2D_hist_sep(basename+'mn_output.h5', basename+'/output/histograms/2D_posterior_hist', 60, gp)
    #glb.plot_barrett_histograms_all(basename+'mn_output.h5', basename+'/output/histograms/posterior_hist.pdf', 60, gp)
    #glb.plot_barrett_histograms_split(basename+'mn_output.h5', basename+'/output/histograms/posterior_hist', 60, gp)

    # first plot all chi^2 values in histogram
    #pc.plot_profile(basename, 'chi2', 0, profile_source, gp)  # BODGE !!



    if profile_source == 'MNoutput':
        pc.weighted_sort_profiles_disc(gp)
    elif profile_source == 'livepoints':
        pc.sort_profiles_disc(gp)

    pc.write_all_disc(basename, gp)

    for t_pop in range(0, gp.ntracer_pops):
        pc.plot_profile(basename, 'chi2_sigz2_vecs', t_pop, profile_source, gp)
        pc.plot_profile(basename, 'chi2_sigRz2_vecs', t_pop, profile_source, gp)

        if gp.nu_model == 'kz_nu':
            pc.plot_profile(basename, 'kz_nu_vecs', t_pop, profile_source, gp)
        pc.plot_profile(basename, 'nu_vecs', t_pop, profile_source, gp)
        pc.plot_profile(basename, 'sigz2_vecs', t_pop, profile_source, gp)
        if gp.tilt:
            pc.plot_profile(basename, 'sigmaRz2_vecs', t_pop, profile_source, gp)




    if gp.darkmattermodel == 'kz_dm':
        pc.plot_profile(basename, 'kz_rho_DM_vec', 0, profile_source, gp)
    pc.plot_profile(basename, 'rho_DM_vec', 0, profile_source, gp)
    pc.plot_profile(basename, 'Sig_DM_vec', 0, profile_source, gp)

    if gp.baryonmodel not in ['simplenu_baryon', 'trivial_baryon', 'obs_baryon', 'simplenu_baryon_gaussian']:
        gh.LOG(1, 'No baryon model, all mass is in DM.')
        return

    pc.plot_profile(basename, 'rho_baryon_vec', 0, profile_source, gp)
    pc.plot_profile(basename, 'Sig_baryon_vec', 0, profile_source, gp)

    pc.plot_profile(basename, 'rho_total_vec', 0, profile_source, gp)
    pc.plot_profile(basename, 'Sig_total_vec', 0, profile_source, gp)

    #Output DM limits
    LocalDM_95hi = pc.M95hi.get_prof('rho_DM_vec', 0)[0]
    LocalDM_68hi = pc.M68hi.get_prof('rho_DM_vec', 0)[0]
    LocalDM_medi = pc.Mmedi.get_prof('rho_DM_vec', 0)[0]
    LocalDM_68lo = pc.M68lo.get_prof('rho_DM_vec', 0)[0]
    LocalDM_95lo = pc.M95lo.get_prof('rho_DM_vec', 0)[0]
    print('Local Dark Matter Limits x10^-3 Msun/pc^3')
    print('    Median = ', LocalDM_medi, ' x10^-3 Msun/pc^3')
    print('    68pc CL = [', LocalDM_68lo, LocalDM_68hi, '] x10^-3 Msun/pc^3' )
    print('    95pc CL = [', LocalDM_95lo, LocalDM_95hi, '] x10^-3 Msun/pc^3' )
    print('    95pc CL width = ', LocalDM_95hi-LocalDM_95lo, 'x10^-3 Msun/pc^3')
    print('    68pc CL width = ', LocalDM_68hi-LocalDM_68lo, 'x10^-3 Msun/pc^3')


    print('Local Dark Matter Limits, GeV/cm^3')
    confac = 37.98E-3 #Conversion factor
    print('    Median = ', LocalDM_medi*confac, ' GeV/cm^3')
    print('    68pc CL = [', LocalDM_68lo*confac, LocalDM_68hi*confac, '] GeV/cm^3' )
    print('    95pc CL = [', LocalDM_95lo*confac, LocalDM_95hi*confac, '] GeV/cm^3' )
    print('    95pc CL width = ', (LocalDM_95hi-LocalDM_95lo)*confac, ' GeV/cm^3')
    print('    68pc CL width = ', (LocalDM_68hi-LocalDM_68lo)*confac, ' GeV/cm^3')

    #print(gp.nbins[0], gp.nbins[0], gp.nbins[0], gp.nbins[0], gp.nbins[0])
    #print(LocalDM_95lo, LocalDM_68lo, LocalDM_medi, LocalDM_68hi, LocalDM_95hi)
    print(str(LocalDM_95lo) + ', ' + str(LocalDM_68lo) + ', ' +  str(LocalDM_medi) + ', ' + str(LocalDM_68hi) + ', ' + str(LocalDM_95hi))


## \fn run(timestamp, basename, gp)
# call all model read-in, and profile-plotting routines
# @param timestamp string
# @param basename string
# @param gp global parameters


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", "--investigation", dest="investigate",
                      default="", help="investigation to plot")
    parser.add_option("-c", "--case", dest="case",
                      default=-1, help="case to plot")
    parser.add_option("-l", "--latest", help='plot latest one', dest =
                      'latest', default = False, action = 'store_true')
    #parser.add_option("-t", "--timestamp", dest="timestamp",
    #                  default=-1, help="timestamp of run to plot. Overrides -l")
    #parser.add_option("-a", "--action", dest="action",
    #                  default="p", help="action to take: p: print, k: kill")
    (options, args) = parser.parse_args()
    gh.LOG(1, 'plot_profiles '+str(options.investigate)+' '+str(options.case)+' '+str(options.latest))
    import select_run as sr
    timestamp, basename, investigate, case, profile_source = sr.run(options.investigate, \
                                 options.case,\
                                 options.latest)

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

    gp = glp.Params(timestamp, investigate, case) # Change back! TODO FIXME !!!
    gp.pops = sr.get_pops(basename)
    glmh.write_mn_info(gp)
    print('working with ', gp.pops, ' populations')

    if profile_source == 'livepoints':
        try:
            open(basename+"phys_live_profiles.save")
        except OSError:
            gh.LOG(0, 'No phys_live_profiles.save file found, generating from livepoints now')
            glmh.paracube_to_profile(basename, "outputphys_live.points", "phys_live_profiles.save", investigate, case, timestamp)
            gh.LOG(0, 'phys_live_profiles.save generation completed.')

    if profile_source == 'MNoutput':
        try:
            open(basename+"phys_MNoutput_profiles.save")
        except OSError:
            gh.LOG(0, 'No phys_MNoutput_profiles.save file found, generating from MultiNest output now')
            glmh.mn_output_to_profile(basename, "output.txt", "phys_MNoutput_profiles.save", investigate, case, timestamp)
            gh.LOG(0, 'phys_MNoutput_profiles.save generation completed.')

        try:
            temp = h5py.File(basename+"mn_output.h5", 'r')
        except OSError:
            gh.LOG(0, 'No HDF5 file mn_output.h5 file found, generating from MultiNest output now')
            glmh.mn_output_to_hdf5(basename, "output.txt", "mn_output.h5", investigate, case, timestamp, gp)
            #gh.LOG(0, 'Adding baryon profiles to HDF5 file now')
            #glmh.mn_h5_baryon_vecs(basename, "output.txt", "mn_output.h5", investigate, case, timestamp, gp)

            gh.LOG(0, 'mn_output.h5 generation completed.')

    run(timestamp, basename, profile_source, gp)
