#!/usr/bin/env python3

## Interfaces with multinest output files

# Hamish Silverwood, GRAPPA, UvA, 15 January 2015

import numpy as np
import gl_helper as gh
import pdb
import pickle
import sys
sys.path.insert(0, 'plotting/')
from plot_profiles import correct_E_error

def paracube_to_profile(basepath_ts, paracube_filename, profile_filename, investigation, case, timestamp):
    # Takes a file with points defined in the multinest cube, output file with profiles
    # as per the pc2.save type


    #Pull in parameters from the run basepath
    import import_path as ip
    # load stored parameters
    ip.insert_sys_path(basepath_ts+'programs/')
    import gl_params as glp
    ip.remove_first()
    gp = glp.Params(timestamp, investigation, int(case))

    import gl_class_cube as glcc
    from gl_loglike import geom_loglike
    import gl_file as glf

    # Import the bin centres used for the run
    gp.dat = glf.get_binned_data_noscale(gp)

    nu_1_data = np.loadtxt(basepath_ts +'/nu/nu_1.txt', skiprows=1)
    gp.z_bincenters = nu_1_data[:,0]
    gp.z_binmins = nu_1_data[:,1]
    gp.z_binmaxs = nu_1_data[:,2]
    gp.z_all_pts = np.append([0.0], gp.z_bincenters)

    correct_E_error(basepath_ts + paracube_filename)

    paracube_data = np.loadtxt(basepath_ts + paracube_filename, dtype='a')

    gp.plotting_flag = True

    for iter in range(0, len(paracube_data[:,0])):
        #print('Iter = ', iter)
        paracube = paracube_data[iter,:]
        paracube = paracube.astype('float')
        #pdb.set_trace()

        try:
            temp_profile = geom_loglike(paracube, gp.ndim, gp.ndim+1, gp)
        except ValueError:
            continue

        temp_profile.x0 = gp.z_bincenters
        temp_profile.xbins = np.hstack([gp.z_binmins, gp.z_binmaxs[-1]])
        with open(gp.files.outdir+profile_filename, 'ab') as fi:
            pickle.dump(temp_profile, fi)


def mn_output_to_profile(basepath_ts, mn_output_filename, profile_filename, investigation, case, timestamp):
    # Takes a file with points defined in the multinest cube, output file with profiles
    # as per the pc2.save type


    #Pull in parameters from the run basepath
    import import_path as ip
    # load stored parameters
    ip.insert_sys_path(basepath_ts+'programs/')
    import gl_params as glp
    ip.remove_first()
    gp = glp.Params(timestamp, investigation, int(case))

    import gl_class_cube as glcc
    from gl_loglike import geom_loglike
    import gl_file as glf

    # Import the bin centres used for the run
    gp.dat = glf.get_binned_data_noscale(gp)

    nu_1_data = np.loadtxt(basepath_ts +'/nu/nu_1.txt', skiprows=1)
    gp.z_bincenters = nu_1_data[:,0]
    gp.z_binmins = nu_1_data[:,1]
    gp.z_binmaxs = nu_1_data[:,2]
    gp.z_all_pts = np.append([0.0], gp.z_bincenters)

    correct_E_error(basepath_ts + mn_output_filename)
    mn_output_data = np.loadtxt(basepath_ts + mn_output_filename, dtype ='a')

    #paracube_data = np.loadtxt(basepath_ts + paracube_filename, dtype='a')

    gp.plotting_flag = True

    for iter in range(0, len(mn_output_data[:,0])):
        #print('Iter = ', iter)
        mn_weight = mn_output_data[iter, 0].astype('float')
        mn_m2lnL = mn_output_data[iter, 1].astype('float')
        paracube = mn_output_data[iter, 2:]

        paracube = paracube.astype('float')
        #pdb.set_trace()

        try:
            temp_profile = geom_loglike(paracube, gp.ndim, gp.ndim+1, gp)
        except ValueError:
            continue

        temp_profile.x0 = gp.z_bincenters
        temp_profile.xbins = np.hstack([gp.z_binmins, gp.z_binmaxs[-1]])
        temp_profile.mn_weight = mn_weight
        with open(gp.files.outdir+profile_filename, 'ab') as fi:
            pickle.dump(temp_profile, fi)







if __name__=="__main__":
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-i", "--investigation", dest="investigation", default="")
    parser.add_option("-c", "--case", dest="case", default=0)
    parser.add_option("-t", "--timestamp", dest="timestamp", default="")
    (options, args) = parser.parse_args()

#    basepath_ts='/home/hsilverw/LoDaM/darcoda/gravimage/DT' + str(options.investigation) + '/' + str(options.case) + '/' + str(options.timestamp) +'/'
    basepath_ts='/home/sofia/darcoda/gravimage/DT' + str(options.investigation) + '/' + str(options.case) + '/' + str(options.timestamp) +'/'
    paracube_to_profile(basepath_ts, "phys_live.point", "phys_live_profiles.save", 'simplenu', options.case, str(options.timestamp))
