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
import time

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

    nu_1_data = np.loadtxt(basepath_ts +'/nu/nu_0.txt', skiprows=1)
    gp.z_bincenters = nu_1_data[:,0]
    gp.z_binmins = nu_1_data[:,1]
    gp.z_binmaxs = nu_1_data[:,2]
    gp.z_all_pts = np.append([0.0], gp.z_bincenters)

    correct_E_error(basepath_ts + paracube_filename)

    paracube_data = np.loadtxt(basepath_ts + paracube_filename, dtype='a')

    gp.plotting_flag = True

    print('Converting ', len(paracube_data[:,0]), ' paracube points to profiles.')
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

    nu_0_data = np.loadtxt(basepath_ts +'/nu/nu_0.txt', skiprows=1)
    gp.z_bincenters = nu_0_data[:,0]
    gp.z_binmins = nu_0_data[:,1]
    gp.z_binmaxs = nu_0_data[:,2]
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



def write_mn_info(gp):
    #Writes info file for multinest post analysis
    # H Silverwood 19/10/2015
    darcoda_path = gh.find_darcoda_path()
    ts = gp.files.timestamp

    info_file_name = darcoda_path +'/gravimage/DT'+ gp.investigate +'/'+ str(gp.case) +'/'+ ts +'/'+ 'output.info'
    info_file = open(info_file_name, 'w+')

    #Run information
    info_file.writelines('# GravImage Multinest Output Info \n')
    info_file.writelines('# Start:' + time.strftime("%d/%m/%Y %H:%M:%S") + '\n')
    info_file.writelines('# Time stamp = ' + ts +'\n')
    info_file.writelines('#\n')
    info_file.writelines('# number dimensions = ' + str(gp.ndim) +'\n')
    info_file.writelines('# number tracer pops = ' + str(gp.ntracer_pops) +'\n')
    info_file.writelines('#\n')

    #sigz constant 'C' parameters
    lc=1 #line count
    info_file.writelines('# sigc constant C parameter(s) \n')
    for t_pop in range(0, gp.ntracer_pops):
        info_file.writelines('lab' + str(lc) + '= C_' + str(t_pop) + '\n')
        lc+=1

    #Dark matter parameters
    info_file.writelines('# Dark Matter Parameters \n')
    if gp.darkmattermodel == 'const_dm':
        info_file.writelines('#     DM model: const_dm \n')
        info_file.writelines('lab' + str(lc) + '= \\rho_{\\rm DM, const} \n')
        lc+=1
    elif gp.darkmattermodel == 'ConstPlusDD':
        info_file.writelines('#     DM model: ConstPlusDD \n')
        info_file.writelines('lab' + str(lc) + '= \\rho_{\\rm DM, const} \n')
        lc+=1
        info_file.writelines('lab' + str(lc) + '= K_{\\rm DD} \n')
        lc+=1
        info_file.writelines('lab' + str(lc) + '= D_{\\rm DD} \n')
        lc+=1
    elif gp.darkmattermodel == 'kz_dm':
        info_file.writelines('#     DM model: kz_dm \n')
        info_file.writelines('lab' + str(lc) + '= k_{z,C} \n')
        lc+=1
        for jter in range(0, gp.nrhonu):
            info_file.writelines('lab' + str(lc) + '= k_{z,' + str(jter) + '} \n')
            lc+=1

    #Baryon mass profile parameters
    info_file.writelines('# Baryon Model Parameters \n')
    for baryon_pop in range(0, gp.nbaryon_pops):
        if gp.baryonmodel == 'simplenu_baryon':
            info_file.writelines('#     Baryon model: simplenu \n')
            info_file.writelines('lab' + str(lc) + '= K_{\\rm baryon} \n')
            lc+=1
            info_file.writelines('lab' + str(lc) + '= D_{\\rm baryon} \n')
            lc+=1

    #Tracer profile parameters: nu_C, kz_nu_C, kz_nu_vector # kz_nu_LS
    info_file.writelines('# Tracer density model parameters \n')
    for tracer_pop in range(0, gp.ntracer_pops):
        if gp.scan_rhonu_space:
            info_file.writelines('#     Tracer model pop ' + str(tracer_pop) + ': scan_rhonu_space \n')
            raise Exception('Still working on this bit')

        elif gp.nu_model=='kz_nu':
            info_file.writelines('#     Tracer model pop ' + str(tracer_pop) + ': kz_nu \n')
            info_file.writelines('lab' + str(lc) + '= \\nu_C \n')
            lc+=1
            info_file.writelines('lab' + str(lc) + '= k_{\\nu, {\\rm C}} \n')
            lc+=1
            for jter in range(0, gp.nbins[tracer_pop]):
                info_file.writelines('lab' + str(lc) + '= k_{\\nu,' + str(jter) +'} \n')
                lc+=1

        elif gp.nu_model=='gaussian_data':
            info_file.writelines('#     Tracer model pop ' + str(tracer_pop) + ': gaussian_data \n')
            info_file.writelines('lab' + str(lc) + '= \\nu_{\\rm C} \n')
            lc+=1
            for jter in range(0, gp.nbins[tracer_pop]):
                info_file.writelines('lab' + str(lc) + '= \\nu_{' + str(jter) +'} \n')
                lc+=1

    # Introducing tilt term:
    if gp.tilt:
        info_file.writelines('# Tilt term model parameters \n')
        for tracer_pop in range(0, gp.ntracer_pops):
            info_file.writelines('lab' + str(lc) + '= A_{' + str(tracer_pop) + ',{\\rm tilt}} \n')
            lc+=1
            info_file.writelines('lab' + str(lc) + '= n_{' + str(tracer_pop) + ',{\\rm tilt}} \n')
            lc+=1
            info_file.writelines('lab' + str(lc) + '= R_{' + str(tracer_pop) + ',{\\rm tilt}} \n')
            lc+=1

    if lc-1 != gp.ndim:
        raise Exception('Incorrect number of labels in mn_info')

    info_file.writelines('# Number of Parameters saved \n')
    info_file.writelines('params_saved = ' + str(gp.ndim))

    info_file.close()
    return




#if __name__=="__main__":
#    from optparse import OptionParser
#    parser = OptionParser()
#    parser.add_option("-i", "--investigation", dest="investigation", default="")
#    parser.add_option("-c", "--case", dest="case", default=0)
#    parser.add_option("-t", "--timestamp", dest="timestamp", default="")
#    (options, args) = parser.parse_args()
#
##    basepath_ts='/home/hsilverw/LoDaM/darcoda/gravimage/DT' + str(options.investigation) + '/' + str(options.case) + '/' + str(options.timestamp) +'/'
#    basepath_ts='/home/sofia/darcoda/gravimage/DT' + str(options.investigation) + '/' + str(options.case) + '/' + str(options.timestamp) +'/'
#    paracube_to_profile(basepath_ts, "phys_live.point", "phys_live_profiles.save", 'simplenu', options.case, str(options.timestamp))
