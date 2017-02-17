#!/usr/bin/env python3

## Interfaces with multinest output files

# Hamish Silverwood, GRAPPA, UvA, 15 January 2015

import numpy as np
import gl_helper as gh
import pdb
import pickle
import sys
sys.path.insert(0, 'plotting/')
sys.path.append("/home/sofia/darcoda")
from plot_profiles import correct_E_error
import time
import barrett.util
import disc.gl_physics as phys

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
    # Loading data from output.txt, which consists of many rows with 18 columns
    #  (for tun with gp.ndim = 15)  (last of the 18 looks a bit odd)
    # The output.txt params are below run through gl_loglike to extract
    #  the physical parameters of each profile. To be saved in ..MN_out...save

    print ('In gl_m_h, m_out_to_prof')
    print ('basepath_ts:',basepath_ts)
    print ('mn_output_filename:',mn_output_filename)
    print ('gp.ndim:',gp.ndim)
    print ('gp.files.outdir:',gp.files.outdir,' profile_filename:',profile_filename)

    #paracube_data = np.loadtxt(basepath_ts + paracube_filename, dtype='a')

    gp.plotting_flag = True

    print ('len(mn_output_data[:,0])=',len(mn_output_data[:,0]))
    counting = 0
    for iter in range(0, len(mn_output_data[:,0])):
        #print('Iter = ', iter)
        if counting == 100:
            print ('iter:',iter)
            counting = 0
        counting = counting+1

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
        # Exporting the profiles to phys_MNoutput_profiles.save


def mn_output_to_hdf5(basepath_ts, mn_output_filename, h5_filename, investigation, case, timestamp, gp):
    #Define headers from the info file
    headers = ['mult','-2lnL'] + gp.param_headers
    units = []
    for jter in range (0, gp.ndim+2):
        units.append('')

    barrett.util.convert_chain([basepath_ts + mn_output_filename], headers, units, basepath_ts + h5_filename, 5000)
    return


def mn_h5_baryon_vecs(basepath_ts, mn_output_filename, h5_filename, investigation, case, timestamp, gp):
    #This appends the baryon density in each bin to the h5 file
    pdb.set_trace()
    bincenters, binmins, binmaxs, nudat, nuerr = gh.readcol5(gp.files.nufiles[0])

    h5file = barrett.data.Chain(basepath_ts + h5_filename)
    z_vec = bincenters
    for j in range(0, len(z_vec)):
        h5file.apply('rho_baryon_%j', '', lambda k, d: phys.rho_baryon_simplenu(z_vec[j], [k, d]), '$K_{\\rm DD}$', '$D_{\\rm DD}$')





def write_mn_info(gp):
    #Writes info file for multinest post analysis
    # H Silverwood 19/10/2015
    darcoda_path = gh.find_darcoda_path()
    ts = gp.files.timestamp

    info_file_name = darcoda_path +'/gravimage/DT'+ gp.investigate +'/'+ str(gp.case) +'/'+ ts +'/'+ 'output.info'
    info_file = open(info_file_name, 'w+')

    param_headers=[]

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

    #Dark matter parameters

    info_file.writelines('# Dark Matter Parameters \n')
    if gp.darkmattermodel == 'const_dm':
        info_file.writelines('#     DM model: const_dm \n')
        #info_file.writelines('lab' + str(lc) + '= \\rho_{\\rm DM, const} \n')
        #param_headers.append('\\rho_{\\rm DM, const}')
        info_file.writelines('lab' + str(lc) + '= \\rho_{\\rm DM, const} \n')
        param_headers.append('$\\rho_{\\rm DM, const}$')
        lc+=1
    elif gp.darkmattermodel == 'ConstPlusDD':
        info_file.writelines('#     DM model: ConstPlusDD \n')
        info_file.writelines('lab' + str(lc) + '= \\rho_{\\rm DM, const} \n')
        param_headers.append('$\\rho_{\\rm DM, const}$')
        lc+=1
        info_file.writelines('lab' + str(lc) + '= K_{\\rm DD} \n')
        param_headers.append('$K_{\\rm DD}$')
        lc+=1
        info_file.writelines('lab' + str(lc) + '= D_{\\rm DD} \n')
        param_headers.append('$D_{\\rm DD}$')
        lc+=1
    elif gp.darkmattermodel == 'gaussian_per_bin':
        info_file.writelines('lab' + str(lc) + '= \\rho_{\\rm DM, C} \n')
        param_headers.append('$\\rho_{\\rm DM, C}$')
        lc+=1
        for jter in range(0, gp.nrhonu-1):
            info_file.writelines('lab' + str(lc) + '= \\rho_{\\rm DM, ' + str(jter) +'} \n')
            param_headers.append('$\\rho_{\\rm DM, ' + str(jter) +'}$')
            lc+=1

    elif gp.darkmattermodel == 'kz_dm':
        info_file.writelines('#     DM model: kz_dm \n')
        info_file.writelines('lab' + str(lc) + '= k_{z,C} \n')
        param_headers.append('$k_{z,C}$')
        lc+=1
        for jter in range(0, gp.nrhonu):
            info_file.writelines('lab' + str(lc) + '= k_{z,' + str(jter) + '} \n')
            param_headers.append('$k_{z,' + str(jter) + '}$')
            lc+=1

    #Baryon mass profile parameters
    info_file.writelines('# Baryon Model Parameters \n')
    for baryon_pop in range(0, gp.nbaryon_pops):
        if gp.baryonmodel == 'simplenu_baryon':
            info_file.writelines('#     Baryon model: simplenu \n')
            info_file.writelines('lab' + str(lc) + '= K_{\\rm baryon} \n')
            param_headers.append('$K_{\\rm baryon}$')
            lc+=1
            info_file.writelines('lab' + str(lc) + '= D_{\\rm baryon} \n')
            param_headers.append('$D_{\\rm baryon}$')
            lc+=1

        elif gp.baryonmodel == 'trivial_baryon':
            info_file.writelines('#     Baryon model: trivial \n')
            info_file.writelines('lab' + str(lc) + '= \\Sigma_{\\rm baryon} \n')
            param_headers.append('$\\Sigma_{\\rm baryon}$') 
            lc+=1

        elif gp.baryonmodel == 'obs_baryon':
            info_file.writelines('#     Baryon model: obs \n')
            info_file.writelines('lab' + str(lc) + '= \\Sigma_{\\rm baryon} \n')
            param_headers.append('$\\Sigma_{\\rm baryon}$') 
            lc+=1
            info_file.writelines('lab' + str(lc) + '= \\Sigma_{\\rm dwarf} \n')
            param_headers.append('$\\Sigma_{\\rm dwarf}$')
            lc+=1
            info_file.writelines('lab' + str(lc) + '= \\rho 0_{\\rm MS} \n')
            param_headers.append('$\\rho 0_{\\rm MS}$')
            lc+=1
            info_file.writelines('lab' + str(lc) + '= h_{\\rm dwarf} \n')
            param_headers.append('$h_{\\rm dwarf}$')
            lc+=1
            info_file.writelines('lab' + str(lc) + '= h_{\\rm 1} \n')
            param_headers.append('$h_{\\rm 1}$')
            lc+=1
            info_file.writelines('lab' + str(lc) + '= h_{\\rm 2} \n')
            param_headers.append('$h_{\\rm 2}$')
            lc+=1
            info_file.writelines('lab' + str(lc) + '= h_{\\rm 3} \n')
            param_headers.append('$h_{\\rm 3}$')
            lc+=1
            info_file.writelines('lab' + str(lc) + '= x_{\\rm thick} \n')
            param_headers.append('$x_{\\rm thick}$')
            lc+=1
            info_file.writelines('lab' + str(lc) + '= \\Sigma_{\\rm HII} \n')
            param_headers.append('$\\Sigma_{\\rm HII}$')
            lc+=1
            info_file.writelines('lab' + str(lc) + '= h_{\\rm HII} \n')
            param_headers.append('$h_{\\rm HII}$')
            lc+=1

        elif gp.baryonmodel == 'simplenu_baryon_gaussian':  # Not supported
            info_file.writelines('# Baryon model: simplenu gaussian \n')
            info_file.writelines('lab' + str(lc) + '= rho_{{\rm baryon},C} \n')
            param_headers.append('$\\rho_{{\\rm baryon},C}$')
            lc+=1
            for jter in range(0, sum(gp.nbins)):
                info_file.writelines('lab' + str(lc) + '= rho_{{\rm baryon},' + str(jter) + '} \n')
                param_headers.append('$\\rho_{{\\rm baryon},' + str(jter) + '}$')
                lc+=1

    #Tracer profile parameters: nu_C, kz_nu_C, kz_nu_vector # kz_nu_LS
    info_file.writelines('# Tracer density model parameters \n')
    for tracer_pop in range(0, gp.ntracer_pops):
        if gp.scan_rhonu_space:
            info_file.writelines('#     Tracer model pop ' + str(tracer_pop) + ': scan_rhonu_space \n')
            raise Exception('Still working on this bit')

        elif gp.nu_model=='kz_nu':  # not supported
            info_file.writelines('#     Tracer model pop ' + str(tracer_pop) + ': kz_nu \n')
            info_file.writelines('lab' + str(lc) + '= \\nu_C \n')
            param_headers.append('$\\nu_C$')
            lc+=1
            #info_file.writelines('lab' + str(lc) + '= k_{\\nu, {\\rm C}} \n')
            #param_headers.append('k_{\\nu, {\\rm C}}')
            info_file.writelines('lab' + str(lc) + '= k_{\\nu,} \n')
            param_headers.append('$k_{\\nu,C}$')
            lc+=1
            for jter in range(0, gp.nbins[tracer_pop]):
                info_file.writelines('lab' + str(lc) + '= k_{\\nu,' + str(jter) +'} \n')
                param_headers.append('$k_{\\nu,' + str(jter) +'}$')
                lc+=1

        elif gp.nu_model=='gaussian_data': # not supported
            info_file.writelines('#     Tracer model pop ' + str(tracer_pop) + ': gaussian_data \n')
            info_file.writelines('lab' + str(lc) + '= \\nu_{\\rm C} \n')
            param_headers.append('$\\nu_{\\rm C}$')
            lc+=1
            for jter in range(0, gp.nbins[tracer_pop]):
                info_file.writelines('lab' + str(lc) + '= $\\nu_{' + str(jter) +'}$ \n')
                param_headers.append('$\\nu_{' + str(jter) + '}$')
                lc+=1


        elif gp.nu_model=='exponential_sum':
            info_file.writelines('#     Tracer model pop ' + str(tracer_pop) + ': sum of exponentials \n')
            for kter in range(0, gp.N_nu_model_exps):
                info_file.writelines('lab' + str(lc) + '= \\nu_{\\rm C, pop' + str(tracer_pop) + ', exp ' + str(kter) + '} \n')
                lc+=1
                info_file.writelines('lab' + str(lc) + '= \\nu_{\\rm h, pop' + str(tracer_pop) + ', exp ' + str(kter) + '} \n')
                lc+=1
                param_headers.append('$\\nu_{\\rm C, pop' + str(tracer_pop) + ', exp ' + str(kter) + '}$')
                param_headers.append('$\\nu_{\\rm h, pop' + str(tracer_pop) + ', exp ' + str(kter) + '}$')



    # Introducing tilt term:
    if gp.tilt:
        info_file.writelines('# Tilt term model parameters \n')
        for tracer_pop in range(0, gp.ntracer_pops):
            info_file.writelines('lab' + str(lc) + '= A_{' + str(tracer_pop) + ',{\\rm tilt}} \n')
            param_headers.append('$A_{' + str(tracer_pop) + ',{\\rm tilt}}$')
            lc+=1
            info_file.writelines('lab' + str(lc) + '= n_{' + str(tracer_pop) + ',{\\rm tilt}} \n')
            param_headers.append('$n_{' + str(tracer_pop) + ',{\\rm tilt}}$')
            lc+=1
            info_file.writelines('lab' + str(lc) + '= R_{' + str(tracer_pop) + ',{\\rm tilt}} \n')
            param_headers.append('$R_{' + str(tracer_pop) + ',{\\rm tilt}}$')
            lc+=1

    if gp.analytic_sigz2 == False:
        info_file.writelines('# sigc constant C parameter(s) \n')
        for t_pop in range(0, gp.ntracer_pops):
            info_file.writelines('lab' + str(lc) + '= C_' + str(t_pop) + '\n')
            param_headers.append('$C_' + str(t_pop) + '$')
            lc+=1

    print ('lc-1:',lc-1,' gp.ndim:',gp.ndim) 
    if lc-1 != gp.ndim:
        raise Exception('Incorrect number of labels in mn_info')

    info_file.writelines('# Number of Parameters saved \n')
    info_file.writelines('params_saved = ' + str(gp.ndim))
    info_file.close()
    gp.param_headers = param_headers
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
