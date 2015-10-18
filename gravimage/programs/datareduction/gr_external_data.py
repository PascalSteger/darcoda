#!/usr/bin/env ipython3

##
# @ file
# generate simple disc data ourselves

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch


import numpy as np
import numpy.random as npr
import pdb
from scipy.integrate import simps

import gl_units as gu
import gl_helper as gh
import gl_mc_errors as gmcer
import numpy.random as rand
import os

#from gl_collection import ProfileCollection#.true_sigz2_func
from plotting.gl_collection import ProfileCollection

def write_disc_output_files(Bincenter, Binmin, Binmax, nudat, nuerr, sigz2dat, sigz2err, gp):
    # write tracer densities 3D
    for pop in range(0, gp.ntracer_pops):
        file_nu = open(gp.files.nufiles[pop], 'w')
        print('Bin centres z [kpc];','Binmin z [kpc];','Binmax z [kpc];','nu(z) [#/pc^3];','error [#/pc^3]', file=file_nu)
        for b in range(gp.nbins[pop]):
            print(Bincenter[pop][b], Binmin[pop][b], Binmax[pop][b], nudat[pop][b], nuerr[pop][b], file=file_nu)
        file_nu.close()
        print ('gp.files.nufiles[pop]:',gp.files.nufiles[pop])

        file_sig = open(gp.files.sigfiles[pop],'w')
        print('Bin centres z [kpc];','Binmin z [kpc];','Binmax [kpc];', 'sigma_z^2(z) [km^2/s^2];', 'error [km^2/s^2]', file=file_sig)
        for b in range(gp.nbins[pop]):
            print(Bincenter[pop][b], Binmin[pop][b], Binmax[pop][b], sigz2dat[pop][b], sigz2err[pop][b], file=file_sig)
        file_sig.close()
        print ('gp.files.sigfiles[pop]:',gp.files.sigfiles[pop])

def write_tilt_output_files(Bincenter, Binmin, Binmax, tilt, tilterr, gp):
    for pop in range(0, gp.ntracer_pops):
        file_tilt = open(gp.files.tiltfiles[pop], 'w')
        print('Bin centres z [kpc];','Binmin z [kpc];','Binmax z [kpc];','sigmaRz(z) [km^2/s^2];','error [km^2/s^2]', file=file_tilt)
        for b in range(gp.nbins[pop]):
            print(Bincenter[pop][b], Binmin[pop][b], Binmax[pop][b], tilt[pop][b], tilterr[pop][b], file=file_tilt)
        file_tilt.close()
        print ('gp.files.tiltfiles[pop]:',gp.files.tiltfiles[pop])

def write_sigRz2_output_files(Bincenter, Binmin, Binmax, sigRz2_dat, sigRz2_err, gp):
    for pop in range(0, gp.ntracer_pops):
        file_sigRz2 = open(gp.files.sigRz2_files[pop], 'w')
        print('Bin centres z [kpc];','Binmin z [kpc];','Binmax z [kpc];','sigmaRz(z) [km^2/s^2];','error [km^2/s^2]', file=file_sigRz2)
        for b in range(gp.nbins[pop]):
            print(Bincenter[pop][b], Binmin[pop][b], Binmax[pop][b], sigRz2_dat[pop][b], sigRz2_err[pop][b], file=file_sigRz2)
        file_sigRz2.close()
        print ('gp.files.sigRz2_files[pop]:',gp.files.sigRz2_files[pop])


def load_simplenu_posvel(gp):
    darcoda_path = gh.find_darcoda_path()
    data_set_folder = darcoda_path + '/Data_Sets/'
    external_file = [data_set_folder + temp for temp in gp.external_data_file]
    print('External file name = ', external_file)
    if gp.tilt:
        external_file_tilt=[data_set_folder + temp for temp in gp.external_data_file_tilt]
        print('External tilt file name = ', external_file_tilt)

    #Load external data
    external_data = [np.loadtxt(file_name) for file_name in external_file]
    if gp.tilt:
        external_data_tilt = [np.loadtxt(file_name_tilt) for file_name_tilt in external_file_tilt]

    #Extract data
    z_data = [external_data[ii][:, 0] for ii in range(0, len(external_data))] #[kpc]
    v_data = [external_data[ii][:, 1] for ii in range(0, len(external_data))] #[km/s]
    if gp.tilt:
        z_data_tilt = [external_data_tilt[ii][:,0] for ii in range(0, len(external_data))]     #[kpc]
        vRz_data_tilt = [external_data_tilt[ii][:,1] for ii in range(0, len(external_data))]   #[km^2/s^2]
    else:
        z_data_tilt = 0.
        vRz_data_tilt = 0.

    # z-cut, use only stars below a certain height for each population
    z_data_used = [z_data[ii][z_data[ii] <gp.data_z_cut[ii]] for ii in range(0, len(external_data))]
    v_data_used = [v_data[ii][z_data[ii] <gp.data_z_cut[ii]] for ii in range(0, len(external_data))]

    return z_data_used, v_data_used, z_data_tilt, vRz_data_tilt



def load_disc_nbody_posvel(gp):
    darcoda_path = gh.find_darcoda_path()
    data_set_folder = darcoda_path + '/Data_Sets/'
    external_file = [data_set_folder + temp for temp in gp.external_data_file]
    if gp.tilt:
        external_file_tilt=[data_set_folder + temp for temp in gp.external_data_file_tilt]

    #Load external data
    print('Using external data file(s): ', gp.external_data_file)
    external_data = [np.loadtxt(file_name) for file_name in external_file]
    z_data = [external_data[ii][:, 0] for ii in range(0, len(external_data))] #[kpc]
    vz_data = [external_data[ii][:, 1] for ii in range(0, len(external_data))] #[km/s]
    vR_data = [external_data[ii][:, 2] for ii in range(0, len(external_data))] #[km/s]

    #Centre in z
    z_data = z_data - np.mean(z_data) #mean for all populations
    #Remove peculiar velocity in z
    vz_data = vz_data - np.mean(vz_data)

    #Mirror negative data upwards - assumes symmetry in the plane
    negative_z_mask = np.sign(z_data) #-1 for z<0, 1 for z>0
    z_data = z_data*negative_z_mask
    vz_data = vz_data*negative_z_mask

    #Cut on z
    z_data_cut = [z_data[ii][z_data[ii] <gp.data_z_cut[ii]] for ii in range(0, len(external_data))]
    vz_data_cut = [vz_data[ii][z_data[ii] <gp.data_z_cut[ii]] for ii in range(0, len(external_data))]
    vR_data_cut = [vR_data[ii][z_data[ii] <gp.data_z_cut[ii]] for ii in range(0, len(external_data))]

    #Calculate vzvR
    vzvR_data_cut = [vz_data_cut[jj]*vR_data_cut[jj] for jj in range(0, len(vz_data_cut))]

    return z_data_cut, vz_data_cut, vR_data_cut


def run(gp):

    if gp.investigate == 'simplenu':
        z_data_used, vz_data_used, z_data_tilt, vRz_data_tilt = load_simplenu_posvel(gp)
    elif gp.investigate == 'disc_nbody':
        z_data_used, vz_data_used, vR_data_used = load_disc_nbody_posvel(gp)

    #vz_data_used = np.array(vz_data_used)
    #vR_data_used = np.array(vR_data_used)

    # Bin data, constant number of tracers per bin.
    # First, calculate binmins, binmaxes, and bin centres
    binmin_pops = []
    binmax_pops = []
    bincentermed_pops = []

    for pop in range(0, gp.ntracer_pops):
        if gp.binning == 'consttr':
            binmin, binmax, bincentermed = gh.bin_r_const_tracers(z_data_used[pop], gp.nbins[pop], np.ones(len(z_data_used[pop])))
            #pdb.set_trace()
            #binmin, binmax, bincentermed = gh.bin_r_const_tracers_weighted(z_data_used[pop], gp.nbins[pop], np.ones(len(z_data_used[pop])))
        elif gp.binning == 'linspace':
            binmin, binmax, bincentermed = gh.bin_r_linear(0., round(max(z_data_used[0]),1), gp.nbins[pop])

        binmin_pops.append(binmin)
        binmax_pops.append(binmax)
        bincentermed_pops.append(bincentermed)

    # Then calculate tracer number density [#stars/kpc^3], [#stars/kpc^3], [km/s], [km/s]
    nu_data=[]
    nu_err_pois=[]
    sigz2_data = []
    sigz2_err_pois = []
    sigRz2_data = []
    sigRz2_err_pois = []
    Ntr_per_bin = []
    tilt2_data = []
    tilt2_data_err = []

    for pop in range(0, gp.ntracer_pops):
        nu_data_tmp, nu_err_pois_tmp, sigz2_data_tmp, sigz2_err_pois_tmp, Ntr_per_bin_tmp = gh.nu_sig_from_bins(binmin_pops[pop], binmax_pops[pop], z_data_used[pop], vz_data_used[pop], np.ones(len(z_data_used[pop])))
        nu_data.append(nu_data_tmp)
        nu_err_pois.append(nu_err_pois_tmp)
        sigz2_data.append(sigz2_data_tmp)
        sigz2_err_pois.append(sigz2_err_pois_tmp)
        Ntr_per_bin.append(Ntr_per_bin_tmp)

        if gp.tilt:
            if gp.investigate == 'simplenu':
                sigRz2_data_tmp, sigRz2_data_err_tmp = gh.sigRz_from_bins_simplenu(binmin_pops[pop], binmax_pops[pop], z_data_tilt[pop], vRz_data_tilt[pop])
            elif gp.investigate == 'disc_nbody':
                sigRz2_data_tmp, sigRz2_data_err_tmp = gh.sigRz_from_bins(binmin_pops[pop], binmax_pops[pop], z_data_used[pop], vz_data_used[pop], vR_data_used[pop])

            sigRz2_data.append(sigRz2_data_tmp)
            sigRz2_err_pois.append(sigRz2_data_err_tmp)



    print ('z vectors: ', bincentermed_pops)
    print ('nu_data:',nu_data)
    print ('nu_err_pois:',nu_err_pois)
    print ('sigz2_data:',sigz2_data)
    print ('sigz2_err_pois:',sigz2_err_pois)
    if gp.tilt:
        print('sigRz2_data = ', sigRz2_data)
        print('sigRz2_err_pois = ', sigRz2_err_pois)


    Ntr_used = np.array([len(z_data_used[ii]) for ii in range(0, len(z_data_used))]) # Total number of tracers used (& binned)

    ## Use MC to estimate errors on nu
    #if gp.investigate == 'simplenu':
    #    z_sampler = ErSamp_gauss_linear_w_z
    #[nu_err_meas, sigz2_err_meas] = gmcer.mc_nu_error(z_sampler, gp.mc_err_N_iters, binmin, binmax, bincentermed)

    ##Combine nu errors in quadrature TODO check that this is the right way to add errors
    #nu_err_tot = np.sqrt(nu_err_pois**2 + nu_err_meas**2)

    #Combine sig (vel disp) errors
    #Ntr = np.array([len(z_data[ii]) for ii in range(0,len(z_data))]) #Total number of tracers avaliable

    #sigz2_err_tot = sigz2_err_pois + np.sqrt(2/Ntr_per_bin)*gp.vz_SDerr_meas**2
    #sigz2_err_tot = np.sqrt(sigz2_err_pois**2 + sigz2_err_meas**2)

    #SS: Running the code without measurement errors   TODO FIXME !!
    nu_err_tot = nu_err_pois     # SS-TODO
    sigz2_err_tot = sigz2_err_pois  # SS-TODO
    sigRz2_err_tot = sigRz2_err_pois

    #Output data to file
    #TO DO MULTIPOPS BODGE

    write_disc_output_files(bincentermed_pops, binmin_pops, binmax_pops, nu_data, nu_err_tot, sigz2_data, sigz2_err_tot, gp)
    if gp.tilt:
        #write_tilt_output_files(bincentermed_pops, binmin_pops, binmax_pops, tilt2_data, tilt2_data_err, gp)
        write_sigRz2_output_files(bincentermed_pops, binmin_pops, binmax_pops, sigRz2_data, sigRz2_err_tot, gp)

    #Set central nu and sigz prior range #TODO: look again at this prior
    nu_0_values = [nu_data[ii][0] for ii in range(0, gp.ntracer_pops)]
    nu_0_errs = [nu_err_tot[ii][0] for ii in range(0, gp.ntracer_pops)]
    sigz_0_values = np.sqrt([sigz2_data[ii][0] for ii in range(0, gp.ntracer_pops)])
    sigz_0_errs = np.sqrt([sigz2_err_tot[ii][0] for ii in range(0, gp.ntracer_pops)])

    #gp.nu_C_max = 10*max(nu_0_values)
    #gp.nu_C_min = 0.1*min(nu_0_values)

    gp.nu_C_max = max(nu_0_values) + 5*max(nu_0_errs) # 5 sigma either way
    gp.nu_C_min = min(nu_0_values) - 5*max(nu_0_errs)

    #gp.sigz_C_max = 10*max(sigz_0_values) #set from data in gr_external_data
    #gp.sigz_C_min = 0.1*min(sigz_0_values)

    gp.sigz_C_max = max(sigz_0_values) + 5*max(sigz_0_errs)
    gp.sigz_C_min = max(min(sigz_0_values) - 5*max(sigz_0_errs), 0.)

    import gr_params #WHAT DOES THIS DO.
    gpr = gr_params.Params(gp)
    if gpr.showplots:
        nuscaleb = nu_zth[np.argmin(np.abs(zth-z0))]
        plt.loglog(zth, nu_zth/nuscaleb, 'b.-')
        nuscaler = nu_dat_bin1[np.argmin(np.abs(zth-z0))]
        plt.loglog(zth, nu_dat_bin1/nuscaler, 'r.-')

    return gp.dat

## \fn run(gp)
# Read in data from JR's simple 1D mock
# Bin it, calculate nu and sigma, output it to file
# @param gp global parameters

if __name__=="__main__":
    import gl_params
    gp = gl_params.Params()
    run(gp)
