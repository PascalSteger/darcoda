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

def write_disc_output_files(Bincenter, Binmin, Binmax, nudat, nuerr, sigz2dat, sigz2err, gp):

    # write tracer densities 3D
    file_nu = open(gp.files.nufiles[0], 'w')
    print('Bin centres z [kpc];','Binmin z [kpc];','Binmax z [kpc];','nu(z) [#/pc^3];','error [#/pc^3]', file=file_nu)
    for b in range(gp.nbins):
        print(Bincenter[b], Binmin[b], Binmax[b], nudat[b], nuerr[b], file=file_nu)
    file_nu.close()

    file_sig = open(gp.files.sigfiles[0],'w')
    print('Bin centres z [kpc];','Binmin z [kpc];','Binmax [kpc];', 'sigma_z^2(z) [km^2/s^2];', 'error [km^2/s^2]', file=file_sig)
    for b in range(gp.nbins):
        print(Bincenter[b], Binmin[b], Binmax[b], sigz2dat[b], sigz2err[b], file=file_sig)
    file_sig.close()


def ErSamp_gauss_linear_w_z():
    fraction_err = 0.05
    datafile = '/home/hsilverw/LoDaM/darcoda/Data_Sets/simplenu/simplenu_sigz_raw_sdz_p05_sdvz_5.dat' #Todo: work out how to un-hard code this
    data = np.loadtxt(datafile)
    z_data = data[:, 0]

    z_sampled = []
    for z_val in z_data:
        z_sampled.append(rand.normal(loc = z_val, scale= z_val*fraction_err))

    return z_sampled


## \fn write_disc_output_files(Rbin, Binmin, Binmax, nudat, nuerr, Sigdat, Sigerr, Mdat, Merr, sigdat, sigerr, scales, gp)
# for permanent and consistent data handling
# @param Rbin
# @param Binmin
# @param Binmax
# @param nudat
# @param nuerr
# @param Sigdat
# @param Sigerr
# @param Mdat
# @param Merr
# @param sigdat
# @param sigerr
# @param scales
# @param gp global parameters


def run(gp):
    if gp.machine == 'lisa_HS_login' or gp.machine == 'lisa_HS_batch':
        external_file='/home/hsilverw/LoDaM/darcoda/Data_Sets/' + gp.external_data_file
    elif gp.machine == 'lisa_SS_login' or gp.machine == 'lisa_SS_batch':
        external_file='/home/sofia/darcoda/gravlite/Data_Sets/' + gp.external_data_file

    external_data = np.loadtxt(external_file)

    z_data = external_data[:, 0] #[kpc]
    v_data = external_data[:, 1] #[km/s]

    #Population Splitting
    if gp.ntracer_pops > 1:
        print('More than one population')
        #DO POPULATION SPLITTING

    # Bin data, constant number of tracers per bin.
    # First, calculate binmins, binmaxes, and bin centres
    if gp.binning == 'consttr':
        binmin, binmax, bincentermed = gh.bin_r_const_tracers(z_data, gp.nbins)
    elif gp.binning == 'linspace':
        binmin, binmax, bincentermed = gh.bin_r_linear(0., round(max(z_data),1), gp.nbins)

    # Then calculate tracer number density [#stars/kpc^3], [#stars/kpc^3], [km/s], [km/s]
    nu_data, nu_err_pois, sigz2_data, sigz2_err_pois, Ntr_per_bin = gh.nu_sig_from_bins(binmin, binmax, z_data, v_data)

    # Use MC to estimate errors on nu
    if gp.investigate == 'simplenu':
        z_sampler = ErSamp_gauss_linear_w_z
    nu_err_meas = gmcer.mc_nu_error(z_sampler, gp.mc_err_N_iters, binmin, binmax, bincentermed)

    #Combine nu errors in quadrature TODO check that this is the right way to add errors
    nu_err_tot = np.sqrt(nu_err_pois**2 + nu_err_meas**2)

    #Combine sig (vel disp) errors
    Ntr = len(z_data) #Number of tracers
    sigz2_err_tot = sigz2_err_pois + np.sqrt(2/Ntr_per_bin)*gp.vz_SDerr_meas**2

    #Output data to file
    write_disc_output_files(bincentermed, binmin, binmax, nu_data, nu_err_tot, sigz2_data, sigz2_err_tot, gp)

    #Set central nu prior range
    gp.nu_C_max = 2.*Ntr/(binmax[0]-binmin[0])
    gp.nu_C_min = 0.1*Ntr/(binmax[-1]-binmin[0])

    #Temp:
    gp.nu_C_max = 3.E4
    gp.nu_C_min = 1.E3

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
