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

#from gl_collection import ProfileCollection#.true_sigz2_func
from plotting.gl_collection import ProfileCollection

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
    velocity_err = 5. #km s^-1

    if gp.machine == 'lisa_HS_login' or gp.machine == 'lisa_HS_batch':
        datafile='/home/hsilverw/LoDaM/darcoda/Data_Sets/simplenu/simplenu_sigz_raw_sdz_p05_sdvz_5.dat' #Todo: work out how to un-hard code this
    elif gp.machine == 'lisa_SS_login' or gp.machine == 'lisa_SS_batch':
        datafile='/home/sofia/darcoda/Data_Sets/simplenu/simplenu_sigz_raw_sdz_p05_sdvz_5.dat'

    data = np.loadtxt(datafile)
    z_data = data[:, 0]
    v_data = data[:, 1]

    z_sampled = []
    v_sampled = []
    for z_val in z_data:
        z_sampled.append(rand.normal(loc = z_val, scale = z_val*fraction_err))
    for v_val in v_data:
        v_sampled.append(rand.normal(loc = v_val, scale = velocity_err))

    return [z_sampled, v_sampled]

# SS: seems to generate a new file with added measurement errors (or equivalently estimate a realization of true values given the
#     measurement with measurement error). Looks only at z_data (i.e. nu).
# HS: takes positions and velocities of starts, uses these as centres for Gaussians, then draws a
#    new set of positions and velocities from those Gaussians.


def run(gp):
    if gp.machine == 'lisa_HS_login' or gp.machine == 'lisa_HS_batch':
        external_file='/home/hsilverw/LoDaM/darcoda/Data_Sets/' + gp.external_data_file
    elif gp.machine == 'lisa_SS_login' or gp.machine == 'lisa_SS_batch':
        external_file='/home/sofia/darcoda/Data_Sets/' + gp.external_data_file

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

    if gp.TheoryData == True:
        G1 = 4.299e-6   # Newton's constant in (km)^2*kpc/(Msun*s^2)
        K=1500.
        F=267.65
        D=0.18
        z0=0.4
        normC = 22.85
        nu0 = 1.
        Ntr = 9516.
        zmax= 1.2

        nu0 = Ntr/(z0*(1-np.exp(-zmax/z0)))
        truen_arr = nu0*z0*(np.exp(-1.*binmin/z0)-np.exp(-1.*binmax/z0))  # true n.o. stars in bins
        truenu_arr = truen_arr/(binmax-binmin)
        nu_data = truenu_arr

        true_sigz2_arr = ProfileCollection(gp.ntracer_pops,gp.nbins).true_sigz2_func(binmin,binmax,1000,gp.nbins,z0,Ntr,nu0,K,D,F,gp)
        sigz2_data = true_sigz2_arr

    # Use MC to estimate errors on nu
    if gp.investigate == 'simplenu':
        z_sampler = ErSamp_gauss_linear_w_z
    [nu_err_meas, sigz2_err_meas] = gmcer.mc_nu_error(z_sampler, gp.mc_err_N_iters, binmin, binmax, bincentermed)

    #Combine nu errors in quadrature TODO check that this is the right way to add errors
    nu_err_tot = np.sqrt(nu_err_pois**2 + nu_err_meas**2)

    #Combine sig (vel disp) errors
    Ntr = len(z_data) #Number of tracers
    #sigz2_err_tot = sigz2_err_pois + np.sqrt(2/Ntr_per_bin)*gp.vz_SDerr_meas**2
    sigz2_err_tot = np.sqrt(sigz2_err_pois**2 + sigz2_err_meas**2)

    #SS: Running the code without measurement errors   TODO FIXME !!
    nu_err_tot = nu_err_pois     # SS-TODO
    sigz2_err_tot = sigz2_err_pois  # SS-TODO

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
