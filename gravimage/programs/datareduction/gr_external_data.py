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
    for pop in range(0, gp.ntracer_pops):
        file_nu = open(gp.files.nufiles[pop], 'w')
        print('Bin centres z [kpc];','Binmin z [kpc];','Binmax z [kpc];','nu(z) [#/pc^3];','error [#/pc^3]', file=file_nu)
        for b in range(gp.nbins):
            print(Bincenter[b], Binmin[b], Binmax[b], nudat[pop][b], nuerr[pop][b], file=file_nu)
        file_nu.close()
        print ('gp.files.nufiles[pop]:',gp.files.nufiles[pop])

        file_sig = open(gp.files.sigfiles[pop],'w')
        print('Bin centres z [kpc];','Binmin z [kpc];','Binmax [kpc];', 'sigma_z^2(z) [km^2/s^2];', 'error [km^2/s^2]', file=file_sig)
        for b in range(gp.nbins):
            print(Bincenter[b], Binmin[b], Binmax[b], sigz2dat[pop][b], sigz2err[pop][b], file=file_sig)
        file_sig.close()
        print ('gp.files.sigfiles[pop]:',gp.files.sigfiles[pop])

def write_tilt_output_files(Bincenter, Binmin, Binmax, tilt, tilterr, gp):
    for pop in range(0, gp.ntracer_pops):
        file_tilt = open(gp.files.tiltfiles[pop], 'w')
        print('Bin centres z [kpc];','Binmin z [kpc];','Binmax z [kpc];','sigmaRz(z) [km^2/s^2];','error [km^2/s^2]', file=file_tilt)
        for b in range(gp.nbins):
            print(Bincenter[b], Binmin[b], Binmax[b], tilt[pop][b], tilterr[pop][b], file=file_tilt)
        file_tilt.close()
        print ('gp.files.tiltfiles[pop]:',gp.files.tiltfiles[pop])

def ErSamp_gauss_linear_w_z():
    fraction_err = 0.05
    velocity_err = 5. #km s^-1

    #datafile='/home/hsilverw/LoDaM/darcoda/Data_Sets/simplenu/simplenu_sigz_raw_sdz_p05_sdvz_5.dat'
    datafile='/home/sofia/darcoda/Data_Sets/simplenu/simplenu_sigz_raw_sdz_p05_sdvz_5.dat'
    #Todo: work out how to un-hard code this. Ultimately this should be a method of a data class.

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
# HS: takes positions and velocities of stars, uses these as centres for Gaussians, then draws a
#    new set of positions and velocities from those Gaussians.


def run(gp):
    if gp.machine == 'lisa_HS_login' or gp.machine == 'lisa_HS_batch':
        external_file=['/home/hsilverw/LoDaM/darcoda/Data_Sets/' + temp for temp in gp.external_data_file]
        if gp.tilt:
            external_file_tilt=['/home/hsilverw/LoDaM/darcoda/Data_Sets/' + temp for temp in gp.external_data_file_tilt]
    elif gp.machine == 'lisa_SS_login' or gp.machine == 'lisa_SS_batch':
        external_file=['/home/sofia/darcoda/Data_Sets/' + temp for temp in gp.external_data_file]
        if gp.tilt:
            external_file_tilt=['/home/sofia/darcoda/Data_Sets/' + temp for temp in gp.external_data_file_tilt]

    #Check number of data files corresponds with population count
    if len(external_file) != gp.ntracer_pops:
        raise Exception('Incorect number of data files for number of tracer populations')
    if gp.tilt:
        if len(external_file_tilt) != gp.ntracer_pops:
            raise Exception('Incorect number of tilt data files for number of tracer populations')

    external_data = [np.loadtxt(file_name) for file_name in external_file]
    if gp.tilt:
        external_data_tilt = [np.loadtxt(file_name_tilt) for file_name_tilt in external_file_tilt]

    z_data = [external_data[ii][:, 0] for ii in range(0, len(external_data))] #[kpc]
    v_data = [external_data[ii][:, 1] for ii in range(0, len(external_data))] #[km/s]

    if gp.tilt:
        z_data_tilt = [external_data_tilt[ii][:,0] for ii in range(0, len(external_data))]     #[kpc]
        vRz_data_tilt = [external_data_tilt[ii][:,1] for ii in range(0, len(external_data))]   #[km^2/s^2]

    # z-cut, use only stars below a certain height for each population
    z_data_used = [z_data[ii][z_data[ii]<gp.data_z_cut[ii]] for ii in range(0, len(external_data))]
    v_data_used = [v_data[ii][z_data[ii]<gp.data_z_cut[ii]] for ii in range(0, len(external_data))]

    [print ('N.o. tracer stars used as data for popn ', t_pop, ': ', len(z_data_used[t_pop])) for t_pop in range(0, len(external_data))]

    # Bin data, constant number of tracers per bin.
    # First, calculate binmins, binmaxes, and bin centres
    binmin_pops = []
    binmax_pops = []
    bincentermed_pops = []

    for pop in range(0, gp.ntracer_pops):
        if gp.binning == 'consttr':
            binmin, binmax, bincentermed = gh.bin_r_const_tracers(z_data_used[pop], gp.nbins)
        elif gp.binning == 'linspace':
            binmin, binmax, bincentermed = gh.bin_r_linear(0., round(max(z_data_used[0]),1), gp.nbins)

        binmin_pops.append(binmin)
        binmax_pops.append(binmax)
        bincentermed_pops.append(bincentermed)


    # Then calculate tracer number density [#stars/kpc^3], [#stars/kpc^3], [km/s], [km/s]
    nu_data=[]
    nu_err_pois=[]
    sigz2_data = []
    sigz2_err_pois = []
    Ntr_per_bin = []
    tilt2_data = []
    tilt2_data_err = []
    for pop in range(0, gp.ntracer_pops):
        nu_data_tmp, nu_err_pois_tmp, sigz2_data_tmp, sigz2_err_pois_tmp, Ntr_per_bin_tmp = gh.nu_sig_from_bins(binmin_pops[pop], binmax_pops[pop], z_data[pop], v_data[pop])
        nu_data.append(nu_data_tmp)
        nu_err_pois.append(nu_err_pois_tmp)
        sigz2_data.append(sigz2_data_tmp)
        sigz2_err_pois.append(sigz2_err_pois_tmp)
        Ntr_per_bin.append(Ntr_per_bin_tmp)

        if gp.tilt:
            tilt2_data_tmp, tilt2_data_err_tmp = gh.tilt_from_bins(binmin_pops[pop], binmax_pops[pop], z_data_tilt[pop], vRz_data_tilt[pop])
            tilt2_data.append(tilt2_data_tmp)
            tilt2_data_err.append(tilt2_data_err_tmp)
    #MULTIPOPS REWRITE WORKING POINT

    #nu_data, nu_err_pois, sigz2_data, sigz2_err_pois, Ntr_per_bin = gh.nu_sig_from_bins(binmin, binmax, z_data, v_data) # faster to instead use .._data_used

    if gp.tilt:
        print ('tilt2_data:',tilt2_data)
        print ('tilt2_data_fit:',(tilt2_data-np.square(0.0087*(1000.*bincentermed)**1.44)/tilt2_data))
    #pdb.set_trace()

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
        #nu_data = truenu_arr

        true_sigz2_arr = ProfileCollection(gp.ntracer_pops,gp.nbins).true_sigz2_func(binmin,binmax,1000,gp.nbins,z0,Ntr,nu0,K,D,F,gp)

        print ('truenu_arr=',truenu_arr)
        print ('true_sigz2_arr=',true_sigz2_arr)

        true_chi2_nu = np.sum((truenu_arr-nu_data)**2./nu_err_pois**2.)
        true_chi2_sigz2 = np.sum((true_sigz2_arr-sigz2_data)**2./sigz2_err_pois**2.)
        print ('Without measurement errors the true result has the Chi2:')
        print ('nu:',true_chi2_nu/20.,'sigz2:',true_chi2_sigz2/20.,'average:',(true_chi2_nu+true_chi2_sigz2)/40.)

        nu_data = truenu_arr
        sigz2_data = true_sigz2_arr

    print ('nu_data:',nu_data)
    print ('sigz2_data:',sigz2_data)
    print ('nu_err_pois:',nu_err_pois)
    print ('sigz2_err_pois:',sigz2_err_pois)

    ## Use MC to estimate errors on nu  # SS: Throws error if file not set correctly, so commented out for now as meas errors are currently not used.
    #if gp.investigate == 'simplenu':
    #    z_sampler = ErSamp_gauss_linear_w_z
    #[nu_err_meas, sigz2_err_meas] = gmcer.mc_nu_error(z_sampler, gp.mc_err_N_iters, binmin, binmax, bincentermed)

    ##Combine nu errors in quadrature TODO check that this is the right way to add errors
    #nu_err_tot = np.sqrt(nu_err_pois**2 + nu_err_meas**2)

    #Combine sig (vel disp) errors
    Ntr = np.array([len(z_data[ii]) for ii in range(0,len(z_data))]) #Total number of tracers avaliable
    Ntr_used = np.array([len(z_data_used[ii]) for ii in range(0, len(z_data))]) # Total number of tracers used (& binned)
    #sigz2_err_tot = sigz2_err_pois + np.sqrt(2/Ntr_per_bin)*gp.vz_SDerr_meas**2
    #sigz2_err_tot = np.sqrt(sigz2_err_pois**2 + sigz2_err_meas**2)

    #SS: Running the code without measurement errors   TODO FIXME !!
    nu_err_tot = nu_err_pois     # SS-TODO
    sigz2_err_tot = sigz2_err_pois  # SS-TODO

    #Output data to file
    #TO DO MULTIPOPS BODGE

    write_disc_output_files(bincentermed_pops[0], binmin_pops[0], binmax_pops[0], nu_data, nu_err_tot, sigz2_data, sigz2_err_tot, gp)
    if gp.tilt:
        write_tilt_output_files(bincentermed, binmin, binmax, tilt2_data, tilt2_data_err, gp)

    #Set central nu prior range #TODO: look again at this prior
    gp.nu_C_max = max(2.*Ntr_used/(binmax[0]-binmin[0]))
    gp.nu_C_min = min(0.1*Ntr_used/(binmax[-1]-binmin[0]))

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
