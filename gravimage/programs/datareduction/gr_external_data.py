#!/usr/bin/env ipython3

##
# @ file
# generate simple disc data ourselves

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch


import numpy as np
import numpy.random as npr
import pdb
from scipy.integrate import simps
from scipy.signal import savgol_filter

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

    return z_data, v_data, z_data_tilt, vRz_data_tilt



def load_disc_nbody_posvel(gp):
    darcoda_path = gh.find_darcoda_path()
    data_set_folder = darcoda_path + '/Data_Sets/'
    external_file = [data_set_folder + temp for temp in gp.external_data_file]
    #if gp.tilt:
    #    external_file_tilt=[data_set_folder + temp for temp in gp.external_data_file_tilt]

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

    return z_data, vz_data, vR_data

def run(gp):
    if gp.investigate == 'simplenu':
        z_data, vz_data, z_data_tilt, vRz_data_tilt = load_simplenu_posvel(gp)
        z_data_zcut = [z_data[ii][z_data[ii] <gp.data_z_cut[ii]] for ii in range(0, len(z_data))]
        #v_data_zcut = [vz_data[ii][z_data[ii] <gp.data_z_cut[ii]] for ii in range(0, len(z_data))]
    elif gp.investigate == 'disc_nbody':
        z_data, vz_data, vR_data = load_disc_nbody_posvel(gp)
        z_data_zcut = [z_data[ii][z_data[ii] <gp.data_z_cut[ii]] for ii in range(0, len(z_data))]
        #vz_data_zcut = [vz_data[ii][z_data[ii] <gp.data_z_cut[ii]] for ii in range(0, len(external_data))]
        #vR_data_zcut = [vR_data[ii][z_data[ii] <gp.data_z_cut[ii]] for ii in range(0, len(external_data))]

    # Bin data, constant number of tracers per bin.
    # First, calculate binmins, binmaxes, and bin centres
    binmin_pops = []
    binmax_pops = []
    bincentermed_pops = []

    for pop in range(0, gp.ntracer_pops):
        if gp.binning == 'consttr':
            binmin, binmax, bincentermed = gh.bin_r_const_tracers(z_data_zcut[pop], gp.nbins[pop],
                                                np.ones(len(z_data_zcut[pop])), gp.data_z_cut[pop])
        elif gp.binning == 'linspace':
            binmin, binmax, bincentermed = gh.bin_r_linear(0., round(max(z_data_zcut[0]),1), gp.nbins[pop])

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
        nu_data_tmp, nu_err_pois_tmp, sigz2_data_tmp, sigz2_err_pois_tmp, Ntr_per_bin_tmp\
            = gh.nu_sig_from_bins(binmin_pops[pop], binmax_pops[pop], z_data[pop],
                                    vz_data[pop], np.ones(len(z_data[pop])))

        nu_data.append(nu_data_tmp)
        nu_err_pois.append(nu_err_pois_tmp)
        sigz2_data.append(sigz2_data_tmp)
        sigz2_err_pois.append(sigz2_err_pois_tmp)
        Ntr_per_bin.append(Ntr_per_bin_tmp)

        if gp.tilt:
            if gp.investigate == 'simplenu':
                sigRz2_data_tmp, sigRz2_data_err_tmp = gh.sigRz2_from_bins_simplenu(binmin_pops[pop],
                                                            binmax_pops[pop], z_data_tilt[pop], vRz_data_tilt[pop])

                #Reverse sigRz2 sign if necessary
                if not gp.positive_sigRz_data_sign:
                    sigRz2_data_tmp = sigRz2_data_tmp * -1.0

            elif gp.investigate == 'disc_nbody':
                sigRz2_data_tmp, sigRz2_data_err_tmp = gh.sigRz2_from_bins(binmin_pops[pop],
                                                            binmax_pops[pop], z_data[pop], vz_data[pop], vR_data[pop])

            sigRz2_data.append(sigRz2_data_tmp)
            sigRz2_err_pois.append(sigRz2_data_err_tmp)

    ##Nu bias correction
    #dh_vec = (binmax-binmin)/0.9
    #correction_vec = 0.5 * dh_vec * (np.exp(dh_vec)+1)/(np.exp(dh_vec)-1)
    #print('nu correction:', correction_vec)
    ##nu_data[0] = nu_data[0] * correction_vec


    print ('z vectors: ', bincentermed_pops)
    print ('nu_data:',nu_data)
    print ('nu_err_pois:',nu_err_pois)
    print ('sigz2_data:',sigz2_data)
    print ('sigz2_err_pois:',sigz2_err_pois)
    if gp.tilt:
        print('sigRz2_data = ', sigRz2_data)
        print('sigRz2_err_pois = ', sigRz2_err_pois)


    Ntr_used = np.array([len(z_data_zcut[ii]) for ii in range(0, len(z_data_zcut))]) # Total number of tracers used (& binned)

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




    prior_sigma = 2.0
    gp.nu_C_max = 2*max(nu_0_values) + prior_sigma*max(nu_0_errs) # 5 sigma either way #10/3 bodge 2*
    gp.nu_C_min = min(nu_0_values) - prior_sigma*max(nu_0_errs)

    gp.sigz_C_max = max(sigz_0_values) + prior_sigma*max(sigz_0_errs)
    gp.sigz_C_min = max(min(sigz_0_values) - prior_sigma*max(sigz_0_errs), 0.)

    #Set priors for the exponential disks #REDO THESE PRIORS
    for jter in range(0, gp.N_nu_model_exps):
        gp.nu_exp_sum_priors[jter][0] = max(nu_0_values)
        gp.nu_exp_sum_priors[jter][1] = 0.1 * max(nu_0_values)


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


















#        #TEMPORARY TEST 25/01 increasing errors on sigRz2
#        sigz2_err_pois_tmp = sigz2_err_pois_tmp * 1.1
#
#        #Temportary test 15/02 increasing errors on nu
#        nu_err_pois_tmp = nu_err_pois_tmp * 1.2

        #TEMPORARY TEST 25/01 smoothing sigz2 data
        #sigz2_data_tmp = savgol_filter(sigz2_data_tmp, 5, 1)
        #nu_data_tmp = savgol_filter(nu_data_tmp, 5, 1)

        ##Set 02 Frank
        #sigz2_data_tmp = np.array([ 1610.61804191,  1665.59877152,  1722.80205969,  1782.57043621,
        #                            1817.24953954,  1859.59282268,  1912.78591356,  1941.68893828,
        #                            2002.24509393,  2018.79823251, 2063.27821223,  2126.09805272,
        #                            2161.57323368,  2220.85999367,  2281.50627278,  2326.66874511,
        #                            2442.85040529,  2525.20151286,  2643.39541919,  2765.76378331])


        ##Set 20 Frank
        #sigz2_data_tmp = np.array([ 1602.2919723 ,  1689.54799504,  1738.46953328,  1778.40284861,
        #                            1830.55026302,  1856.55686437,  1900.62203068,  1936.93570973,
        #                            2015.7846952 ,  1995.81055793, 2064.55710872,  2093.60970765,
        #                            2148.40595152,  2205.95909336,  2252.28646799,  2339.97492179,
        #                            2417.97394751,  2512.25825913,  2646.4444893 ,  2835.21902204])


        ##bins 11-14, 16, 19 transplant
        #sigz2_data_tmp = np.array([ 1602.2919723 ,  1689.54799504,  1738.46953328,  1778.40284861,
        #                            1830.55026302,  1856.55686437,  1900.62203068,  1936.93570973,
        #                            2015.7846952 ,  1995.81055793,  2063.27821223,
        #                            2100.0,  #11
        #                            2123.0,  #12
        #                            2216.0,  #13
        #                            2245.0,  #14
        #                            2326.66874511,
        #                            2411.0,  #16,
        #                            2525.20151286,  2643.39541919,
        #                            2806.0]) #19


        ##bins 11-14, 19 transplant
        #sigz2_data_tmp = np.array([ 1602.2919723 ,  1689.54799504,  1738.46953328,  1778.40284861,
        #                            1830.55026302,  1856.55686437,  1900.62203068,  1936.93570973,
        #                            2015.7846952 ,  1995.81055793,  2063.27821223,
        #                            2100.0,  #11
        #                            2123.0,  #12
        #                            2216.0,  #13
        #                            2245.0,  #14
        #                            2326.66874511, 2442.85040529, 2525.20151286,  2643.39541919,
        #                            2806.0]) #19

        #Bin 16 transplant
        #sigz2_data_tmp = np.array([ 1602.2919723 ,  1689.54799504,  1738.46953328,  1778.40284861,
        #                            1830.55026302,  1856.55686437,  1900.62203068,  1936.93570973,
        #                            2015.7846952 ,  1995.81055793,  2063.27821223,  2126.09805272,
        #                            2161.57323368,  2220.85999367,  2281.50627278,  2326.66874511,
        #                            2411.0, #16
        #                            2525.20151286,  2643.39541919,  2765.76378331])


        ##Bin 19 transplant
        #sigz2_data_tmp = np.array([ 1602.2919723 ,  1689.54799504,  1738.46953328,  1778.40284861,
        #                            1830.55026302,  1856.55686437,  1900.62203068,  1936.93570973,
        #                            2015.7846952 ,  1995.81055793,  2063.27821223,  2126.09805272,
        #                            2161.57323368,  2220.85999367,  2281.50627278,  2326.66874511,
        #                            2442.85040529,  2525.20151286,  2643.39541919,
        #                            2806.0])

        ##Bin 16 & 19 transplant transplant
        #sigz2_data_tmp = np.array([ 1602.2919723 ,  1689.54799504,  1738.46953328,  1778.40284861,
        #                            1830.55026302,  1856.55686437,  1900.62203068,  1936.93570973,
        #                            2015.7846952 ,  1995.81055793,  2063.27821223,  2126.09805272,
        #                            2161.57323368,  2220.85999367,  2281.50627278,  2326.66874511,
        #                            2411.0, #16
        #                            2525.20151286,  2643.39541919,
        #                            2807.0]) #19



        #Test 17/02
#        pdb.set_trace()
#        sigz2_err_pois_tmp = sigz2_err_pois_tmp*0.1
#        sigz2_err_pois_cumu = [sigz2_err_pois_tmp[0]]
#        for jter in range(1, len(sigz2_err_pois_tmp)):
#            nu_factor_i = (nu_data_tmp[jter-1]/nu_data_tmp[jter])
#            #nu_factor_i = 1.0
#            print('jter = ', jter, ', nu_factor_i = ', nu_factor_i)
#            sigz2_err_pois_cumu.append(sigz2_err_pois_tmp[jter] + nu_factor_i * sigz2_err_pois_cumu[jter-1])
#            print('jter = ', jter, ', sigz2_err_pois_tmp = ', sigz2_err_pois_tmp[jter], sigz2_err_pois_cumu[jter])
#
#        print('sigz2_err_pois_tmp = ', sigz2_err_pois_tmp)
#        print('sigz2_err_pois_cumu = ', sigz2_err_pois_cumu)
#        pdb.set_trace()
#        #sigz2_err_pois_cumu = np.array([(sigz2_err_pois_tmp[ii] + sigz2_err_pois_cumu[ii-1]*nu_data_tmp[ii-1]/nu_data_tmp[ii]) for ii in range(1, len(sigz2_err_pois_tmp))])
#        #sigz2_err_pois_cumu = np.array([(nu_data_tmp[ii-1]/nu_data_tmp[ii]) for ii in range(1, len(sigz2_err_pois_tmp))])
#
#
#
#        sigz2_err_pois_tmp = np.append(sigz2_err_pois_0, sigz2_err_pois_cumu)
#        pdb.set_trace()


        #1 SD up and down test

        #sigz2_data_tmp = np.array([1623.25139087,  1688.37444732,  1743.97698327,  1791.7950385,   1834.83250905,
        #                        1874.45496845,  1912.79200978,  1950.16273803,  1988.93871516,  2028.12197511,
        #                        2070.10042631,  2114.74191993,  2163.33053234,  2217.17667516,  2277.46685436,
        #                        2346.35939292,  2427.9963064,  2526.79694338,  2652.82556956,  2829.90130777])


        #sigz2_data_tmp = np.array([ 1602.23761271,  1666.21632173,  1721.17726028,  1768.4715971,   1810.82516345,
        #                            1850.10655084,  1887.86568523,  1924.76016596,  1962.50205248,  2001.94726972,
        #                            2043.04089448,  2086.85851691,  2134.98187854,  2188.05048591,  2247.5452998,
        #                            2315.84554045,  2395.95875184,  2493.67936848,  2618.15790231,  2793.6288009])


        #sigz2_data_tmp = np.array([1633.7582799493587, 1699.4535101136034, 1755.3768447575219, 1803.4567592091764,
        #                            1846.8361818550111, 1886.6291772629531, 1925.2551720466179, 1962.8640240668881,
        #                            2002.157046506449, 2041.2093277957306, 2083.6301922310295, 2128.6836214397135,
        #                            2177.5048592381504, 2231.7397697874244, 2292.4276316462337, 2361.6163191561204,
        #                            2444.0150836842595, 2543.3557308243558, 2670.1594031782442, 2848.037561206846])

        #sigz2_data_tmp = np.array([1591.7307236339855, 1655.1372589329853, 1709.777398791415, 1756.8098763914988,
        #                            1798.8214906517383, 1837.9323420300918, 1875.4025229602846, 1912.058879918701,
        #                            1949.2837211313911, 1988.8599170338089, 2029.5111285660084, 2072.9168153971295,
        #                            2120.8075516442682, 2173.48739127811, 2232.584522511971, 2300.5886142165477,
        #                            2379.9399745572236, 2477.120581036183, 2600.8240686920881, 2775.4925474592505])
