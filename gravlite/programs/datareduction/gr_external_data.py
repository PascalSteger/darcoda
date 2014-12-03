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

def write_disc_output_files(Bincenter, Binmin, Binmax, nudat, nuerr, sigdat, sigerr, gp):

    # write tracer densities 3D
    file_nu = open(gp.files.nufiles[0], 'w')
    print('Bin centres z [kpc];','Binmin z [kpc];','Binmax z [kpc];','nu(z) [#/pc^3];','error [#/pc^3]', file=file_nu)
    for b in range(gp.nbins):
        print(Bincenter[b], Binmin[b], Binmax[b], nudat[b], nuerr[b], file=file_nu)
    file_nu.close()

    file_sig = open(gp.files.sigfiles[0],'w')
    print('Bin centres z [kpc];','Binmin z [kpc];','Binmax [kpc];', 'sigma_z(z) [km/s];', 'error [km/s]', file=file_sig)
    for b in range(gp.nbins):
        print(Bincenter[b], Binmin[b], Binmax[b], sigdat[b], sigerr[b], file=file_sig)
    file_sig.close()

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
    global K,C,D,F, zth, zp_kz, zmin, zmax, z0, z02


    external_file='/home/hsilverw/LoDaM/darcoda/gravlite/Data_Sets/simplenu/simplenu_sigz_raw.dat'
    external_data = np.loadtxt(external_file)
    z_data = external_data[:, 0] #[kpc]
    v_data = external_data[:, 1] #[km/s]

    #Population Splitting
    if gp.ntracer_pops > 1:
        print('More than one population')
        #DO POPULATION SPLITTING

    pdb.set_trace()
    # Bin data, constant number of tracers per bin.
    # First, calculate binmins, binmaxes, and bin centres
    binmin, binmax, bincentermed = gh.bin_r_const_tracers(z_data, gp.nbins)

    #JR binning
    bincentermed = np.array(range(1,24))*0.05
    binmin=np.array(range(1,24))*0.05 - 0.025
    binmin[0]=0.0
    binmax=np.array(range(1,24))*0.05 +0.025


    # Then calculate tracer number density [#stars/kpc^3], [#stars/kpc^3], [km/s], [km/s]
    nu_data, nu_err_data, sig_data, sig_err_data = gh.nu_sig_from_bins(binmin, binmax, z_data, v_data)

    write_disc_output_files(bincentermed, binmin, binmax, nu_data, nu_err_data, sig_data, sig_err_data, gp)

    import gr_params #WHAT DOES THIS DO.
    gpr = gr_params.Params(gp)
    if gpr.showplots:
        nuscaleb = nu_zth[np.argmin(np.abs(zth-z0))]
        plt.loglog(zth, nu_zth/nuscaleb, 'b.-')
        nuscaler = nu_dat_bin1[np.argmin(np.abs(zth-z0))]
        plt.loglog(zth, nu_dat_bin1/nuscaler, 'r.-')
        # pdb.set_trace()

    return gp.dat

## \fn run(gp)
# Read in data from JR's simple 1D mock
# Bin it, calculate nu and sigma, output it to file
# @param gp global parameters

if __name__=="__main__":
    import gl_params
    gp = gl_params.Params()
    run(gp)
