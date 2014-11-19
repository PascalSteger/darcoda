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

def write_disc_output_files(Rbin, Binmin, Binmax, nudat, nuerr, Sigdat, Sigerr, Mdat, Merr, sigdat, sigerr, scales, gp):
#    pdb.set_trace()
    for pop in range(gp.pops+1):
        # write scales
        crscale = open(gp.files.get_scale_file(pop), 'w')
        print('# Zscale in [pc], surfdens_central (=dens0) in [Munit/Zscale^2], and in [Munit/pc^2], and totmass_tracers [Munit], and max(sigma_LOS) in [km/s]', file=crscale)
        for b in range(5):
            # Zscale, Dens0pc, totmass_tracers, nu0pc, maxsiglos
            # Dens0pc is *not* density at halflight radius as in spherical
            # code, but rather 3D tracer density at center
            print(scales[pop][b], file=crscale)
        crscale.close()

        # write tracer densities
        f_Sig = open(gp.files.Sigfiles[pop], 'w')
        print('Rbin [Zscale];','Binmin [Zscale];','Binmax [Zscale];',\
              'Sig(R) [Munit/pc^2];','error [1]', file=f_Sig)
        for b in range(gp.nipol):
            print(Rbin[b], Binmin[b], Binmax[b], \
                  Sigdat[pop][b], Sigerr[pop][b], file=f_Sig)
        f_Sig.close()

        # write tracer densities 3D
        f_nu = open(gp.files.nufiles[pop], 'w')
        print('Rbin [Zscale];','Binmin [Zscale];','Binmax [Zscale];',\
              'nu(R) [Munit/pc^3];','error [1]', file=f_nu)
        for b in range(gp.nipol):
            print(Rbin[b], Binmin[b], Binmax[b], \
                  nudat[pop][b], nuerr[pop][b], file=f_nu)
        f_nu.close()

        # write velocity dispersion
        f_sig = open(gp.files.sigfiles[pop],'w')
        print('R [Zscale];','Binmin [Zscale];','Binmax [Zscale];',\
              'sigma_p(R) [maxsiglos];','error [km/s]', file=f_sig)
        for b in range(gp.nipol):
            print(Rbin[b], Binmin[b], Binmax[b], \
                  sigdat[pop][b], sigerr[pop][b], file=f_sig)
        f_sig.close()

        # write enclosed mass
        f_mass = open(gp.files.massfiles[pop],'w')
        print('R [Zscale];','Binmin [Zscale];','Binmax [Zscale];',\
              'M(<Binmax) [Munit];','error [Munit]', file=f_mass)
        for b in range(gp.nipol):
            print(Rbin[b], Binmin[b], Binmax[b], \
                  Mdat[pop][b], Merr[pop][b], file=f_mass)
        f_mass.close()

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
    pdb.set_trace()

    external_file='/home/hsilverw/LoDaM/darcoda/gravlite/simplenu/simplenu_sigz_raw.dat'
    external_data = np.loadtxt(external_file)
    z_data = external_data[:, 0]
    v_data = abs(external_data[:, 1])

    # Bin data, constant number of tracers per bin.
    # First, calculate binmins, binmaxes, and bin centres
    binmin, binmax, bincentermed = gh.bin_r_const_tracers(z_data, gp.nbins)
    # Then calculate tracer number density
    nu_data, nu_err_data, sig_data, sig_err_data = gh.nu_sig_from_bins(binmin, binmax, z_data, v_data)








    #HS Working Line
    #Everthing below this line is old, and hasn't been considered for keeping
    #deletion, or modification
    # ----------------------------------------------------------------------

    import gr_params
    gpr = gr_params.Params(gp)
    if gpr.showplots:
        nuscaleb = nu_zth[np.argmin(np.abs(zth-z0))]
        plt.loglog(zth, nu_zth/nuscaleb, 'b.-')
        nuscaler = nu_dat_bin1[np.argmin(np.abs(zth-z0))]
        plt.loglog(zth, nu_dat_bin1/nuscaler, 'r.-')
        # pdb.set_trace()

    Sig_dat_bin1 = np.cumsum(nu_dat_bin1)
    Sig_dat_err_bin1 = np.sqrt(Sig_dat_bin1)
    Mrdat1 = np.cumsum(Sig_dat_bin1)
    Mrerr1 = Mrdat1*Sig_dat_err_bin1/Sig_dat_bin1

    scales=[[],[],[]]
    scales[1].append(z0) # [pc]
    scales[1].append(Sig_dat_bin1[0])
    scales[1].append(Mrdat1[-1])
    scales[1].append(nu_dat_bin1[0])
    scales[1].append(max(sig_dat_bin1))

    # start analysis of "all stars" with only component 1,
    # append to it later if more populations required
    z_dat0 = z_dat1   # [pc]
    vz_dat0 = vz_dat1 # [km/s]

    if gp.pops == 2:
        # enforce observational constraints on z<z_max
        sel = (zstar2 < zmax)
        z_dat2  = zstar2[sel]
        vz_dat2 = vzstar2[sel]

        # cut zero velocities:
        sel = (abs(vz_dat2) > 0)
        z_dat2  = z_dat2[sel]
        vz_dat2 = vz_dat2[sel]

        # Calulate binned data (for plots/binned analysis):
        binmin2, binmax2, z_dat_bin2, sig_dat_bin2, count_bin2 = gh.binsmooth(z_dat2, vz_dat2, \
                                                                              zmin, zmax, gp.nipol, 0.)
        sig_dat_err_bin2 = np.sqrt(sig_dat_bin2) # Poissonian errors

        nu_dat_bin2, nu_dat_err_bin2 = gh.bincount(z_dat2, binmax2)
        nu_dat_bin2 /= (binmax2-binmin2)
        nu_dat_err_bin2 /= (binmax2-binmin2)

        Sig_dat_bin2 = np.cumsum(nu_dat_bin2)
        Sig_dat_err_bin2 = np.sqrt(Sig_dat_bin2)
        Mrdat2 = np.cumsum(nu_dat_bin2)
        Mrerr2 = np.sqrt(Mrdat2)

        scales[2].append(z02) # [pc]
        scales[2].append(Sig_dat_bin2[0])
        scales[2].append(Mrdat2[-1])
        scales[2].append(nu_dat_bin2[0]) # normalize by max density of first bin, rather
        scales[2].append(max(sig_dat_bin2))

        # calculate properties for all pop together with stacked values
        z_dat0 = np.hstack([z_dat1, z_dat2])
        vz_dat0 = np.hstack([vz_dat1, vz_dat2])

    # Calulate binned data (for plots/binned anal.). old way, linear spacings, no const #particles/bin
    binmin0, binmax0, z_dat_bin0, sig_dat_bin0, count_bin0 = gh.binsmooth(z_dat0, vz_dat0, \
                                                                          zmin, zmax, gp.nipol, 0.)
    sig_dat_err_bin0 = np.sqrt(sig_dat_bin0)
    # binmin, binmax, z_dat_bin = gh.bin_r_const_tracers(z_dat, gp.nipol) # TODO: enable, get sig2

    nu_dat_bin0, nu_dat_err_bin0 = gh.bincount(z_dat0, binmax0)
    nu_dat_bin0 /= (binmax0-binmin0)
    nu_dat_err_bin0 /= (binmax0-binmin0)

    Sig_dat_bin0 = np.cumsum(nu_dat_bin0)
    Sig_dat_err_bin0 = np.sqrt(Sig_dat_bin0)
    #renorm0 = max(nu_dat_bin0)

    xip = np.copy(z_dat_bin0)                        # [pc]
    Mrdat0   = K*xip/np.sqrt(xip**2.+D**2.) / (2.0*np.pi*gu.G1__pcMsun_1km2s_2)
    Mrerr0   = Mrdat0*nu_dat_err_bin0/nu_dat_bin0

    scales[0].append(D)  # [pc]
    scales[0].append(Sig_dat_bin0[0])
    scales[0].append(Mrdat0[-1])
    scales[0].append(nu_dat_bin0[0])
    scales[0].append(max(sig_dat_bin0))

    rmin = binmin0/scales[0][0] # [pc]
    rbin = xip/scales[0][0]     # [pc]
    rmax = binmax0/scales[0][0] # [pc]

    # store parameters for output
    # normalized by scale values
    nudat = []
    nudat.append(nu_dat_bin0/scales[0][3])       # [Msun/pc^3]
    nudat.append(nu_dat_bin1/scales[1][3])
    if gp.pops == 2:
        nudat.append(nu_dat_bin2/scales[2][3])

    nuerr = []
    nuerr.append(nu_dat_err_bin0/scales[0][3])   # [Msun/pc^3]
    nuerr.append(nu_dat_err_bin1/scales[1][3])
    if gp.pops == 2:
        nuerr.append(nu_dat_err_bin2/scales[2][3])

    Mrdat = []
    Mrdat.append(Mrdat0/scales[0][2])            # [Msun]
    Mrdat.append(Mrdat1/scales[1][2])
    if gp.pops == 2:
        Mrdat.append(Mrdat2/scales[2][2])

    Mrerr = []
    Mrerr.append(Mrerr0/scales[0][2])            # [Msun]
    Mrerr.append(Mrerr1/scales[1][2])
    if gp.pops == 2:
        Mrerr.append(Mrerr2/scales[2][2])

    Sigdat = []
    Sigdat.append(Sig_dat_bin0/scales[0][1])
    Sigdat.append(Sig_dat_bin1/scales[1][1])
    if gp.pops == 2:
        Sigdat.append(Sig_dat_bin2/scales[2][1])

    Sigerr = []
    Sigerr.append(Sig_dat_err_bin0/scales[0][1])
    Sigerr.append(Sig_dat_err_bin1/scales[1][1])
    if gp.pops == 2:
        Sigerr.append(Sig_dat_err_bin2/scales[2][1])

    sigdat = []
    sigdat.append(sig_dat_bin0/scales[0][4])     # [km/s]
    sigdat.append(sig_dat_bin1/scales[1][4])
    if gp.pops == 2:
        sigdat.append(sig_dat_bin2/scales[2][4])

    sigerr = []
    sigerr.append(sig_dat_err_bin0/scales[0][4]) # [km/s]
    sigerr.append(sig_dat_err_bin1/scales[1][4])
    if gp.pops == 2:
        sigerr.append(sig_dat_err_bin2/scales[2][4])

#    pdb.set_trace()

    write_disc_output_files(rbin, rmin, rmax, nudat, nuerr, \
                            Sigdat, Sigerr, Mrdat, Mrerr,\
                            sigdat, sigerr, scales, gp)

    return gp.dat

## \fn run(gp)
# generate disc data from analytic form, return 3D densities, delta (=tilt)
# write it to files
# @param gp global parameters

if __name__=="__main__":
    import gl_params
    gp = gl_params.Params()
    run(gp)
