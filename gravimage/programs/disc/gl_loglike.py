#!/usr/bin/env ipython3

## @file
# define log likelihood function to be called by pyMultinest
# disc version

import numpy as np
import pdb
import gl_helper as gh
from gl_class_profiles import Profiles
from gl_priors import check_bprior, check_tilt
from gl_chi import calc_chi2
import gl_physics as phys
from pylab import *
#ion()
import time

def geom_loglike(cube, ndim, nparams, gp):
    tmp_profs = Profiles(gp.ntracer_pops, gp.nbins)#, gp.nbaryon_pops, gp.nbaryon_params)#, gp.nrhonu, )

    # Export z values to the profile
    for pop in range(0, gp.ntracer_pops):
        tmp_profs.set_prof('z_vecs', gp.z_bincenter_vecs[pop], pop, gp)
    tmp_profs.set_prof('z_vecs_comb_w0', gp.z_all_pts_sorted, pop, gp)



    #[tmp_profs.z_vecs[pop] = gp.z_bincenter_vecs[pop] for pop in range(0, gp.ntracer_pops)]
    #Normalisation constant C for sigz calculation
    off = 0
    offstep = gp.ntracer_pops
    norm_C = cube[off:off+offstep]
    for t_pop in range(0, gp.ntracer_pops):
        tmp_profs.norm_C[t_pop] = norm_C[t_pop]
    off += offstep
    #Dark Matter rho parameters (rho_C, kz_C, kz_vector)


    if gp.scan_rhonu_space:
        tmp_rho_DM_allz = np.array(cube[off:off+offstep])
    elif gp.darkmattermodel == 'const_dm':
        offstep = 1
        rho_DM_params = np.array(cube[off:off+offstep])
        rho_DM_C = rho_DM_params[0]
        tmp_rho_DM_allz = rho_DM_C * np.ones(gp.nrhonu)
    elif gp.darkmattermodel == 'ConstPlusDD':
        offstep = 3
        rho_DM_params = np.array(cube[off:off+offstep])
        rho_DM_C = rho_DM_params[0]
        tmp_rho_DM_allz = phys.rho_dm_simplenu(gp.z_all_pts_sorted, rho_DM_params)
    else:
        offstep = gp.nrhonu + 1
        rho_DM_params = np.array(cube[off:off+offstep])
        rho_DM_C = rho_DM_params[0] #rho_C
        kz_rho_DM_allz = rho_DM_params[1:] #kz for rho across all z points [0, bin_centres]
        if gp.monotonic_rho:         #outputs rho across all points:
            tmp_rho_DM_allz = phys.rho(gp.z_all_pts_sorted, abs(kz_rho_DM_allz), rho_DM_C)
        else:
            tmp_rho_DM_allz = phys.rho(gp.z_all_pts_sorted, kz_rho_DM_allz, rho_DM_C)

        tmp_profs.kz_rho_DM_C = kz_rho_DM_allz[0]
        tmp_profs.set_prof('kz_rho_DM_vec', kz_rho_DM_allz[1:], 0, gp)

    tmp_rho_total_allz = tmp_rho_DM_allz*1.0 # Add DM to total mass density

    tmp_profs.rho_DM_C = tmp_rho_DM_allz[0]
    tmp_profs.set_prof('rho_DM_vec', tmp_rho_DM_allz[1:], 0, gp)
    off += offstep

    #Baryons
    tmp_rho_totbaryon_allz = np.zeros(len(tmp_rho_total_allz))

    for baryon_pop in range(0, gp.nbaryon_pops):
        offstep = gp.nbaryon_params
        baryon_params = np.array(cube[off:off+offstep])
        if gp.baryonmodel == 'simplenu_baryon':
            tmp_rho_baryon_allz = phys.rho_baryon_simplenu(gp.z_all_pts_sorted, baryon_params)
        elif gp.baryonmodel == 'kz_baryon':
            print('Not implemented yet')
        tmp_profs.rho_baryon_C = tmp_rho_baryon_allz[0]
        tmp_profs.set_prof('rho_baryon_vec', tmp_rho_baryon_allz[1:], baryon_pop, gp)

        tmp_rho_totbaryon_allz += tmp_rho_baryon_allz # add this pop's density to baryon total
        tmp_rho_total_allz += tmp_rho_baryon_allz # add the baryon density to the total density

        off += offstep


    tmp_profs.rho_total_C = tmp_rho_total_allz[0]
    tmp_profs.set_prof('rho_total_vec', tmp_rho_total_allz[1:], 0, gp)

    #Tracer params, nu_C, kz_nu_C, kz_nu_vector
    tmp_nu_allz=[[].append(None) for ii in range(0,gp.ntracer_pops)]

    for t_pop in range(0, gp.ntracer_pops):
        offstep = gp.nbins[t_pop] + 1 + 1 #kz on bincenters, and zC=0, and nu_C
        if gp.scan_rhonu_space:
            tmp_nu_allz = np.array(cube[off:off+offstep])
        else:
            tracer_params = np.array(cube[off:off+offstep])
            nu_C = tracer_params[0]
            kz_nu_allz = tracer_params[1:] #kz for rho across all z points [0, bin_centres]
            z_points_tmp = np.append(0., gp.z_bincenter_vecs[t_pop])
            if gp.monotonic_nu:  #outputs nu across z points for that population:
                tmp_nu_allz[t_pop] = phys.rho(z_points_tmp, abs(kz_nu_allz), nu_C)
            else:
                tmp_nu_allz[t_pop] = phys.rho(z_points_tmp, kz_nu_allz, nu_C)

            tmp_profs.kz_nu_C[t_pop] = kz_nu_allz[0]
            tmp_profs.set_prof('kz_nu_vecs', kz_nu_allz[1:], t_pop, gp)

        tmp_profs.nu_C[t_pop] = tmp_nu_allz[t_pop][0]
        tmp_profs.set_prof('nu_vecs', tmp_nu_allz[t_pop][1:], t_pop, gp)
        off += offstep


    # Tilt term
    tilt_params=[[].append(None) for ii in range(0,gp.ntracer_pops)]
    if gp.tilt:
        for t_pop in range(0, gp.ntracer_pops):
            offstep = gp.ntilt_params
            tilt_params[t_pop] = np.array(cube[off:off+offstep])
            off += offstep

    #Hyperparameters
    if gp.hyperparams == True:
        for t_pop in range(0, gp.ntracer_pops):
            offstep = 1
            #hyper_nu = cube[off:off+offstep][0]/gp.dat.meannuerr # 0.1->10
            tmp_profs.hyper_nu = cube[off:off+offstep]
            off += offstep
            offstep = 1
            #hyper_sigz2 = cube[off:off+offstep][0]/gp.dat.meansigz2err
            tmp_profs.hyper_sigz2 = cube[off:off+offstep]
            off += offstep
            ## To set the hyperparameters to effectively the same value,
            ##  uncomment and comment out tmp.profs... = ...  above
            #hyper_average = (hyper_nu+hyper_sigz2)/2. # normalized number: 0.1->10
            #tmp_profs.hyper_nu = [gp.dat.meannuerr*hyper_average]
            #tmp_profs.hyper_sigz2 = [gp.dat.meansigz2err*hyper_average]
            #print ('hyper_average:',hyper_average)

    if off != gp.ndim:
        gh.LOG(1,'wrong subscripts in gl_class_cube')
        print ('in loglike',off,gp.ndim)
        raise Exception('wrong subscripts in gl_class_cube')



    #Calculate Sigma (surface density)
    Sig_DM_allz = phys.Sig(gp.z_all_pts_sorted, tmp_rho_DM_allz) # Dark matter Sigma
    Sig_baryon_allz = phys.Sig(gp.z_all_pts_sorted, tmp_rho_totbaryon_allz) # total baryon Sigma
    Sig_total_allz = phys.Sig(gp.z_all_pts_sorted, tmp_rho_total_allz) # total Sigma

    tmp_profs.Sig_DM_C = Sig_DM_allz[0]
    tmp_profs.set_prof('Sig_DM_vec', Sig_DM_allz[1:], 0, gp)

    tmp_profs.Sig_baryon_C = Sig_baryon_allz[0]
    tmp_profs.set_prof('Sig_baryon_vec', Sig_baryon_allz[1:], 0, gp)

    tmp_profs.Sig_total_C = Sig_total_allz[0]
    tmp_profs.set_prof('Sig_total_vec', Sig_total_allz[1:], 0, gp)


    #Calculate tilt for each population
    sigmaRz_allz = [[].append(None) for ii in range(0,gp.ntracer_pops)]
    tilt_allz = [[].append(None) for ii in range(0,gp.ntracer_pops)]
    for t_pop in range(0, gp.ntracer_pops):
        z_points_tmp = np.append(0., gp.z_bincenter_vecs[t_pop])
        if gp.tilt:
            A = tilt_params[t_pop][0]
            n = tilt_params[t_pop][1]
            R = tilt_params[t_pop][2]
            if n<0:
                print('z_points_tmp = ', z_points_tmp, ', A = ', A, ', n = ', n, ', R = ', R)

            sigmaRz_allz[t_pop] = A*(z_points_tmp*1000.)**n #all z for this population
            tmp_profs.set_prof('sigmaRz_vecs', sigmaRz_allz[t_pop][1:], t_pop, gp)
            tilt_allz[t_pop] = sigmaRz_allz[t_pop]*(1./gp.Rsun - 2./R)
        else:
            tilt_allz[t_pop] = np.zeros(gp.nbins[t_pop]+1)
        tmp_profs.set_prof('tilt_vecs', tilt_allz[t_pop][1:], t_pop, gp)
    #print ('gp.tilt:',gp.tilt,' tilt_allz:',tilt_allz)

    #Calculate sigma (velocity dispersion)
    for t_pop in range(0, gp.ntracer_pops):
        Sig_total_popz_tmp = Sig_total_allz[gp.z_vec_masks[t_pop]] #Sig_total at bin centres for this pop
        z_points_tmp = np.append(0., gp.z_bincenter_vecs[t_pop])
        try:
            sigz2_vec = phys.sigz2(z_points_tmp, Sig_total_popz_tmp, tilt_allz[t_pop], tmp_nu_allz[t_pop], norm_C[t_pop])
        except ValueError:
            raise ValueError('negative value in sig2 array')
            print('tracer pop = ', t_pop)
            print('tilt params = ', tilt_params[t_pop])
            print('Sig_total_allz = ', Sig_total_allz)
            print('tilt_all_z[t_pop] = ', tilt_allz[t_pop])
            print('tmp_nu_allz[t_pop] = ', tmp_nu_allz[t_pop])
            print('norm = ', norm)
            return
        tmp_profs.set_prof('sigz2_vecs', sigz2_vec[1:], t_pop, gp)

    if (np.isinf(sigz2_vec)).any():
        print('Negative sigz2_vec')
        print('tracer pop = ', t_pop)
        print('tilt params = ', tilt_params[t_pop])
        print('Sig_total_allz = ', Sig_total_allz)
        print('tilt_all_z[t_pop] = ', tilt_allz[t_pop])
        print('tmp_nu_allz[t_pop] = ', tmp_nu_allz[t_pop])
        print('norm = ', norm)
        pdb.set_trace()

    # determine log likelihood
    chi2 = calc_chi2(tmp_profs, gp)

    gh.LOG(1, '   log L = ', -chi2/2.)
    tmp_profs.chi2 = chi2

    return tmp_profs   # from   likelihood L = exp(-\chi^2/2), want log of that
## \fn geom_loglike(cube, ndim, nparams, gp)
# define log likelihood function to be called by pyMultinest and plot_profiles
# disc version
# @param cube parameter cube as defined by gl_class_cube, in physical space
# @param ndim number of dimensions, needed as argument by multinest
# @param nparams number of parameters, needed as argument by multinest
# @param gp global parameters
