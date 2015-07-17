#!/usr/bin/env ipython3

##
# @file
# all functions called directly from gravimage
# disc version, done

# (c) ETHZ 2013 Pascal Steger, psteger@phys.ethz.ch

#from types import *
import pdb
import numpy as np
import gl_helper as gh
import sys
from mpi4py import MPI

def chi2red(model, data, sig, dof):
    # if Degrees Of Freedom = 1, return non-reduced chi2
    model = np.array(model)
    data  = np.array(data)
    sig   = np.array(sig)
    chired = np.sum(((model-data)**2./sig**2.)/dof)
    return chired
## \fn chi2red(model, data, sig, dof)
# determine 'reduced chi2'
# @param model
# @param data
# @param sig spread
# @param dof Degrees Of Freedom


def calc_chi2(profs, gp):

    #hwmess = "calc_chi2 running on process %d of %d on %s.\n"
    #myrank = MPI.COMM_WORLD.Get_rank()
    #nprocs = MPI.COMM_WORLD.Get_size()
    #procnm = MPI.Get_processor_name()
    #sys.stdout.write(hwmess % (myrank, nprocs, procnm))

    if gp.map_priors:
        return 1.0

    #chi2_nu = [[].append(None) for ii in range(0,gp.ntracer_pops)]
    #chi2_sigz2 = [[].append(None) for ii in range(0,gp.ntracer_pops)]
    #chi2_tilt = [[].append(None) for ii in range(0,gp.ntracer_pops)]
    chi2_nu = 0
    chi2_sigz2 = 0
    chi2_tilt = 0


    #Tracer population comparison
    for pop in range(0, gp.ntracer_pops):
        gh.LOG(1, ' pop = ', pop)
        #Check Monotonicity
        if gp.monotonic_rho:
            kz_rho_DM_vec_temp = profs.get_prof('kz_rho_DM_vec', pop)
            if min(kz_rho_DM_vec_temp) < 0:
                gh.LOG(1,'Negative kz_rho_DM in gl_chi')
                raise ValueError('negative value in kz_rho_DM array')

        if gp.monotonic_nu:
            kz_nu_vec_temp = profs.get_prof('kz_nu_vec', pop)
            if min(kz_nu_vec_temp) < 0:
                gh.LOG(1,'Negative kz_nu_DM in gl_chi')
                raise ValueError('negative value in kz_nu array')

        #If monotonicity passed calculate chi2
        nudat    = gp.dat.nu[pop]
        nuerr    = gp.dat.nuerr[pop]+profs.hyper_nu  # adding hyperparam to error
        numodel  = profs.get_prof('nu_vec', pop)
        chi2_nu_tmp = chi2red(numodel, nudat, nuerr, 1.) #reduced dof = gp.nbins
        chi2_nu += chi2_nu_tmp
        #print('pop = ', pop)
        #print('nudat = ', nudat)
        #print('numodel = ', numodel)
        #print('nuerr = ', nuerr)
        gh.LOG(1, ' chi2_nu = ', chi2_nu_tmp)

        if not gp.chi2_nu_converged and not gp.plotting_flag:
            continue # with pop loop

        sigz2dat    = gp.dat.sigz2[pop]    # [km/s]
        sigz2err    = gp.dat.sigz2err[pop]+profs.hyper_sigz2  # [km/s]
        sigz2_model = profs.get_prof('sigz2_vec', pop)
        chi2_sigz2_tmp  = chi2red(sigz2_model, sigz2dat, sigz2err, 1.) #reduced dof = gp.nbins
        if chi2_sigz2_tmp == np.inf:
            print('chi2_sig has become infinite')
        chi2_sigz2 += chi2_sigz2_tmp
        gh.LOG(1, '  chi2_sigz2  = ', chi2_sigz2_tmp)

        if gp.tilt:
            sigmaRzdat  = gp.dat.tilt
            sigmaRzerr = gp.dat.tilterr
            sigmaRz_model = profs.get_prof('sigmaRz_vec', pop)
            chi2_tilt_tmp = chi2red(sigmaRz_model, sigmaRzdat, sigmaRzerr, 1.) #reduced dof = gp.nbins
            chi2_tilt += chi2_tilt_tmp
            #print ('sigmaRz2dat:',sigmaRz2dat)
            #print ('sigmaRz2_model:',sigmaRz2_model)
            #print ('sigmaRz2err:',sigmaRz2err)
            chi2 = (chi2_nu+chi2_sigz2+chi2_tilt)#/3.
            #print ('chi2:',chi2,'chi2_tilt:',chi2_tilt)
            gh.LOG(1, '  chi2_tilt  = ', chi2_tilt_tmp)
        else:
            chi2 = (chi2_nu+chi2_sigz2)#/2.

    # switch to chi2_sig calculation too, if converged on Sig
    if not gp.chi2_nu_converged:
        chi2 = chi2_nu*1.
        chi2 *= 10
        if chi2 < gp.chi2_switch:
            gp.chi2_switch_counter +=1
            print('chi2 less than switch found, ', gp.chi2_switch_counter, ' out of ', gp.chi2_switch_mincount, ' needed.')
        if gp.chi2_switch_counter>= gp.chi2_switch_mincount:
            gh.LOG(1, 'nu burn-in finished, switching on sigma')
            gp.chi2_nu_converged = True

        #if chi2 < gp.chi2_switch:
        #    pdb.set_trace()
        #    gh.LOG(1, 'nu burn-in finished, switching on sigma')
        #    gp.chi2_nu_converged = True

    return chi2
## \fn calc_chi2(profs)
# Calculate chi^2
# @param profs profiles for rho, M, nu_i, beta_i, sig_i, kap_i
