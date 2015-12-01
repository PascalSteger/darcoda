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
    chi2_vec = ((model-data)**2./sig**2.)/dof
    chired = np.sum(chi2_vec)
    return chired, chi2_vec
## \fn chi2red(model, data, sig, dof)
# determine 'reduced chi2'
# @param model
# @param data
# @param sig spread
# @param dof Degrees Of Freedom


def calc_chi2(profs, gp):

    #hwmess = "calc_chi2 running on process %d of %d on %s.\n"
    myrank = MPI.COMM_WORLD.Get_rank()
    nprocs = MPI.COMM_WORLD.Get_size()
    procnm = MPI.Get_processor_name()
    #sys.stdout.write(hwmess % (myrank, nprocs, procnm))



    #chi2_nu = [[].append(None) for ii in range(0,gp.ntracer_pops)]
    #chi2_sigz2 = [[].append(None) for ii in range(0,gp.ntracer_pops)]
    #chi2_tilt = [[].append(None) for ii in range(0,gp.ntracer_pops)]
    chi2_nu = 0
    chi2_sigz2 = 0
    chi2_sigRz2 = 0

    chi2_nu_vecs=[]
    chi2_sigz2_vecs=[]
    chi2_sigRz2_vecs=[]

    if gp.map_priors:
        chi2 = 1.0
        for pop in range(0, gp.ntracer_pops):
            chi2_nu_vecs.append( np.zeros(gp.nbins[pop]))
            chi2_sigz2_vecs.append( np.zeros(gp.nbins[pop]))
            chi2_sigRz2_vecs.append( np.zeros(gp.nbins[pop]))

        return chi2, chi2_nu_vecs, chi2_sigz2_vecs, chi2_sigRz2_vecs


    #Tracer population comparison
    for pop in range(0, gp.ntracer_pops):
        if gp.ntracer_pops > 1:
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


        #Calculate tracer density chi2 for population no. pop
        if gp.nu_model == 'kz_nu':
            nudat    = gp.dat.nu[pop]
            nuerr    = gp.dat.nuerr[pop]+profs.hyper_nu  # adding hyperparam to error
            numodel  = profs.get_prof('nu_vecs', pop)
            chi2_nu_tmp, chi2_nu_vec = chi2red(numodel, nudat, nuerr, 1.) #reduced dof = gp.nbins
            chi2_nu += chi2_nu_tmp
            chi2_nu_vecs.append(chi2_nu_vec)
            gh.LOG(1, ' chi2_nu = ', chi2_nu_tmp)
        else:
            chi2_nu_tmp = 0.
            chi2_nu = 0.
            chi2_nu_vecs.append(np.zeros(gp.nbins[pop]))
        #print('pop = ', pop)
        #print('nudat = ', nudat)
        #print('numodel = ', numodel)
        #print('nuerr = ', nuerr)


        #Calculate z-velocity dispersion chi2 for population no. pop
        sigz2dat    = gp.dat.sigz2[pop]    # [km/s]
        sigz2err    = gp.dat.sigz2err[pop]+profs.hyper_sigz2  # [km/s]
        sigz2_model = profs.get_prof('sigz2_vecs', pop)
        chi2_sigz2_tmp, chi2_sigz2_vec  = chi2red(sigz2_model, sigz2dat, sigz2err, 1.)
        if chi2_sigz2_tmp == np.inf:
            print('chi2_sig has become infinite')
        chi2_sigz2 += chi2_sigz2_tmp
        chi2_sigz2_vecs.append(chi2_sigz2_vec)
        gh.LOG(1, '  chi2_sigz2  = ', chi2_sigz2_tmp)

        #Calculate Rz-velocity dispersion chi2 for population no. pop
        if gp.tilt:
            sigmaRz2dat  = gp.dat.sigRz2[pop]
            sigmaRz2err = gp.dat.sigRz2err[pop]
            sigmaRz2_model = profs.get_prof('sigmaRz2_vecs', pop)
            chi2_sigRz2_tmp, chi2_sigRz2_vec = chi2red(sigmaRz2_model, sigmaRz2dat, sigmaRz2err, 1.)
            #print ('sigmaRz2dat:',sigmaRz2dat)
            #print ('sigmaRz2_model:',sigmaRz2_model)
            #print ('sigmaRz2err:',sigmaRz2err)
            #print ('chi2:',chi2,'chi2_tilt:',chi2_tilt)
            gh.LOG(1, '  chi2_sigRz2  = ', chi2_sigRz2_tmp)
        else:
            chi2_sigRz2_tmp =0.
            chi2_sigRz2_vec = np.zeros(gp.nbins[pop])
        chi2_sigRz2 += chi2_sigRz2_tmp
        chi2_sigRz2_vecs.append(chi2_sigRz2_vec)

    #Combine chi2 for nu, sigz, and sigRz for all populations
    chi2 = chi2_nu + chi2_sigz2 + chi2_sigRz2

    #print('P', myrank, ': chi2 = ', chi2)

    return chi2, chi2_nu_vecs, chi2_sigz2_vecs, chi2_sigRz2_vecs
## \fn calc_chi2(profs)
# Calculate chi^2
# @param profs profiles for rho, M, nu_i, beta_i, sig_i, kap_i
