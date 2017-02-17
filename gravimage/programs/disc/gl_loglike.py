## SS write new loglike file which only has the actually used settings in it
# (25 april -16), due to lack of overwiev.
# I.e. disables a lot of parameter choices in gl_params

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
    tmp_profs = Profiles(gp.ntracer_pops, gp.nbins, gp.no_Sigrho_bins)
    if gp.print_flag: print ('In geom_loglike in gl_loglike start -----------------------')
    #print ('bin centers:',gp.z_bincenter_vecs)
    #print ('extended bin centers:',gp.extend_z_binc_vecs)

    # FIRST READ IN PARAMETER VALUES:
    off = 0
    # Dark matter:
    offstep = gp.nDM_params
    DM_params = np.array(cube[off:off+offstep])
    rho_const_DM = DM_params[0]
    off += offstep

    # Baryons:
    offstep = gp.nbaryon_params
    baryon_params = np.array(cube[off:off+offstep])
    off += offstep

    # Tracers:
    tracer_params=[[].append(None) for ii in range(0,gp.ntracer_pops)]
    for pop in range(0, gp.ntracer_pops):
        offstep = gp.N_nu_model_exps*2
        tracer_params[pop] = np.array(cube[off:off+offstep])
        #print ('pop:',pop,' tracer_params[pop]',tracer_params[pop])
        off += offstep

    # Tilt term
    sigmaRz2_vecs = [[].append(None) for ii in range(0,gp.ntracer_pops)]
    tilt_vecs     = [[].append(None) for ii in range(0,gp.ntracer_pops)]
    tilt_params=[[].append(None) for ii in range(0,gp.ntracer_pops)]
    if gp.tilt:
        for t_pop in range(0, gp.ntracer_pops):
            offstep = gp.ntilt_params
            tilt_params[t_pop] = np.array(cube[off:off+offstep])
            off += offstep
    if gp.analytic_sigz2 == False:
        #Normalisation constant C for sigz calculation
        offstep = gp.ntracer_pops
        norm_C = cube[off:off+offstep]  # referred to as IntC in gl_class_cube(?)
        off += offstep

    if off != gp.ndim:
        gh.LOG(1,'wrong subscripts in gl_class_cube')
        print ('in loglike',off,gp.ndim)
        raise Exception('wrong subscripts in gl_class_cube')

    # Export Sig and rho profiles (at gp.z_vec_Sigrho locations) for future plotting
    if gp.darkmattermodel == 'const_dm':
        rho_DM_vec = rho_const_DM * np.ones(len(gp.z_vec_Sigrho)) 
        Sig_DM_vec = 2.* rho_const_DM * gp.z_vec_Sigrho
    elif gp.darkmattermodel == 'ConstPlusDD':
        rho_DM_vec = phys.rho_dm_simplenu(gp.z_vec_Sigrho, DM_params)
        Sig_DM_vec = phys.Sig_dm_simplenu(gp.z_vec_Sigrho, DM_params)
    # Code assumes only one baryon pop:
    if gp.baryonmodel == 'simplenu_baryon': 
        rho_baryon_vec = phys.rho_baryon_simplenu(gp.z_vec_Sigrho, baryon_params)
        Sig_baryon_vec = phys.Sig_baryon_simplenu(gp.z_vec_Sigrho, baryon_params)
    elif gp.baryonmodel == 'obs_baryon': 
        rho_baryon_vec = phys.rho_baryon_obs(gp.z_vec_Sigrho, baryon_params)
        Sig_baryon_vec = phys.Sig_baryon_obs(gp.z_vec_Sigrho, baryon_params)
    elif gp.baryonmodel == 'trivial_baryon':
        rho_baryon_vec = np.zeros(len(gp.z_vec_Sigrho))
        Sig_baryon_vec = baryon_params[0]*np.ones(len(gp.z_vec_Sigrho))
    rho_total_vec = rho_DM_vec + rho_baryon_vec
    Sig_total_vec = Sig_DM_vec + Sig_baryon_vec

    #print ('rho_DM_vec:',rho_total_vec)
    #print ('rho_baryon_vec:',rho_total_vec)
    #print ('rho_total_vec:',rho_total_vec)
    #print ('Sig_DM_vec:',Sig_total_vec)
    #print ('Sig_baryon_vec:',Sig_total_vec)
    #print ('Sig_total_vec:',Sig_total_vec)

    pop = 0  # pop not used in set_prof for rho and Sig stuff
    tmp_profs.set_prof('z_vec_Sigrho', gp.z_vec_Sigrho, pop, gp)
    tmp_profs.set_prof('rho_DM_vec', rho_DM_vec, pop, gp)
    tmp_profs.set_prof('rho_baryon_vec', rho_baryon_vec, pop, gp)
    tmp_profs.set_prof('rho_total_vec', rho_total_vec, pop, gp)
    tmp_profs.set_prof('Sig_DM_vec', Sig_DM_vec, pop, gp)
    tmp_profs.set_prof('Sig_baryon_vec', Sig_baryon_vec, pop, gp)
    tmp_profs.set_prof('Sig_total_vec', Sig_total_vec, pop, gp)


    # Export the population dependent profiles
    for pop in range(0, gp.ntracer_pops):
        tmp_profs.set_prof('z_vecs', gp.z_bincenter_vecs[pop], pop, gp)
        tmp_profs.set_prof('nu_z_vecs', gp.nu_z_bincenter_vecs[pop], pop, gp)
        #print ('nu_z_vec:',gp.nu_z_bincenter_vecs[pop])

        # nu vector at z_binc locations and at nu_z_binc loc
        tmp_nu_of_z = np.zeros(len(gp.z_bincenter_vecs[pop]))
        tmp_nu_of_nu_z = np.zeros(len(gp.nu_z_bincenter_vecs[pop]))
        for jter in range(0, 2*gp.N_nu_model_exps,2):
            #print ('pop:',pop,' A:',tracer_params[pop][jter],' h:',tracer_params[pop][jter+1])
            #tmp_nu_of_z +=  tracer_params[pop][jter] * np.exp(-gp.z_bincenter_vecs[pop]/tracer_params[pop][jter+1])
            #tmp_nu_of_nu_z +=  tracer_params[pop][jter] * np.exp(-gp.nu_z_bincenter_vecs[pop]/tracer_params[pop][jter+1])
            z_vec_temp = gp.z_bincenter_vecs[pop] - gp.nu_z_bincenter_vecs[pop][0]
            tmp_nu_of_z    +=  tracer_params[pop][jter] * np.exp(-z_vec_temp/tracer_params[pop][jter+1])
            z_vec_temp = gp.nu_z_bincenter_vecs[pop] - gp.nu_z_bincenter_vecs[pop][0]
            tmp_nu_of_nu_z +=  tracer_params[pop][jter] * np.exp(-z_vec_temp/tracer_params[pop][jter+1])

        #print ('tmp_nu_of_z:',tmp_nu_of_z)
        #print ('tmp_nu_of_nu_z:',tmp_nu_of_nu_z)

        tmp_profs.set_prof('nu_vecs', tmp_nu_of_nu_z, pop, gp) 
        # store nu_vecs given at the nu data locations

        #Calculate tilt for each population 
#        z_points_tmp = np.append(0., gp.z_bincenter_vecs[t_pop])
        if gp.tilt:
            A_tilt = tilt_params[pop][0]
            n_tilt = tilt_params[pop][1]
            #R_tilt = tilt_params[pop][2]
            k_tilt = tilt_params[pop][2]

            sigmaRz2_vecs[pop] = A_tilt*(gp.z_bincenter_vecs[pop])**n_tilt #all z for this pop
            tmp_profs.set_prof('sigmaRz2_vecs', sigmaRz2_vecs[pop], pop, gp)
            #tilt_vecs[pop] = sigmaRz2_vecs[pop]*(1./gp.Rsun - 2./R_tilt)
            tilt_vecs[pop] = sigmaRz2_vecs[pop]*(1./gp.Rsun - k_tilt)
        else:
            tilt_vecs[pop] = np.zeros(gp.nbins[pop][1])
        tmp_profs.set_prof('tilt_vecs', tilt_vecs[pop], pop, gp)

        #Calculate sigma (velocity dispersion)
        if gp.print_flag:  print ('Before sigz2_vec: pop=',pop)
        # Need the Sig values at the locations relevant for current sigz calc:
        if gp.darkmattermodel == 'const_dm':
            Sig_DM_vec = 2.* rho_const_DM * gp.z_bincenter_vecs[pop]
        elif gp.darkmattermodel == 'ConstPlusDD':
            Sig_DM_vec = phys.Sig_dm_simplenu(gp.z_bincenter_vecs[pop], DM_params)

        if gp.baryonmodel == 'simplenu_baryon':
            Sig_baryon_vec = phys.Sig_baryon_simplenu(gp.z_bincenter_vecs[pop], baryon_params)
        elif gp.baryonmodel == 'obs_baryon':
            Sig_baryon_vec = phys.Sig_baryon_obs(gp.z_bincenter_vecs[pop], baryon_params)
        elif gp.baryonmodel == 'trivial_baryon':
            Sig_baryon_vec = baryon_params[0]*np.ones(len(gp.z_bincenter_vecs[pop]))


        # Only for debugging, remove block below !!! 
        #print ('C:',norm_C[pop])
        #sigz2_vec = phys.sigz2(gp.z_bincenter_vecs[pop], Sig_DM_vec+Sig_baryon_vec, tilt_vecs[pop], tmp_nu_of_z, norm_C[pop])
        #print ('sigz2_vec:',sigz2_vec)    
   

        if gp.analytic_sigz2:
            sigz2_vec = phys.sigz2_obsbary(gp.z_bincenter_vecs[pop], baryon_params, tracer_params[pop][1], rho_const_DM)
            if gp.tilt:
                #A_tilde = A_tilt*(1./gp.Rsun - 2./R_tilt)
                A_tilde = A_tilt*(1./gp.Rsun - k_tilt)
                sigz2_tilt_vec = phys.sigz2_tilt(gp.z_bincenter_vecs[pop], tracer_params[pop][1], A_tilde, n_tilt) 
                #print ('sigz2_vec:',sigz2_vec)
                #print ('sigz2_tilt_vec:',sigz2_tilt_vec)
                #print ('----------------------------')
                sigz2_vec += sigz2_tilt_vec
        else:
            try:  # np.zeros(len(gp.extend_z_binc_vecs[pop]))
                sigz2_vec = phys.sigz2(gp.z_bincenter_vecs[pop], Sig_DM_vec+Sig_baryon_vec, tilt_vecs[pop], tmp_nu_of_z, norm_C[pop])
            except ValueError:
                if gp.print_flag : 
                    print ('Negative value in sigz2 array (back in gl_l..), pop=',pop)
                #print('norm_C[pop] = ', norm_C[pop])
            #raise ValueError('negative value in sig2 array')
                    raise ValueError('negative value in, or falling, sig2 array')
                return

        if gp.fit_to_sigz2 == False:
            sigz2_vec = np.sqrt(sigz2_vec)  # Store sigz (instead of sigz2)
            if gp.print_flag: print ('After:pop=',pop)
        tmp_profs.set_prof('sigz2_vecs', sigz2_vec, pop, gp)

        if (np.isinf(sigz2_vec)).any():
            print ('Infinite sigz2_vec,  pop=',pop)
            print ('sigz2_vec',sigz2_vec)
            print ('Tracer params:',tracer_params[pop])
            print ('tmp_nu_of_z:',tmp_nu_of_z)
            print ('Sig_baryon_vec:',Sig_baryon_vec)
            print ('Sig_DM_vec:',Sig_DM_vec)
            print ('tilt_vecs:',tilt_vecs[pop][1:])
            print ('gp.z_bincenter_vecs:',gp.z_bincenter_vecs[pop])
            pdb.set_trace()

    if gp.print_flag: print ('Before entering chi2 calculation (calc_chi2)')
    # determine log likelihood
    chi2, chi2_nu_vecs_tmp, chi2_sigz2_vecs_tmp, chi2_sigRz2_vecs_tmp = calc_chi2(tmp_profs, gp)
    if gp.print_flag: print ('chi2 in gl_loglike:',chi2)
    for pop in range(0, gp.ntracer_pops):
        tmp_profs.set_prof('chi2_nu_vecs', chi2_nu_vecs_tmp[pop], pop, gp)
        tmp_profs.set_prof('chi2_sigz2_vecs', chi2_sigz2_vecs_tmp[pop], pop, gp)
        tmp_profs.set_prof('chi2_sigRz2_vecs', chi2_sigRz2_vecs_tmp[pop], pop, gp)

    gh.LOG(1, '   log L = ', -chi2/2.)
    tmp_profs.chi2 = chi2

    if gp.print_flag: print ('THE END, chi2:',chi2,' |||||||||||||||||||||||')

    return tmp_profs   # from   likelihood L = exp(-\chi^2/2), want log of that
## \fn geom_loglike(cube, ndim, nparams, gp)
# define log likelihood function to be called by pyMultinest and plot_profiles
# disc version
# @param cube parameter cube as defined by gl_class_cube, in physical space
# @param ndim number of dimensions, needed as argument by multinest
# @param nparams number of parameters, needed as argument by multinest
# @param gp global parameters
