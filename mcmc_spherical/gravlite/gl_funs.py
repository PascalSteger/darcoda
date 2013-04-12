#!/usr/bin/python2.7
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
from types import *
import pdb
import numpy.random as npr

import gl_params as gp
import gl_plot as gpl
import gl_file as gfile
from gl_class_params import *
if gp.geom == 'disc':
    import physics_disc as phys
elif gp.geom == 'sphere':
    import physics_sphere as phys
from gl_analytic import *








def get_new_parameters():
    'Wiggle the parameters -> new parameters:'
    ranarr = Params(0)

    # change parameters in one of the following ways:
    # ranarr.setuniformrandom() # if we want to wiggle all parameters as one
    # ranarr.wiggle_dens() # if we want to wiggle only the density/mass/surfden array
    # ranarr.wigglebin()       # if we want to wiggle only the first bin at a time
    # ranarr.wigglepar()       # if we want to wiggle only one parameter at a time
    # ranarr.wiggle_delta()


    ranarr.setuniformrandom()
    ranarr.scale_prop_chi2()

    if not gp.logprior:
        ranarr.mul(gp.parstep)
        gp.parst = ranarr.add(gp.parst)

    else:
        # take one step in logspace:
        gp.lparst = gp.parst.getlog()
        gp.lparst.add(  ranarr.mul( gp.lparstep ) )
        gp.parst.nu1 = np.power( 10.0, gp.lparst.nu1  )
        gp.parst.dens = np.power( 10.0, gp.lparst.dens )
        if gp.cprior<1e29:
            gp.parst.M[0] = cprior
        else:
            # linear sampling for delta (since it can go negative for physical reasons):
            gp.parst.set_delta(gp.lparst.extra_delta()) 

        # linear sampling for Mslope:
        gp.parst.Msl = gp.lparst.Msl    # [Msun/pc]
        gp.parst.set_sigsl(gp.lparst.get_sigsl())



    if gp.deltaprior:
        if gp.investigate == 'walker':
            delta1,delta2 = betawalker(gp.xipol) # [1]
            gp.parst.set_delta([phys.invdelta(delta1), phys.invdelta(delta2)])
        else:
            gp.parst.set_delta([phys.invdelta(gp.delta0),phys.invdelta(gp.delta0)])



    return gp.parst













def calc_M_nu_sig():
    gp.LOG.debug( 'Calculate M / nu / sigma values')
    if gp.geom == 'sphere':
        import physics_sphere as phys
        if not gp.checksigma:
            '''normal case'''
            gp.nu1_x  = phys.nu(gp.parst.nu1) # [munit/pc^3]
            gp.dens_x = phys.dens(gp.xipol, gp.parst.dens)         # [munit/pc^3]

            gp.M_x    = phys.Mr3D(gp.xipol, gp.dens_x) # [munit,3D]
            gp.d1_x   = phys.delta(gp.parst.delta1)
            gp.sig1_x = phys.sig_los(1, gp.xipol, gp.M_x, gp.nu1_x, gp.d1_x)
            
            if gp.pops == 2:
                gp.nu2_x  = phys.nu(gp.parst.nu2)
                gp.d2_x   = phys.delta(gp.parst.delta2)
                gp.sig2_x = phys.sig_los(2, gp.xipol, gp.M_x, gp.nu2_x, gp.d2_x)





        else:
            '''check integration routine'''
            gp.nu1_x  = gp.ipol.nudat1 # [dens0]
            if gp.investigate == 'hernquist' and gp.analytic:
                gp.nu1_x = rho_anf(gp.xipol) # [TODO]
                gp.M_x   = M_anf(gp.xipol)



            # gp.dens_x = gp.ipol.densdat # for bary only
            # gp.dens_x = phys.densdefault(gp.parst.dens) # [Munit/pc^2]

            gp.dens_x = rhowalkertot_3D(gp.xipol) # [msun/pc^3]

            # tmp = rhowalker_3D(gp.xipol)
            # gp.dens_x = tmp[1]+tmp[2]          # theoretical baryonic contribution only
            gp.M_x    = phys.Mr3D(gp.ipol.Mx, gp.dens_x) # [munit,3D]
            
            # gp.totmass[0] = max(gp.M_x) # renormalization, if necessary


            # gp.M_x    = gp.ipol.Mdat   # [totmass] baryonic only
            gp.dens_x = phys.calculate_dens(gp.ipol.Mx, gp.M_x) # [munit/lunit^3]


            delta_r = walker_delta(1)   # [1], actual delta
            
            gp.sig1_x = phys.sig_los(1,gp.ipol.sigx1, gp.M_x, gp.nu1_x, delta_r) # [km/s]
            # takes [1], [pc], [munit, 3D], [munit/pc^3], [1]
            # normalization of nu does not matter, is divided out

            
            if gp.pops == 2:
                gp.nu2_x = gp.ipol.nudat2
                delta_r = walker_delta(2)
                gp.sig2_x = phys.sig_los(2,gp.ipol.sigx2,gp.M_x,gp.nu2_x,delta_r) # [km/s]


            # gp.sig1_x= sig_los_anf(gp.xipol) # cheat to check sig_los plot
            # sig1  = gp.ipol.sigdat1 # to check: plot is really done as it should






    elif gp.geom == 'disc':
        import physics_disc as phys
        # TODO: generalize for 2 populations
        # gp.nu1_x  = phys.nu_decrease(gp.xipol,gp.xipol,gp.parst.nu1)
        gp.nu1_x  = phys.nu(gp.parst.nu1)
        gp.dens_x = gp.parst.dens # Kzpars
            # TODO: naming
        Rsun = 8.; hr = 3.0; hsig = 3.0
        gp.sig1_x = phys.sigmaz(gp.xipol,gp.xipol,gp.parst.dens,gp.parst.norm,\
                                    gp.blow,gp.parst.nu1,gp.parst.delta1,\
                                    [Rsun,hr,hsig])
        if gp.checksigma:
            gp.nu1_x  = gp.ipol.nudat1
            Rsun = 8.; hr = 3.0; hsig = 3.0 # irrelevant, as long as gp.deltaprior = True
            gp.parst.norm = 21.**2
            gp.sig1_x = phys.sigmaz(gp.xipol,gp.xipol,gp.ipol.densdat,  gp.blow, gp.ipol.nudat1,gp.parst.norm,gp.parst.delta1,[Rsun,hr,hsig])

    return





def calc_likelihood():
    if gp.uselike:
        # Errors:
        if gp.adderrors:
            # Solve error convolution: # can be jumped
            prob_z = np.zeros(len(z_dat))
            for jj in range(len(z_dat)):
                zintmin = -abs(z_dat_err[jj])*3.0
                zintmax = abs(z_dat_err[jj])*3.0
                zintpnts = 25
                zint = np.arange(zintpnts) * (zintmax-zintmin)/(zintpnts-1.) + zintmin + z_dat[jj]
                perr_z = errorz(z_dat[jj] - zint, z_dat_err[jj])
                nu_int = nu(abs(zint),zp_kz,nupars)
                p_z = nu_int
                prob_z[jj] = simps(perr_z * p_z,zint)
                if prob_z[jj] < 0:
                    print 'Ooops a negative probability! Something is wrong ... stopping.'
                    exit(1)
            nu_pz = phys.nu(abs(z_dat),zp_kz,nupars)
            prob_z = nu_pz
            
            # FOR TESTING! 
            # test = gp.nu1_x
            # gpl.plot(z_dat,test,psym=3)
            # gpl.plot(z_dat,prob_z,color=2,psym=3)
            # gpl.show(); exit(1)
            
            sig_sum = sqrt(sig_z**2. + vz_dat_err**2.)
            aprob_sigz = alog(1.0/(sqrt(2.0*np.pi) * sig_sum)) - (vz_dat-vz_mean)**2./2./sig_sum**2.
            
            # *** WARNING *** TI< with errors not yet supported ... ! 
            prob_tilt = 1.0
        else:
            # Calcualte likelihood [N.B. Log[Li] can be +ve!]:
            prob_z = nu_z
            aprob_sigz = np.log(1.0/(sqrt(2.0*np.pi) * sig_z)) - (vz_dat-vz_mean)**2./2./sig_z**2.
            
            if not gp.deltaprior:
                wid_Rz = sigma_rz(abs(z_dat),zp_kz,tparswt)
                prob_tilt = 1.0/np.pi/wid_Rz * BESELK(abs(vz_dat*vR_dat-sig_Rz)/wid_Rz,0)
            else:
                prob_tilt = 1.0

        prob_t = np.sum(np.log(prob_z) + aprob_sigz + np.log(prob_tilt))
        if math.isnan(prob_t):
            print 'Ooops prob_t is NaN ... '
            exit(1)

        # Calculate the likelihood ratio:
        fnewoverf = np.exp(prob_t - prob)
    return gp.fnewoverf






def compare_nu(pop, dat, err):
    '''return nu required for comparison to data'''
    # pop \in {1,2} gives population
    # if dat == True: return data; else: model (possibly projected)
    # if err == True: return data error instead
    if not dat and err:
        print 'wrong use of compare_nu, error only given for data, not model'
        exit(1)
    ret = []
    if dat:
        if gp.geom == 'disc':
            ret = gp.dat.nudat1 if pop == 1 else gp.dat.nudat2
            if err:
                ret = gp.dat.nuerr1 if pop == 1 else gp.dat.nuerr2
        elif gp.geom == 'sphere':
            ret = gp.dat.nudat1_2D if pop == 1 else gp.dat.nudat2_2D
            if err:
                ret = gp.dat.nuerr1_2D if pop == 1 else gp.dat.nuerr2_2D
    else:
        if gp.geom == 'disc':
            ret = gp.ipol.nudat1 if pop == 1 else gp.ipol.nudat2
            if err:
                ret = gp.ipol.nuerr1 if pop == 1 else gp.ipol.nuerr2
        elif gp.geom == 'sphere':
            ret = int_surfden(gp.ipol.nudat1) if pop == 1 else int_surfden(gp.ipol.nudat2)
            if err:
                # TODO: does that make sense? what do we want to compare to here?
                ret = int_surfden(gp.ipol.nuerr1) if pop == 1 else int_surfden(gp.ipol.nuerr2)
    return ret




def calc_chi2():
    '''Calculate chi squared:'''

    numodel1 = int_surfden(gp.nipol,gp.xipol,gp.nu1_x) if gp.geom=='sphere' else gp.nu1_x
    nudata1  = compare_nu(1,True,False)
    nuerr1   = compare_nu(1,True,True) # gp.ipol.nuerr1 # old, used 3D nu in spherical case
    gp.chi2t_nu1      = ((numodel1 - nudata1)**2./nuerr1**2.).sum()
    gp.chi2t_nu       = gp.chi2t_nu1
    if gp.analytic:
        gp.chi2t_sig1     = ((gp.sig1_x- rho_anf(gp.xipol))**2./gp.ipol.sigerr1**2.).sum() #gp.ipol.sigdat1
    else:
        gp.chi2t_sig1     = ((gp.sig1_x- gp.ipol.sigdat1)**2./gp.ipol.sigerr1**2.).sum()
    gp.chi2t_sig      = gp.chi2t_sig1

    if not gp.deltaprior and gp.uselike:
        sig_Rz = phys.sigma_rz(gp.xipol, gp.xipol, gp.parst.delta1)
        prob_t = prob_t + np.sum((sig_Rz - sigRz_dat)**2./sigRz_dat_err**2.)

    gp.chi2t1         = gp.chi2t_nu1 + gp.chi2t_sig1
    gp.chi2t          = gp.chi2t1
    if gp.pops == 2:
        numodel2 = int_surfden(gp.nipol,gp.xipol,gp.nu2_x) if gp.geom=='sphere' else gp.nu2_x
        nudata2  = compare_nu(2,True,False)
        nuerr2   = compare_nu(2,True,True)
        gp.chi2t_nu2  = ((numodel2 - nudata2)**2./nuerr2**2.).sum()
        gp.chi2t_nu  += gp.chi2t_nu2
        gp.chi2t_sig2 = ((gp.sig2_x- gp.ipol.sigdat2)**2/gp.ipol.sigerr2**2.).sum()
        gp.chi2t_sig += gp.chi2t_sig2
        gp.chi2t2     = gp.chi2t_nu2 + gp.chi2t_sig2
        gp.chi2t     += gp.chi2t2

    # gp.LOG.info('Calculate the f-function')
    # gp.LOG.warning(['chi2  = ',gp.chi2,'chi2t = ',gp.chi2t])
    gp.fnewoverf = np.exp(gp.chi2/2.0-gp.chi2t/2.0)
    return gp.fnewoverf











def accept_reject(n):
    gp.LOG.info( 'Accept the new parameters?')
    ran = npr.rand()
    if ran < gp.fnewoverf:
        gp.acccount = gp.acccount + 1.
        gp.pars.set(gp.parst)
        gp.chi2 = gp.chi2t
        
        gfile.store_working_pars(n, gp.pars, gp.chi2, gp.parstep)
        if npr.rand() < gp.fplot:
            gpl.update_plot()

            
        gp.LOG.warning([ n,gh.pretty([gp.chi2]),\
#                         np.sum(abs(1-gp.pars.dens/rhowalkertot_3D(gp.xipol))),\
                         gh.pretty([gp.acccount/(gp.rejcount+1.)])])
        # to check difference from analytic mass profile
        
        # Decide whether to end initphase:
        if gp.endgame and gp.initphase:
            gp.endcount -= 1
            print 'gp.endcount = ',gp.endcount
            if (gp.endcount <= 0 or gp.chi2 < gp.chi2tol/2.):
                print( '*** initialization phase over ***')
                print( '*********************************')
                gp.initphase = False
                gp.pars.dens    = phys.densdefault(gp.pars.dens)
                # gp.ipol.densdat = phys.densdefault(gp.ipol.densdat)
                gp.parst.dens   = phys.densdefault(gp.parst.dens)
                # gp.parst.dens   = np.array(phys.calculate_dens(gp.xipol,M_anf(gp.xipol)))
                # ^-- cheat: start off, from near analytic 1 pop result
                gp.parstep.dens = gp.parst.dens/20.
                gp.parstep.adaptall(1./gp.scaleafterinit)
                gp.stepafterrunaway = 1.
                gp.poly = False

    else:
        gp.rejcount = gp.rejcount + 1.
        # jump back to last known good point
        faraway = gp.farinit if gp.initphase else gp.farover
        if gp.chi2/gp.chi2t < faraway:
            gp.LOG.warning(' too far off, setting back to last known good point')
            gfile.get_working_pars()
    return









def adapt_stepsize():
    gp.LOG.info( 'Adapt stepsize during initialisation phase: ')
    if gp.initphase:
        if gp.acccount>0 and gp.rejcount>0:
            if (gp.acccount/gp.rejcount < gp.accrejtollow\
                or gp.acccount/gp.rejcount > gp.accrejtolhigh):
                if gp.logprior:
                    gp.lparstep.adaptworst(gp.stepcorr)
                else:
                    gp.parstep.adaptworst(gp.stepcorr)
            else:
                    gp.parstep.adaptall(1./gp.stepcorr)
            if gp.chi2 < gp.chi2tol:
                gp.endgame = True
    return
