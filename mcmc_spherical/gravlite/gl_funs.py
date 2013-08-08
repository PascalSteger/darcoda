#!/usr/bin/python
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
from types import *
import pdb
import numpy.random as npr

import gl_params as gp
if gp.geom == 'disc':
    import physics_disc as phys
else:
    import physics_sphere as phys
import gl_plot as gpl
import gl_file as gfile
from gl_class_params import *
from gl_analytic import *








def get_new_parameters():
    'Wiggle the parameters -> new parameters:'
    ranarr = Params(0)

    # change parameters in one of the following ways:
    # ranarr.wiggle_dens() # if we want to wiggle only the density/mass/surfden array
    # ranarr.wigglebin()       # if we want to wiggle only the first bin at a time
    # ranarr.wigglepar()       # if we want to wiggle only one parameter at a time
    # ranarr.wiggle_delta()

    ranarr.setuniformrandom()   # wiggle all parameters as one

    # TODO: check outcome here. was used when working with adaptstepwait=1
    # if gp.initphase:
    #     ranarr.scale_prop_chi2()    # change proportional to error on chi2
    
    ranarr.mul(gp.parstep)
    ranarr.add(gp.pars)
    gp.parst.assign(ranarr)

    # keep deltaprior during init phase, if set:
    if gp.deltaprior and gp.initphase:  # TODO: rename
        if gp.investigate == 'walker':           
            delta1,delta2 = betawalker(gp.xipol) # [1]
            gp.parst.set_delta([phys.invdelta(delta1), phys.invdelta(delta2)])
        else:
            gp.parst.set_delta([phys.invdelta(gp.delta0),phys.invdelta(gp.delta0)])

    return gp.parst













def calc_M_nu_sig_kap_sphere():
    # Calculate M / nu / sigma values for spherical case
    if gp.geom != 'sphere':
        print 'wrong branch, use calc_M_nu_sig_disc'
        pdb.set_trace()

    import physics_sphere as phys
    if not gp.checksigma:
        '''normal case'''
        gp.nu1_x  = phys.nu(gp.parst.nu1) # [munit/pc^3]
        gp.dens_x = phys.dens(gp.xipol, gp.parst.dens)         # [munit/pc^3]
        
        gp.M_x    = phys.Mr3D(gp.xipol, gp.dens_x) # [munit,3D]
        gp.d1_x   = phys.delta(gp.parst.delta1)
        gp.sig1_x = phys.sig_los(gp.M_x, gp.nu1_x, gp.d1_x)
        gp.kap1_x = phys.kap_los(gp.M_x, gp.nu1_x, gp.d1_x)

        if gp.pops == 2:
            gp.nu2_x  = phys.nu(gp.parst.nu2)
            gp.d2_x   = phys.delta(gp.parst.delta2)
            gp.sig2_x = phys.sig_los(gp.M_x, gp.nu2_x, gp.d2_x)
            gp.kap2_x = phys.kap_los(gp.M_x, gp.nu2_x, gp.d2_x)

        # to debug almost fitting curves:
        # if gp.chi2 < 60:
        #     pdb.set_trace()
        

            
    else:
        '''checksigma: check integration routine'''
        if gp.investigate == 'hernquist':
            ### set nu to data, or analytic value
            gp.nu1_x  = gp.ipol.nudat1 # [Msun/pc^3]
            if gp.pops == 2:
                gp.nu2_x = gp.ipol.nudat2 # [Msun/pc^3]
                    
            if gp.analytic:
                gp.nu1_x = rho_anf(gp.xipol)  # [Msun/pc^3]


            # set dens
            # gp.dens_x = gp.ipol.densdat # for bary only
            # gp.dens_x = phys.densdefault(gp.parst.dens) # [Munit/pc^2]
            
            gp.dens_x = rho_anf(gp.xipol) # [Msun/pc^3]
            # attention: same as gp.nu1_x in this case! fine if potential comes from 'tracers' only
            # gp.dens_x = phys.calculate_dens(gp.ipol.Mx, gp.M_x) # [munit/lunit^3]
            
            # set M
            # gp.M_x    = phys.Mr3D(gp.ipol.Mx, gp.dens_x) # [munit,3D]
            gp.M_x   = M_anf(gp.xipol)    # [Msun]
            
            # set delta
            gp.d1_x  = np.zeros(gp.nipol) # [1]
            if gp.pops == 2:
                gp.d2_x = walker_delta(2)

            gp.sig1_x = phys.sig_los(gp.M_x, gp.nu1_x, gp.d1_x)
            gp.kap1_x = phys.kap_los(gp.M_x, gp.nu1_x, gp.d1_x)
            if gp.pops == 2:
                gp.sig2_x = phys.sig_los(gp.M_x, gp.nu2_x, gp.d2_x)
                gp.kap2_x = phys.kap_los(gp.M_x, gp.nu2_x, gp.d2_x)
            
        # still checksigma, but for walker now
        elif gp.investigate == 'walker':
            rhodm, rhostar1, rhostar2 = rhowalker_3D(gp.xipol)
            gp.nu1_x = rhostar1
            if gp.pops == 2:
                gp.nu2_x = rhostar2 # gp.ipol.nudat2 # data
                    
            gp.dens_x = rhowalkertot_3D(gp.xipol) # [msun/pc^3]
            gp.M_x = phys.Mr3D(gp.ipol.Mx, gp.dens_x) # [msun]

            gp.d1_x = walker_delta(1) # [1], actual delta
            if gp.pops == 2:
                gp.d2_x = walker_delta(2) # [1]


            # set sigma_LOS, by calculating expected values
            # for both (spherical) cases simultaneously
            gp.sig1_x = phys.sig_los(gp.M_x, gp.nu1_x, gp.d1_x) # [km/s]
            # takes [munit, 3D], [munit/pc^3], [1]
            # normalization of nu does not matter, is divided out
            gp.kap1_x = phys.kap_los(gp.M_x, gp.nu1_x, gp.d1_x) # [km/s]
            
            # gp.sig1_x= sig_los_anf(gp.xipol) # cheat to check sig_los plot
            if gp.pops == 2:
                gp.sig2_x = phys.sig_los(gp.M_x,gp.nu2_x,gp.d2_x) # [km/s]
                gp.kap2_x = phys.kap_los(gp.M_x,gp.nu2_x,gp.d2_x) # [km/s]









def calc_M_nu_sig_disc():
    gp.LOG.debug( 'Calculate M / nu / sigma values for disc case')
    import physics_disc as phys

    gp.nu1_x  = phys.nu(gp.parst.nu1)
    # gp.nu1_x  = phys.nu_decrease(gp.xipol,gp.xipol,gp.parst.nu1)
    gp.d1_x   = phys.delta(gp.parst.delta1)
    gp.dens_x = phys.dens(gp.xipol, gp.parst.dens) # means Kz, *not* surface density
    # gp.dens_x = gp.parst.dens+gp.blow # Kzpars with baryonic lower limit

    gp.M_x    = phys.Sigmaz(gp.dens_x)
    # gp.M_x    = phys.kz(gp.xipol, gp.xipol, gp.dens_x, gp.blow) # from kz, added baryons min

    # TODO: marginalise over lower bound, see gl_priors, kappa_DM < 0
    # attention: gives sometime unphysical decreasing surface density like this
    gp.blow   = gp.dat.Mdat # npr.normal(gp.dat.Mdat,gp.dat.Merr/4.,gp.nipol)

    
    # TODO: naming
    Rsun = 8000.; hr = 3.0; hsig = 3.0
    gp.sig1_x = phys.sigmaz(gp.xipol, gp.dens_x, gp.nu1_x, gp.parst.norm1, gp.parst.delta1, [Rsun,hr,hsig])
    gp.kap1_x = gp.xipol*0. # set to 0 consistently

    
    if gp.pops == 2:
        gp.nu2_x  = phys.nu(gp.parst.nu2)
        # gp.nu2_x  = phys.nu_decrease(gp.xipol,gp.xipol,gp.parst.nu2)
        gp.d2_x   = phys.delta(gp.parst.delta2)
        gp.sig2_x = phys.sigmaz(gp.xipol, gp.dens_x, gp.nu2_x, gp.parst.norm2, gp.parst.delta2, [Rsun,hr,hsig])
        gp.kap2_x = gp.xipol*0.

    if gp.checksigma:
        gp.nu1_x  = gp.ipol.nudat1
        Rsun = 8.; hr = 3.0; hsig = 3.0 # irrelevant, as long as gp.deltaprior = True
        gp.parst.norm1 = 17.**2
        gp.sig1_x = phys.sigmaz(gp.xipol, gp.ipol.densdat, gp.blow, gp.ipol.nudat1,gp.parst.norm1,gp.parst.delta1, [Rsun,hr,hsig])
        gp.kap1_x = gp.xipol*0.

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
    '''return nu required for comparison to interpolated data'''
    # pop \in {1,2} gives population
    # if dat == True: return data; else: model (possibly projected)
    # if err == True: return data error instead
    if (not dat) and err:
        print 'wrong use of compare_nu, error only given for data, not model'
        exit(1)
    ret = []
    if dat:
        if gp.geom == 'disc':
            ret = gp.ipol.nudat1 if pop == 1 else gp.dat.nudat2
            if err:
                ret = gp.ipol.nuerr1 if pop == 1 else gp.ipol.nuerr2
        elif gp.geom == 'sphere':
            ret = gp.ipol.nudat1_2D if pop == 1 else gp.dat.nudat2_2D
            if err:
                ret = gp.ipol.nuerr1_2D if pop == 1 else gp.ipol.nuerr2_2D
    else:                               # [model]
        if gp.geom == 'disc':
            ret = gp.ipol.nudat1 if pop == 1 else gp.ipol.nudat2
            if err:
                ret = gp.ipol.nuerr1 if pop == 1 else gp.ipol.nuerr2
        elif gp.geom == 'sphere':
            ret = int_surfden(gp.xipol,gp.ipol.nudat1) if pop == 1 else int_surfden(gp.ipol.nudat2)
            if err:
                # TODO: does that make sense? what do we want to compare to here?
                ret = int_surfden(gp.xipol,gp.ipol.nuerr1) if pop == 1 else int_surfden(gp.ipol.nuerr2)
    return ret






def chi2red(model, data, sig, dof):
    # reduced chi2
    # if Degrees Of Freedom = 1, return non-reduced chi2
    return sum(((model-data)**2./sig**2)/dof)



def calc_chi2():
    '''Calculate new \\chi^2, stored in gp.chi2t'''
    numodel1 = int_surfden(gp.xipol,gp.nu1_x) if gp.geom=='sphere' else gp.nu1_x
    nudata1  = compare_nu(1,True,False)
    nuerr1   = compare_nu(1,True,True) # gp.ipol.nuerr1 # old, used 3D nu in spherical case
    gp.chi2t_nu1      = chi2red(numodel1, nudata1, nuerr1, gp.dof)
    gp.chi2t_nu       = gp.chi2t_nu1
    if gp.analytic:
        gp.chi2t_sig1     = chi2red(gp.sig1_x,rho_anf(gp.xipol),gp.ipol.sigerr1,1.,gp.dof) #gp.ipol.sigdat1
        # TODO: kap1 chi2 in analytic case
    else:
        gp.chi2t_sig1     = chi2red(gp.sig1_x,gp.ipol.sigdat1,gp.ipol.sigerr1,gp.dof)
        gp.chi2t_kap1     = chi2red(gp.kap1_x,gp.ipol.kapdat1,gp.ipol.kaperr1,gp.dof)
    gp.chi2t_sig      = gp.chi2t_sig1
    gp.chi2t_kap      = gp.chi2t_kap1

    if not gp.deltaprior and gp.uselike:
        sig_Rz = phys.sigma_rz(gp.xipol, gp.xipol, gp.parst.delta1)
        prob_t = prob_t + np.sum((sig_Rz - sigRz_dat)**2./sigRz_dat_err**2.)

    gp.chi2t1         = gp.chi2t_nu1 + gp.chi2t_sig1 + gp.chi2t_kap1
    gp.chi2t          = gp.chi2t1
    if gp.pops == 2:
        numodel2 = int_surfden(gp.xipol,gp.nu2_x) if gp.geom=='sphere' else gp.nu2_x
        nudata2  = compare_nu(2,True,False)
        nuerr2   = compare_nu(2,True,True)
        gp.chi2t_nu2  = chi2red(numodel2, nudata2, nuerr2, gp.dof)
        gp.chi2t_nu  += gp.chi2t_nu2
        gp.chi2t_sig2 = chi2red(gp.sig2_x, gp.ipol.sigdat2, gp.ipol.sigerr2,gp.dof)
        gp.chi2t_sig += gp.chi2t_sig2
        gp.chi2t_kap2 = chi2red(gp.kap2_x, gp.ipol.kapdat2, gp.ipol.kaperr2,gp.dof)
        gp.chi2t_kap += gp.chi2t_kap2
        gp.chi2t2     = gp.chi2t_nu2 + gp.chi2t_sig2 + gp.chi2t_kap2
        gp.chi2t     += gp.chi2t2

    if np.isnan(gp.chi2t):
        print 'NaN occured! go search where it happened!'
        pdb.set_trace()

    gp.fnewoverf = np.exp(gp.chi2/2.0-gp.chi2t/2.0)
    return gp.fnewoverf











def accept_reject(n):
    gp.LOG.info( 'Accept the new parameters?')
    if npr.rand() < gp.fnewoverf:
        gp.accrate.update(True)
        gp.pars.assign(gp.parst)
        gp.chi2 = gp.chi2t
        gp.lasterr = 'None'
        gp.d1wild = False; gp.d2wild = False; gp.dens2wild = False
        gp.b2wild = False; gp.sig2wild = False; gp.nu2wild = 1000
        gfile.store_working_pars(n, gp.pars, gp.chi2, gp.parstep)
        # fplot
        if npr.rand() < max(0.01, (1.*gp.chi2-gp.chi2tol)/gp.chi2tol/100.)\
               or (gp.initphase and gp.adaptstepwait == 1):
            gpl.update_plot()

        np.set_printoptions(precision=3)

        if gp.pops == 1:
            print 'n:',n, ' chi2:',gh.pretty(gp.chi2,1),\
                  ' rate:',gh.pretty(100*gp.accrate.rate(),2),\
                  ' nu1:',gh.pretty(100*abs(np.median((phys.nu(gp.pars.nu1+gp.parstep.nu1)\
                                                       -phys.nu(gp.pars.nu1))/\
                                                      phys.nu(gp.pars.nu1))),3),\
                  ' d1:',gh.pretty(100*abs(np.median(gp.parstep.delta1/gp.pars.delta1)),3),\
                  ' dens:',gh.pretty(100*abs(np.median((phys.densdefault(gp.parstep.dens+\
                                                                         gp.pars.dens)\
                                                        -phys.densdefault(gp.pars.dens))/\
                                                       phys.densdefault(gp.pars.dens))),3),\
                  ' norm1:',gh.pretty(100*abs(np.median(gp.parstep.norm1/gp.pars.norm1)),3)
        else:
            
            print 'n:',n, ' chi2:',gh.pretty(gp.chi2,1),\
                  ' rate:',gh.pretty(100*gp.accrate.rate(),2),\
                  ' nu1:',gh.pretty(100*abs(np.median((phys.nu(gp.pars.nu1+gp.parstep.nu1)\
                                                       -phys.nu(gp.pars.nu1))/\
                                                      phys.nu(gp.pars.nu1))),3),\
                  ' nu2:',gh.pretty(100*abs(np.median((phys.nu(gp.pars.nu2+gp.parstep.nu2)\
                                                       -phys.nu(gp.pars.nu2))/\
                                                      phys.nu(gp.pars.nu2))),3),\
                ' d1:',gh.pretty(100*abs(np.median(gp.parstep.delta1/gp.pars.delta1)),3),\
                ' d2:',gh.pretty(100*abs(np.median(gp.parstep.delta2/gp.pars.delta2)),3),\
                ' dens:',gh.pretty(100*abs(np.median((phys.densdefault(gp.parstep.dens+gp.pars.dens)\
                                                      -phys.densdefault(gp.pars.dens))/\
                                                     phys.densdefault(gp.pars.dens))),3)


        adapt_stepsize()
        end_initphase()

    else:
        gp.accrate.update(False)
        # jump back to last known good point
        faraway = gp.farinit if gp.initphase else gp.farover
        if gp.chi2t > gp.chi2 * faraway:
            gp.LOG.warning(' too far off, setting back to last known good point')
            gfile.get_working_pars(gp.initphase)
    return









def adapt_stepsize():
    '''Adapt stepsize during initialisation phase: '''
    if gp.initphase:
        if gp.adaptstepwait > 1:
            gp.adaptstepwait -= 1
            return
        print 'adapt stepsize!'
        gp.adaptstepwait = gp.rollsize
        if (not gp.accrate.rightrate()): # and (not gp.accrate.getsbetter()): # too bad, increase even if getting better. TODO: use this second part if adaptstepwait is small
            print ' >> increase'
            gp.parstep.adaptworst(gp.stepcorr)
        else:
            gp.parstep.adaptall(1./gp.stepcorr)
            print ' << decrease'
        if gp.chi2 < gp.chi2tol:
            gp.endgame = True
    return



def end_initphase():
    '''Decide whether to end initphase:'''
    if not gp.initphase: return
    if not gp.endgame: return

    if gp.chi2 < gp.chi2tol and gp.chi2t_nu < gp.chi2t_sig:
        # TODO: include kap condition?
        gp.endcount -= 1;
        print 'gp.endcount = ',gp.endcount

    if gp.chi2 < gp.chi2tol/10.:
        print 'really really small chi2! stopping init phase'
    else:       # if gp.chi2 really really small: stop NOW!
        # do not accept too high stepsizes (must be < 5% in dens)
        if abs(np.median(gp.parstep.dens/gp.pars.dens)) > 0.05:
            return
        # do not stop unless gp.endcount counted down
        if gp.endcount > 0:
            return

    # if none of the above conditions hold:
    print( '*** initialization phase over ***')
    print( '*********************************')
    gp.initphase = False
    if gp.denslog:
        # use parstep if possible
        # gp.parstep.dens = abs(np.log10(phys.densdefault(gp.pars.dens+gp.parstep.dens))-\
        #                       np.log10(phys.densdefault(gp.pars.dens)))
        # or get 10% step
        # gp.parstep.dens = abs(np.log10(phys.densdefault(gp.pars.dens*1.05))-\
        #                       np.log10(phys.densdefault(gp.pars.dens)))

        # 2) use constant step, from largest step, factor-wise
        mul = np.median(gp.parstep.dens[0]/gp.pars.dens[0])
        # # gp.parstep.dens = mul*np.log10(phys.densdefault(gp.pars.dens)) # not proportional
        # # gp.parstep.dens = abs(np.log10(phys.densdefault(gp.pars.dens))*mul) # prop, wrong for arg = 0
        gp.parstep.dens = np.ones(gp.nipol)*abs(np.log10(\
            phys.densdefault(gp.pars.dens)[-1])*mul)
        gp.pars.dens = np.log10(phys.densdefault(gp.pars.dens))
        gp.parst.dens = np.log10(phys.densdefault(gp.parst.dens))
        # constant stepsize over all xipol, from highest contrib from
        # last point, same sign everywhere (even if first point is log10>0.)
        # end 2)


        # 3) use real space density change
        # den        = phys.densdefault(gp.pars.dens)
        # denerr     = phys.densdefault(gp.pars.dens+gp.parstep.dens/2)-phys.densdefault(gp.pars.dens)
        # gp.parstep.dens = np.array(np.log10(den+denerr)-np.log10(den))
        # boundbelow = [max(x, max(gp.parstep.dens)/10) for x in gp.parstep.dens] # give last bins a chance
        # gp.parstep.dens = np.array(boundbelow)
        # gp.pars.dens    = np.array(np.log10(den))
        # gp.parst.dens = np.array(np.log10(phys.densdefault(gp.parst.dens)))
        # end 3)
    else:
        gp.parstep.dens = abs(phys.densdefault(gp.pars.dens+gp.parstep.dens)-\
                              phys.densdefault(gp.pars.dens))
        gp.pars.dens    = phys.densdefault(gp.pars.dens)
        gp.parst.dens   = phys.densdefault(gp.parst.dens)
        
    # gp.parst.dens = np.array(phys.calculate_dens(gp.xipol,M_anf(gp.xipol)))
    # ^-- cheat: start off, from near analytic 1 pop result
    # gp.parstep.dens = gp.parst.dens/20. 
    gp.parstep.adaptall(gp.scaleafterinit)
    gp.stepafterrunaway = 1.
    gp.poly = False
    gp.safepars.assign(gp.pars)
    gp.safeparstep.assign(gp.parstep)
    gp.safechi2 = gp.chi2

    # create folder to hold all output files
    import os; import os.path
    if not os.path.exists(gp.files.outdir):
        os.makedirs(gp.files.outdir)
        os.makedirs(gp.files.outdir+'/programs/')

    # copy across all programs to get working copy for easy repetition
    import os
    os.system('cp '+ gp.files.progdir+'/*.py '+gp.files.outdir+'/programs/')

    gfile.adump()                       # only do after endinit
    return
