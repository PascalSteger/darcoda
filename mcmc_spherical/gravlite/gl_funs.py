#!/usr/bin/python
from types import *
import pdb
import numpy.random as npr

import gl_params as gp
import gl_plot as gplot
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
    if gp.chisqt_nu > gp.chisqt_sig:
        ranarr.wiggle_nu()              # first wiggle nu, has most error on it in beginning
        ranarr.wiggle_dens()
        # ranarr.setuniformrandom() # if we want to wiggle all parameters as one
    else:
        ranarr.wiggle_dens() # if we want to wiggle only the density/mass/surfden array
        # ranarr.setuniformrandom() # if we want to wiggle all parameters as one
        # ranarr.wigglebin()       # if we want to wiggle only the first bin at a time
        # ranarr.wigglepar()       # if we want to wiggle only one parameter at a time
    ranarr.wiggle_delta()

    if not gp.logprior:
        ranarr.mul(gp.parstep)
        gp.parst = ranarr.add(gp.parst)
        if gp.deltaprior:
            if gp.investigate == 'walker':
                delta1,delta2 = betawalker(gp.xipol)
                gp.parst.set_delta([delta1, delta2])
            else:
                gp.parst.set_delta([gp.delta0,gp.delta0])

    else:
        # take one step in logspace:
        gp.lparst = gp.parst.getlog()
        gp.lparst.add(  ranarr.mul( gp.lparstep ) )
        gp.parst.nu1 = np.power( 10.0, gp.lparst.nu1  )
        gp.parst.dens = np.power( 10.0, gp.lparst.dens )
        if gp.cprior<1e29:
            gp.parst.M[0] = cprior
        if gp.deltaprior:
            if gp.investigate == 'walker':
                delta1,delta2 = betawalker(gp.xipol)
                gp.parst.set_delta([delta1, delta2])
            else:
                gp.parst.set_delta([gp.delta0,gp.delta0])
        else:
            # linear sampling for delta (since it can go negative for physical reasons):
            gp.parst.set_delta(gp.lparst.get_delta()) 


        # linear sampling for Mslope:
        gp.parst.Msl = gp.lparst.Msl
        gp.parst.set_sigsl(gp.lparst.get_sigsl())
    return gp.parst

def calc_nu_sig():
    gp.LOG.debug( 'Calculate nu / sigma values')
    if gp.geom == 'sphere':
        import physics_sphere as phys
        gp.nu1_x = phys.nu(gp.parst.nu1)
        gp.sig1_x = phys.sig_los(1)
        if (gp.pops == 2) :
            gp.nu2_x = phys.nu(gp.parst.nu2)
            gp.sig2_x = phys.sig_los(2)
        if gp.checksigma:
            gp.nu1_x = gp.ipol.nudat1
            gp.sig1_x = phys.sig_los(1)
            if gp.pops == 2:
                gp.nu2_x = gp.ipol.nudat2
                gp.sig2_x = phys.sig_los(2)
            gp.M_x   = gp.ipol.Mdat
            if gp.analytic:
                gp.nu1_x = rho_anf(gp.xipol)
                gp.M_x   = M_anf(gp.xipol)
            # gp.sig1_x= sig_los_anf(gp.xipol) # cheat to check sig_los plot
    elif gp.geom == 'disc':
        import physics_disc as phys
        # TODO: generalize for 2 populations
        # gp.nu1_x  = phys.nu_decrease(gp.xipol,gp.xipol,gp.parst.nu1)
        gp.nu1_x  = phys.nu(gp.parst.nu1)
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
            # gplot.plot(z_dat,test,psym=3)
            # gplot.plot(z_dat,prob_z,color=2,psym=3)
            # gplot.show(); exit(1)
            
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

def calc_chisqt():
    gp.LOG.info('Calculate chi-squared:')
    gp.chisqt_nu1      = ((gp.nu1_x - gp.ipol.nudat1)**2./gp.ipol.nuerr1**2.).sum()
    gp.chisqt_nu       = gp.chisqt_nu1
    if gp.analytic:
        gp.chisqt_sig1     = ((gp.sig1_x- rho_anf(gp.xipol))**2./gp.ipol.sigerr1**2.).sum() #gp.ipol.sigdat1
    else:
        gp.chisqt_sig1     = ((gp.sig1_x- gp.ipol.sigdat1)**2./gp.ipol.sigerr1**2.).sum()
    gp.chisqt_sig      = gp.chisqt_sig1

    if not gp.deltaprior and gp.uselike:
        sig_Rz = phys.sigma_rz(gp.xipol, gp.xipol, gp.parst.delta1)
        prob_t = prob_t + np.sum((sig_Rz - sigRz_dat)**2./sigRz_dat_err**2.)

    gp.chisqt1         = gp.chisqt_nu1 + gp.chisqt_sig1
    gp.chisqt          = gp.chisqt1
    if gp.pops == 2:
        gp.chisqt_nu2  = ((gp.nu2_x - gp.ipol.nudat2)**2./(gp.ipol.nuerr2)**2.).sum()
        gp.chisqt_nu  += gp.chisqt_nu2
        gp.chisqt_sig2 = ((gp.sig2_x- gp.ipol.sigdat2)**2./(gp.ipol.sigerr2)**2.).sum()
        gp.chisqt_sig += gp.chisqt_sig2
        gp.chisqt2     = gp.chisqt_nu2 + gp.chisqt_sig2
        gp.chisqt     += gp.chisqt2

    # gp.LOG.info('Calculate the f-function')
    # gp.LOG.warning(['chisq  = ',gp.chisq,'chisqt = ',gp.chisqt])
    gp.fnewoverf = np.exp(gp.chisq/2.0-gp.chisqt/2.0)
    return gp.fnewoverf

def accept_reject(n):
    gp.LOG.info( 'Accept the new parameters?')
    ran = npr.rand()
    if ran < gp.fnewoverf:
        gp.acccount = gp.acccount + 1.
        gp.pars.set(gp.parst)
        gp.chisq = gp.chisqt
        
        gfile.store_working_pars(n, gp.pars, gp.chisq, gp.parstep)
        if npr.rand() < 1./gp.nplot:
            gplot.update_plot()

        gp.LOG.warning([ n,gp.chisq,np.sum(abs(1-phys.Mzdefault(gp.pars.dens)/M_anf(gp.xipol))),\
                         gp.acccount/(gp.rejcount+1.)])
        # np.sum( abs(1-Mzdefault( gp.pars.dens )/ M_anf(gp.xipol) ) ) # to check difference from analytic mass profile
        
        # Decide whether to end initphase:
        if gp.endgame and gp.initphase == 'start':
            gp.endcount -= 1
            print 'gp.endcount = ',gp.endcount
            if (gp.endcount <= 0 or gp.chisq < gp.chisqtol/2.):
                print( '*** initialization phase over ***')
                print( '*********************************')
                gp.initphase = 'over'
                gp.pars.dens    = phys.densdefault(gp.pars.dens)
                # gp.ipol.densdat = phys.densdefault(gp.ipol.densdat)
                gp.parst.dens   = phys.densdefault(gp.parst.dens)
                # gp.parst.dens   = np.array(phys.calculate_dens(M_anf(gp.xipol),gp.xipol)) # cheat: start of from near analytic 1 pop result
                gp.parstep.dens = gp.parst.dens/20.
                gp.parstep.adaptall(1./gp.scaleafterinit)
                gp.stepafterrunaway = 1.
                gp.poly = False

    else:
        gp.rejcount = gp.rejcount + 1.
        # jump back to last known good point
        faraway = gp.farinit if gp.initphase == 'init' else gp.farover
        if gp.chisq/gp.chisqt < faraway:
            gp.LOG.warning(' too far off, setting back to last known good point')
            gfile.get_working_pars()
    return

def adapt_stepsize():
    gp.LOG.info( 'Adapt stepsize during initialisation phase: ')
    if (gp.initphase == 'start'):
        if (gp.acccount>0) and (gp.rejcount>0) :
            if (gp.acccount/gp.rejcount < gp.accrejtollow\
                or gp.acccount/gp.rejcount > gp.accrejtolhigh):
                if gp.logprior:
                    gp.lparstep.adaptworst(gp.stepcorr)
                else:
                    gp.parstep.adaptworst(gp.stepcorr)
            else:
                    gp.parstep.adaptall(1./gp.stepcorr)
            if (gp.chisq < gp.chisqtol):
                gp.endgame = True
    return
