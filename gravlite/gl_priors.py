#!/usr/bin/env ipython

##
# @file
# check all parameters for prior constraints
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import pdb
import numpy as np
from scipy.interpolate import splrep, splev

import gl_params as gp
if gp.geom == 'disc':
    import physics_disc as phys
else:
    import physics_sphere as phys
import gl_plot as gpl
import gl_file as gf
from gl_project import rho_INT_Rho
from gl_int import g

## check that density is not jumping wildly
def check_density():
    if gp.rprior:
        nu1 = phys.nu(gp.parst.nu1)
        rightnu = nu1[1:];        leftnu  = nu1[:-1]
        if sum(rightnu/leftnu > gp.nutol)>0:
            # print('nutol = ',gp.nutol,' < nu1: ',max(rightnu/leftnu))
            if gp.lasterr == 'nu': gp.nu2wild -= 1
            if gp.nu2wild <= 0:
                print('*** nutol_1 runaway condition detected,\
                setting parameters back to start values')
                # gp.pars.assign(gp.safepars); gp.parstep.assign(gp.safeparstep)
                # gp.chi2 = gp.safechi2;
                # TODO: check which one of nu1 or nu2 is correct handling of runaway
                # TODO: check possible handling of exception
                gp.pars.nu1 = gp.safepars.nu1
                gp.parstep.nu1 = gp.safeparstep.nu1
                gp.pars.dens = gp.safepars.dens
                gp.parstep.dens = gp.safeparstep.dens
                gp.nu2wild = 1000
            else:
                gf.get_working_pars() # change gp.pars to gp.init_configs
            gp.lasterr = 'nu'
            return True
        if gp.pops==2:
            nu2 = phys.nu(gp.parst.nu2)
            rightnu = nu2[1:];           leftnu  = nu2[:-1]
            if sum(rightnu/leftnu > gp.nutol)>0:
                # print('nutol = ',gp.nutol,' < nu2: ',max(rightnu/leftnu))
                if gp.lasterr == 'nu': gp.nu2wild -= 1
                if gp.nu2wild <= 0:
                    gp.pars.assign(gp.safepars); gp.parstep.assign(gp.safeparstep)
                    gp.chi2 = gp.safechi2;       gp.nu2wild = 1000
                else:
                    gf.get_working_pars(False)    # change gp.parst to gp.init_configs
                gp.lasterr = 'nu'    
                return True

    # now checking dens > gprior
    if gp.gprior > 0 and max(denscheck) > gp.gpriorconv:
        print('gprior')
        return True

    return False


## check that mass is not jumping wildly
# TODO: work with exceptions
def check_mass():
    if gp.geom == 'sphere':
        # check rising mass prior
        denscheck = phys.dens(gp.xipol, gp.parst.dens)

        # if max((denscheck[1:]-denscheck[:-1])/denscheck[:-1])>0.5:
        for i in range(len(denscheck)-1):
            if (denscheck[i+1]-denscheck[i])/denscheck[i] > gp.ktol:
                if not gp.dens2wild:
                    print('rising dens prior')
                    gp.dens2wild = True
                    gf.get_working_pars()
                # gp.parst.dens[i+1] *= 0.9
                # gp.parst.dens *= 1./np.sqrt(np.arange(1.,gp.nipol+1)[::-1])
                return True

        # check rho \propto r^(-\alpha), \alpha>0.:
        # first, wanted to enforce alpha>3, 
        # but does not have to be at these small scales
        lre = np.log(gp.xipol[-1]);      lrm = np.log(gp.xipol[-2])
        lrhoe = np.log(denscheck[-1]);   lrhom = np.log(denscheck[-2])
        if (lrhoe-lrhom)/(lre-lrm) > 0:
            return True

    elif gp.geom == 'disc':
        denscheck = phys.Sigmaz(phys.densdefault(gp.parst.dens))
        # TODO: prior as described in paper
        if min(denscheck[1:]-denscheck[:-1])<0.0:
            if not gp.dens2wild:
                print('Surface density decrease found')
                gp.dens2wild = True
            gf.get_working_pars()
            return True

    gp.LOG.debug('check that observed tracer mass is less than total mass')
    if gp.bprior:
        for jj in range(gp.nipol):
            if denscheck[jj] < gp.blow[jj]:
                if not gp.b2wild:
                    print('bprior')
                gp.b2wild = True
                return True
            
    if gp.geom == 'disc':
        if min(gp.parst.dens < 0.):
            print('kappa < 0')
            gf.get_working_pars()
            return True
        kappa_DM = gp.parst.dens - phys.kappa(gp.xipol, -gp.blow*2.*np.pi*gp.G1)
        if min(kappa_DM) < 0.:
            print('kappa_DM < 0')
            gf.get_working_pars()
            return True
        if max(abs((kappa_DM-np.mean(kappa_DM))/np.mean(kappa_DM)))>1.0:
            print('kappa_DM too wild')
            # pdb.set_trace()
            gf.get_working_pars()
            return True

    # for i in range(1,gp.nipol-1):
    #     alphalong = (kappa_DM[i+1]-kappa_DM[i-1])/(gp.xipol[i+1]-gp.xipol[i-1])
    #     alphashort = (kappa_DM[i]-kappa_DM[i-1])/(gp.xipol[i]-gp.xipol[i-1])
    #     if alphashort > 0.001 and alphashort > 2.*alphalong:
            
    gp.LOG.debug(' last bin prior:')
    if (gp.lbprior) :  
        totmlastb = np.sum(denscheck[0:gp.nipol-2]) + np.sum(gp.blow[0:gp.nipol-2])
        lastb = Mpars[gp.nipol-1] + gp.blow[gp.nipol-1]
        if lastb / totmlastb > gp.lbtol :
            print('lbprior')
            return True

    return False



# check that beta is bound and not jumping
def check_delta():
    # now checking delta <= 1
    d1 = phys.delta(gp.parst.delta1)
    if max(d1)>1.:
        # correction on the fly
        # for jj in range(gp.nipol):
        #     if d1[jj] >1.:
        #         gp.parst.delta1[jj] /= 2.
        #         print('d1 too high, entry ',jj,' set to half its value')
        return True
    # prior: if last delta >= 0, the stars are not radially biased
    # if d1[-1]<0:
    #     return True
    if gp.pops == 2:
        d2 = phys.delta(gp.parst.delta2)
        if max(d2)>1.:
            # correction on the fly
            # for jj in range(gp.nipol):
            #     if d2[jj] >1.:
            #         gp.parst.delta2[jj] /= 2.
            #         print('d2 too high, entry ',jj,' set to half its value')
            return True
        # if d2[-1]<0:
        #     return True

    # now checking physical kappa: g(rvar, rfix, beta, dbetadr) >= 0
    r0 = gp.xipol
    dR = r0[1:]-r0[:-1]
    r0extl = np.array([r0[0]/6.,r0[0]/5.,r0[0]/4.,r0[0]/3.,r0[0]/2.,r0[0]/1.5])
    
    # extrapolation to the right (attention, could overshoot)
    dr0 = (r0[-1]-r0[-2])/8.
    r0extr = np.hstack([r0[-1]+dr0,r0[-1]+2*dr0,r0[-1]+3*dr0,r0[-1]+4*dr0])

    r0nu = np.hstack([r0extl,r0,dR/2.+r0[:-1],r0extr])
    r0nu.sort()

    tck0 = splrep(r0,d1*(r0**2+np.median(r0)**2),k=3) #previous: k=2, s=0.1
    betanu = splev(r0nu,tck0)/(r0nu**2+np.median(r0)**2)

    if max(betanu)>1:
        print('beta_interpolated overshooting!')
        return True
    drspl = splev(r0nu,tck0, der=1)

    dbetanudr = (drspl-betanu*2*r0nu)/(r0nu**2+np.median(r0)**2)
    for i in range(len(r0nu)-4):
        for j in range(i+1,len(r0nu)):
            if g(r0nu[j], r0nu[i], betanu[j], dbetanudr[j]) < 0:
                # print('g() < 0')
                return True


        
    # gp.LOG.debug( 'now checking delta smoothness in both directions')
    # dcheck = abs(gp.parst.delta1[1:]-gp.parst.delta1[:-1])
    # for i in range(gp.nipol-1):
    #     if dcheck[i]>gp.deltol:
    #         if not gp.d1wild:
    #             # print('delta1 too wild!') # TODO: enable print
    #             gp.d1wild = True
    #         # gf.get_working_pars()
    #         # correct: smooth out, by assigning mean value of left/right points
    #         # gp.pars.delta1[i]/=2.
    #         return True
    # if gp.pops==2:
    #     dcheck = abs(gp.parst.delta2[1:]-gp.parst.delta2[:-1])
    #     for i in range(gp.nipol-1):
    #         if dcheck[i]>gp.deltol:
    #             if not gp.d2wild:
    #                 # print('delta2 too wild!') # TODO: enable print
    #                 gp.d2wild = True
    #             # gf.get_working_pars()
    #             # gp.pars.delta2[i]/=2.
    #             return True

    return False









## check that sigma is positive
def check_sigma():
    if (not gp.mirror) :
        # now checking: Ensure positivity --> monotonicity constraint:
        if gp.parst.has_negative():
            print('not mirror: parst has negative')
            gf.get_working_pars()
            return True

    # Apply "cprior" and "gprior" on M function:
    # if (cpriorconv>0) : 
    #   if (abs(denscheck[0])>cpriorconv) : return True # TODO: determine mass
    # else:  
    #   denscheck[0] = Mminarr[0] # TODO: determine mass
    #   parst[2*nipol] = Mminarr[0]


    # TODO disc
    gp.LOG.debug('check positivity of tilt')
    if gp.geom == 'disc' and gp.mirror == False:
        for jj in range(gp.nipol):
            if gp.parst.delta[jj] < 0.:
                return True

    # Reject models with NaN sig_z:
    import math
    # Another problem with linear sampling (best avoided). 
    # If I do the "right" thing here and reject NaN models, 
    # the MCMC gets stuck. So we set sig_z = 0 where it
    # is NaN and assume that this will be penalised by the 
    # data. This assumption appears to be very good, 
    # but still it's not ideal ):
    
    small = 0.0 # min(gp.sig1_x[(gp.sig1_x > 0)])
    for jj in range(gp.nipol):
        # check for NaN
        if math.isnan(gp.sig1_x[jj]):
            gp.sig1_x[jj] = small
    # end TODO


    # check new disc code priors
    # TODO: check gp.kzmin < kz < gp.kzmax
    # TODO: fix repeated occurencies, scopes,...
    # for jj in range(gp.nipol):
    #     if gp.monotonic:
    #         nuparsu = np.zeros(gp.nipol)
    #         nuparsu[0] = nupars[0]
    #         for kk in range(1,len(nuparsu)):
    #             nuparsu[kk] = nuparsu[kk-1] + nupars[kk]
    #         nuparsu = nuparsu[::-1]
    #     else:
    #         nuparsu = nupars
    #     if nuparsu[jj] < numin: return True
    #     if nuparsu[jj] > numax: return True
   
    gp.LOG.debug( 'S-prior: ensure sigma_z(z) rises (in disc case only):')
    if gp.sprior:
        for jj in range(1,gp.nipol):
            if gp.sig1_x[jj]/gp.sig1_x[jj-1]<0.5:
                if not gp.sig2wild:
                    print('falling sigma_z')
                    gp.sig2wild = True
                return True
    return False
