#!/usr/bin/python2.7
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
'''check all parameters for prior constraints'''
import gl_params as gp
import gl_file as gf
import pdb
import gl_plot as gpl
if gp.geom == 'disc':
    import physics_disc as phys
else:
    import physics_sphere as phys

def check_priors():
    gp.LOG.debug(' check rising mass prior:')
    denscheck = phys.dens(gp.xipol, gp.parst.dens)
    # if max((denscheck[1:]-denscheck[:-1])/denscheck[:-1])>0.5:
    for i in range(len(denscheck)-1):
        if (denscheck[i+1]-denscheck[i])/denscheck[i] > gp.ktol and gp.geom == 'sphere':
            print 'rising dens prior, more than 200% up'

            gf.get_working_pars()
            # gp.parst.dens[i+1] *= 0.9
            # gp.parst.dens *= 1./np.sqrt(np.arange(1.,gp.nipol+1)[::-1])
            return True

    gp.LOG.debug('check that observed tracer mass is less than total mass')
    if gp.bprior:
        for jj in range(gp.nipol):
            if MM[jj] < gp.blow[jj]:
                print 'bprior'
                return True
    
    gp.LOG.debug(' last bin prior:')
    if (gp.lbprior) :  
        totmlastb = np.sum(denscheck[0:gp.nipol-2]) + np.sum(gp.blow[0:gp.nipol-2])
        lastb = Mpars[gp.nipol-1] + gp.blow[gp.nipol-1]
        if lastb / totmlastb > gp.lbtol :
            print 'lbprior'
            return True

    if (not gp.mirror) :
        gp.LOG.debug('now checking: Ensure positivity --> monotonicity constraint:')
        if gp.parst.has_negative():
            print 'not mirror: parst has negative'
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
    if gp.geom == 'dics' and gp.mirror == False:
        for jj in range(gp.nipol):
            if gp.parst.delta[jj] < 0.:
                return True

    # Reject models with NaN sig_z:
    import math
    if gp.logprior:
        if math.isnan(max(gp.sig1_x)):
            print 'nan sig_z'
            return True
        small = min(sig_z[(gp.sig1_x > 0)])
        for jj in range(gp.nipol):
            if math.isnan(gp.sig1_x[jj]):
                gp.sig1_x[jj] = small
    else:
        # Another problem with linear sampling (best avoided). 
        # If I do the "right" thing here and reject NaN models, 
        # the MCMC gets stuck. So we set sig_z = 0 where it
        # is NaN and assume that this will be penalised by the 
        # data. This assumption appears to be very good, 
        # but still it's not ideal ): 
        small = min(gp.sig1_x[(gp.sig1_x > 0)])
        for jj in range(gp.nipol):
            # check for NaN
            if math.isnan(gp.sig1_x[jj]):
                gp.sig1_x[jj] = small

    # end TODO

    gp.LOG.debug('check for extreme M slopes after rpmax')
    if (abs(gp.parst.Msl)>2. and gp.geom == 'sphere') : # TODO: include boundary for disc
        print 'Mslopepars too high'
        gp.parst.Msl *=0.9
        return True # TODO: keep in mind, change Msl during MCMC

    gp.LOG.debug( 'now checking delta <= 1')
    if max(gp.parst.delta1)>1.:
        for jj in range(gp.nipol):
            if gp.parst.delta1[jj] >1.:
                gp.parst.delta1[jj] = 0.95
                print 'delta1 too high, corrected entry ',jj,' to 0.95'
        # return True
    if gp.pops == 2 and max(gp.parst.delta2)>1.:
        for jj in range(gp.nipol):
            if gp.parst.delta2[jj] >1.:
                gp.parst.delta2[jj] = 0.95
                print 'delta2 too high, corrected entry ',jj,' to 0.95'
        # return True


    gp.LOG.debug( 'now checking delta smoothness')
    for i in range(len(gp.parst.delta1)-1):
        if abs(gp.parst.delta1[i+1]-gp.parst.delta1[i])>2./gp.nipol:
            print 'delta1 too wild!'
            gp.parst.delta1 /= 2. # TODO: isn't this overridden anyhow?
            return True
    if gp.pops==2:
        for i in range(len(gp.parst.delta2)-1):
            if abs(gp.parst.delta2[i+1]-gp.parst.delta2[i])>2./gp.nipol:
                print 'delta2 too wild!'
                gp.parst.delta2/=2.
                return True

    gp.LOG.debug( 'now checking dens > gprior')
    if gp.gprior > 0 and max(denscheck) > gp.gpriorconv:
        print 'gprior'
        return True

    gp.LOG.debug( 'Reject models with NaN sig_z: ')
    if (gp.logprior) :  
        gp.LOG.debug('now checking Nan condition')
        for jj in range(len(gp.sig1_x)):
            if (math.isnan(gp.sig1_x[jj])) :
                print 'NaN sig_z prior'
                return True
        if (gp.pops == 2) : 
            for jj in range(len(gp.sig1_x)): gp.sig1_x[jj] = 0.
            for jj in range(len(gp.sig2_x)):
                if (math.isnan(gp.sig2_r[jj])) : gp.sig2_x[jj] = 0.

    gp.LOG.debug( 'check new disc code priors' )
    # TODO: fix repeated occurencies, scopes,...
    # for jj in range(gp.nipol):
    #     # if kzparsu[jj] < kzminarr[jj]/100.: return True
    #     # if denarr[jj] < kzminarr[jj]: return True
    #     # if denarr[jj] > kzmaxarr[jj]: return True
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
   
    gp.LOG.debug(' now check regularisation priors')
    if gp.rprior:
        for jj in range(1,gp.nipol):
            if abs(gp.parst.nu1[jj] - gp.parst.nu1[jj-1])/gp.parst.nu1[jj] > gp.nutol:
                return True
    
    gp.LOG.debug( 'S-prior: ensure sigma_z(z) rises (in disc case only):')
    if gp.sprior:
        for jj in range(1,gp.nipol):
            if gp.sig1_x[jj] < gp.sig1_x[jj-1]:
                return True

    return False
