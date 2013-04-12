#!/usr/bin/python2.7
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
'''set up initial parameters'''

import gl_params as gp
from gl_class_params import *
from gl_analytic import *
import gl_file as gf
import numpy as np
import numpy.random as npr
if gp.geom == 'sphere':
    import physics_sphere as phys
elif gp.geom == 'disc':
    import physics_disc as phys







def mcmc_init():
    # Default Initial seed guess parameters:

    ### nu
    # set all nu to known data plus some offset
    nupars1 = gp.ipol.nudat1
    # * (1.+ npr.uniform(-1.,1.,gp.nipol)/10.)+gp.ipol.nudat1[-1] # [munit/pc^3]
    nuparstep1    = nupars1/10.
    # + gp.ipol.nuerr1    # [munit/pc^3], /20 earlier on, was too high
    if gp.pops == 2:
        nupars2  = gp.ipol.nudat2
        # * (1.+ npr.uniform(-1.,1.,gp.nipol)/10.)+gp.ipol.nudat2[-1] # [munit/pc^3]
        nuparstep2 = nupars2/10.#+gp.ipol.nuerr2     # [munit/pc^3]

    #if gp.geom == 'disc':
        #numin = 1.e-3; numax = 1.
        #nupars1 = np.zeros(gp.nipol) + numin * 2.0
        #nuparstep1 = 0.1*nupars1


    ### delta
    deltapars1 = np.zeros(gp.nipol)
    deltaparstep1 = deltapars1 + 0.01
    mdelta1 = []; mdelta2 = []
    if gp.model:
        print 'TODO: disable model for delta!'
        mdelta1, mdelta2 = betawalker(gp.xipol)
        deltapars1 = phys.invdelta(mdelta1)
        deltaparstep1 = deltapars1/20.
        if gp.deltaprior:
            deltaparstep1 = np.zeros(gp.nipol)

    if gp.pops == 2:
        deltapars2 = np.zeros(gp.nipol)
        deltaparstep2 = deltapars2 + 0.01
        if gp.model:
            deltapars2 = phys.invdelta(mdelta2)
            deltaparstep2 = deltapars2/20.
            if gp.deltaprior:
                deltaparstep2 = np.zeros(gp.nipol)

    if gp.geom == 'disc':
        if not gp.deltaprior: # tparsRin[0] > 0:
            # deltapars1 = np.zeros(gp.nipol) + 50.
            deltapars1 = np.zeros(gp.nipol)



    ### density
    denspars = np.zeros(gp.nipol)
    if gp.poly:
        denspars[0] = gp.densstart # starting offset, set in gl_params
        # this is added with (radius)**0 to all other densities
        for i in range(1,gp.nipol):
            denspars[i] = (gp.scaledens)**i/i**gp.scalepower
        # scale high order dens stepsizes s.t. they change remarkably as well
        densparstep = denspars/30. * (np.arange(1,gp.nipol+1))**0.75
        
    else:
        denspars = nupars1/max(nupars1) # set to normalized density falloff
        if gp.model:
            denspars = rhowalkertot_3D(gp.xipol)   # [munit/pc^3]
            denspars = denspars * (1.+ npr.uniform(-1.,1.,gp.nipol)/15.)\
                       +denspars[-1]/2. # [munit/pc^3]
        densparstep = denspars/30.

    if gp.geom == 'disc':
        # Set up kzmin/max arrays:

        # Min/max Kz for MCMC search
        # If positive assume constant; if negative take fraction
        # of local baryonic value for that bin: 
        # kzmin = 0.0075 * (4.0 * np.pi * np.G1) * 1000.**3.
        kzminarr = np.zeros(gp.nipol) + gp.kzmin
        kzmaxarr = np.zeros(gp.nipol) + gp.kzmax
            
        # Default Initial seed guess parameters [assume flat]:
        kzpars = np.zeros(gp.nipol) + 1.2*(4*np.pi*gp.G1)*1000**3./100.
        # kzpars[0] = kzmin
        kzpars[0] = (4*np.pi*gp.G1)*1000**3./100
        # kzpars[0] = 1.2*(4*np.pi*gp.G1)*1000**3.*2.0

        if gp.kzsimpfix:
            from gl_disc_simple import get_kzpars
            kzpars = get_kzpars()

        if gp.poly:
            denspars = phys.calculate_dens(gp.xipol,kzpars)
        else:
            denspars = np.array(kzpars)
            # densparstep = np.zeros(gp.nipol) + 13.
        densparstep = denspars/100.


    ### Mslope
    Mslopepars = 0.1     if (gp.mprior<0)    else gp.mprior   # [1]
    Mslopeparstep = Mslopepars/20. if gp.mprior<0 else 0.     # [1]

    ### sigma
    # TODO: is this used anywhere in the MCMC? ask Dave, why it was introduced
    sigmaslopepars1 = 0. if gp.sigmaprior1<0 else gp.sigmaprior1
    sigmaslopeparstep1 = 0.1
    if gp.pops == 2:
        sigmaslopepars2 = 0. if gp.sigmaprior2<0 else gp.sigmaprior2
        sigmaslopeparstep2 = 0.1

    # generate parameter class instances out of these variables
    if gp.pops==1:
        gp.pars    = Params(gp.pops,nupars1,denspars,deltapars1,Mslopepars,\
                                sigmaslopepars1)
        gp.parstep = Params(gp.pops,nuparstep1,densparstep,deltaparstep1,\
                                Mslopeparstep,sigmaslopeparstep1)
    elif gp.pops==2:
        gp.pars    = Params(gp.pops,nupars1,denspars,deltapars1,Mslopepars,\
                                sigmaslopepars1,nupars2,deltapars2,sigmaslopepars2)
        gp.parstep = Params(gp.pops,nuparstep1,densparstep,deltaparstep1,\
                                Mslopeparstep,sigmaslopeparstep1,\
                                nuparstep2,deltaparstep2,sigmaslopeparstep2)

    # Logarithmic prior: 
    if gp.logprior:
        gp.lpars = gp.pars.getlog()

        # fill with same float values for two populations
        if gp.pops==1:
            # 1 component, nu, dens, delta
            gp.lparstep = Params(-1, 0.1, 0.1, 0.1, \
                                 Mslopeparstep, sigmaslopeparstep1)
        elif gp.pops==2:
            gp.lparstep = Params(-2, 0.04, 0.2, 0.02, \
                                 Mslopeparstep, sigmaslopeparstep1,\
                                 0.04, 0.04, sigmaslopeparstep2)

        if gp.geom == 'disc':
            lparstep.set_nu([ 0.6]) # 0.013
            lparstep.set_dens([1.2])
            if investigate == 'simple':  
                if gp.kzsimpfix : lparstep.set_dens(0.)
                if gp.nusimpfix : lparstep.set_nu(0.)
                
            if not gp.deltaprior:
                ltpars  = np.zeros(len(tpars)) + np.log10(200)
                ltparsw = np.zeros(len(tpars)) + np.log10(200)
    
    print 'mcmc set up'
    return
