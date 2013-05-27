#!/usr/bin/python
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
'''set up initial parameters'''

import gl_params as gp
from gl_class_params import *
from gl_analytic import *
import gl_file as gf
import gl_helper as gh
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
    nupars1 = gp.ipol.nudat1*1.05

    # * (1.+ npr.uniform(-1.,1.,gp.nipol)/10.)+gp.ipol.nudat1[-1] # [munit/pc^3]
    nuparstep1    = nupars1/10.

    # if nu is taken in log, will not want direct proportionality
    # but rather, what a 1/20. change of nu gives in log space
    if gp.nulog: 
        np1 = gp.ipol.nudat1
        nupars1 = np.log10(np1)
        nuparstep1 = abs(np.log10(np1*1.05)-np.log10(np1))
    # + gp.ipol.nuerr1    # [munit/pc^3], /20 earlier on, was too high
    if gp.pops == 2:
        nupars2  = gp.ipol.nudat2

        # * (1.+ npr.uniform(-1.,1.,gp.nipol)/10.)+gp.ipol.nudat2[-1] # [munit/pc^3]
        nuparstep2 = nupars2/10. # [munit/pc^3]
        if gp.nulog:
            np2 = gp.ipol.nudat2
            nupars2 = np.log10(np2)
            nuparstep2 = abs(np.log10(np2*1.1)-np.log10(np2))

    #if gp.geom == 'disc':
        #numin = 1.e-3; numax = 1.
        #nupars1 = np.zeros(gp.nipol) + numin * 2.0
        #nuparstep1 = 0.1*nupars1


    ### delta
    deltapars1 = np.zeros(gp.nipol)
    deltaparstep1 = deltapars1 + 0.05
    mdelta1 = []; mdelta2 = []
    if gp.model:
        print 'TODO: disable model for delta!'
        mdelta1, mdelta2 = betawalker(gp.xipol)
        deltapars1 = phys.invdelta(mdelta1)
        deltaparstep1 = deltapars1*0. + 0.05
        if gp.deltaprior:
            deltaparstep1 = np.zeros(gp.nipol)

    if gp.pops == 2:
        deltapars2 = np.zeros(gp.nipol)
        deltaparstep2 = deltapars2 + 0.05
        if gp.model:
            deltapars2 = phys.invdelta(mdelta2)
            deltaparstep2 = deltapars2*0. + 0.05
            if gp.deltaprior:
                deltaparstep2 = np.zeros(gp.nipol)

    if gp.geom == 'disc':
        # set tilt to zero, first approximation
        if not gp.deltaprior: # and tparsRin[0] > 0:
            # deltapars1 = np.zeros(gp.nipol) + 50.
            deltapars1 = np.zeros(gp.nipol)
            deltaparstep1 = np.zeros(gp.nipol)



    ### density (3D density in spherical case, K_z parameters in disc case)
    denspars = np.zeros(gp.nipol)
    if gp.poly:
        denspars[0] = gp.densstart # starting offset, set in gl_params
        # this is added with (radius)**0 to all other densities
        for i in range(1,gp.nipol):
            denspars[i] = (gp.scaledens)**i/i**gp.scalepower
        # scale high order dens stepsizes s.t. they change remarkably as well

        densparstep = denspars/20. * (np.arange(1,gp.nipol+1))**0.75
    else:
        denspars = nupars1/max(nupars1) # set to normalized density falloff
        if gp.model:
            denspars = rhowalkertot_3D(gp.xipol)   # [munit/pc^3]
            denspars = denspars * (1.+ npr.uniform(-1.,1.,gp.nipol)/15.)\
                       +denspars[-1]/2. # [munit/pc^3]
        densparstep = denspars/100.

    if gp.geom == 'disc':
        Kz = -gp.dat.Mdat*2.*np.pi*gp.G1 # [(km/s)^2/kpc] = 3.24e-14m/s^2
        # kzpars = np.abs(gh.deriv(Kz, gp.xipol)) # /sqrt(3) for first offset
        r0 = gp.xipol
        kzpars = -np.hstack([Kz[0]/r0[0], (Kz[1:]-Kz[:-1])/(r0[1:]-r0[:-1])])
        if max(kzpars<0.):
            print 'negative kappa'
            pdb.set_trace()
        # which is the inverse to phys.dens()
        if gp.denslog: kzpars = np.log10(kzpars)
        
        if gp.kzsimpfix:
            from gl_disc_simple import get_kzpars
            kzpars = get_kzpars()

        denspars = kzpars[:] # no sign error, as wanted in paper, kappa>=0
        densparstep = denspars/50.

        if gp.poly:
            denspars[0] = gp.densstart # starting offset, set in gl_params
            # this is added with (radius)**0 to all other densities
            for i in range(1,gp.nipol):
                denspars[i] = (gp.scaledens)**i/i**gp.scalepower

            # scale high order dens stepsizes s.t. they change remarkably as well
            densparstep = denspars/40. * (np.arange(1,gp.nipol+1))**0.75


    ### norm
    normpars = 17**2.                   # TODO: meaning? why 17**2?
    normparstep = normpars/30.          # tunable knob


    # generate parameter class instances out of these variables
    if gp.pops==1:
        gp.pars    = Params(gp.pops, nupars1,denspars,deltapars1,normpars)
        gp.parstep = Params(gp.pops, nuparstep1,densparstep,deltaparstep1,normparstep)
    elif gp.pops==2:
        gp.pars    = Params(gp.pops,nupars1, denspars, deltapars1,\
                                normpars, nupars2, deltapars2)
        gp.parstep = Params(gp.pops, nuparstep1, densparstep, deltaparstep1,\
                                normparstep, nuparstep2, deltaparstep2)
    gp.parst = Params(0)
    gp.parst.assign(gp.pars)
    print 'mcmc set up'
    return
