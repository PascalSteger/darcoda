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
    # set all nu to known 3D data (plus some offset, possibly)
    # nupars1 = npr.normal(gp.ipol.nudat1,gp.ipol.nuerr1/2.,gp.nipol)
    nupars1 = gp.ipol.nudat1

    # * (1.+ npr.uniform(-1.,1.,gp.nipol)/10.)+gp.ipol.nudat1[-1] # [munit/pc^3]
    nuparstep1    = gp.ipol.nuerr1            # nupars1/20.
    if gp.geom == 'disc': nuparstep1[0] = 0.0 # first point stays 1 :)

    if gp.nulog: 
        nuparstep1 = np.log10(nupars1+nuparstep1)-np.log10(nupars1) # set first before nupars<0 :)
        nupars1 = np.log10(nupars1)
        # nuparstep1 = nupars1*0.+nupars1[0]/20. # abs(np.log10(np1*1.05)-np.log10(np1))
        #                                 # -1 means /10 in linear space # TODO: propto error
        if gp.geom == 'disc': nuparstep1[0] = 0.0
    if gp.pops == 2:
        # nupars2  = npr.normal(gp.ipol.nudat2,gp.ipol.nuerr2/2.,gp.nipol)
        nupars2  = gp.ipol.nudat2
        nuparstep2 = gp.ipol.nuerr1     # nupars2/20. # [munit/pc^3]
        if gp.nulog:
            nuparstep2 = np.log10(nupars2+nuparstep2)-np.log10(nupars2)
            nupars2 = np.log10(nupars2)
            # nuparstep2 = nupars2*0.+nupars2[0]/20.


    ### delta
    if gp.geom == 'sphere':
        deltapars1 = np.zeros(gp.nipol)
        deltaparstep1 = deltapars1 + 0.03
        mdelta1 = []; mdelta2 = []
        if gp.model:
            print 'TODO: disable model for delta for observed dwarfs!'
            if gp.investigate == 'walker':
                mdelta1, mdelta2 = betawalker(gp.xipol)
            elif gp.investigate == 'triaxial':
                mdelta1 = betatriax(gp.xipol)
            deltapars1 = phys.invdelta(mdelta1)
            deltaparstep1 = deltapars1*0. + 0.03
            # if gp.deltaprior:   # TODO: rename
            #     deltaparstep1 = np.zeros(gp.nipol)

        if gp.pops == 2:
            deltapars2 = np.zeros(gp.nipol)
            deltaparstep2 = deltapars2 + 0.03
            if gp.model:
                deltapars2 = phys.invdelta(mdelta2)
                deltaparstep2 = deltapars2*0. + 0.03
                # if gp.deltaprior:
                #     deltaparstep2 = np.zeros(gp.nipol)

    elif gp.geom == 'disc':
        # set tilt to zero, first approximation
        if not gp.deltaprior: # and tparsRin[0] > 0:
            # deltapars1 = np.zeros(gp.nipol) + 50.
            deltapars1 = np.zeros(gp.nipol)
            deltaparstep1 = np.zeros(gp.nipol)
            if gp.pops == 2:
                deltapars2 = np.zeros(gp.nipol)
                deltaparstep2 = np.zeros(gp.nipol)



    ### density (3D density in spherical case, K_z parameters in disc case)
    denspars = np.zeros(gp.nipol)
    if gp.poly:
        denspars[0] = gp.densstart # starting offset, set in gl_params
        # this is added with (radius)**0 to all other densities
        for i in range(1,gp.nipol):
            denspars[i] = (gp.scaledens)**i/i**gp.scalepower
        # scale high order dens stepsizes s.t. they change remarkably as well

        densparstep = denspars/100. * (np.arange(1,gp.nipol+1))**1.0
    else:
        denspars = nupars1/max(nupars1) # set to normalized density falloff
        if gp.model:
            print 'TODO: disable model for dens for observed dwarfs'
            denspars = rhowalkertot_3D(gp.xipol)   # [munit/pc^3]
            denspars = denspars * (1.+ npr.uniform(-1.,1.,gp.nipol)/15.)+denspars[-1]/2. # [munit/pc^3]
        if gp.denslog:
            denspars = np.log10(denspars)
        densparstep = denspars/20.
        

    if gp.geom == 'disc':
        Kz = -(gp.Mmodel)*2.*np.pi*gp.G1 # [(km/s)^2/kpc] = 3.24e-14m/s^2 # from data of overall surface density.
        # kzpars = np.abs(gh.deriv(Kz, gp.xipol)) # /sqrt(3) for first offset
        kzpars = phys.kappa(gp.xipol, Kz)        # which is the inverse to phys.dens()
        if max(kzpars<0.):
            print 'negative kappa, go check!'
            pdb.set_trace()

        if gp.kzsimpfix:
            from gl_disc_simple import get_kzpars
            kzpars = get_kzpars()

        denspars = kzpars[:] # no sign error, as wanted in paper, kappa>=0
        densparstep = denspars/30.
        # densparstep[0] = 0.0

        if gp.poly:
            denspars[0] = gp.densstart # starting offset, set in gl_params
            # this is added with (radius)**0 to all other densities
            for i in range(1,gp.nipol):
                denspars[i] = (gp.scaledens)**i/i**gp.scalepower

            # scale high order dens stepsizes s.t. they change remarkably as well
            densparstep = denspars/30. * (np.arange(1,gp.nipol+1))**0.75

        if gp.denslog:          # still in disc case
            denspars = np.log10(denspars)
            densparstep = denspars*0.+denspars[0]/30.
            

    ### norm
    normpars1 = 17.**2                   # TODO: meaning? why 17**2? needed for normalized nu=1 max.
    normparstep1 = normpars1/30.
    if gp.pops == 2:
        normpars2 = 10.**2
        normparstep2 = normpars2/30.

    # generate parameter class instances out of these variables
    if gp.pops==1:
        gp.pars    = Params(gp.pops, nupars1, denspars, deltapars1, normpars1)
        gp.parstep = Params(gp.pops, nuparstep1, densparstep, deltaparstep1, normparstep1)
    elif gp.pops==2:
        gp.pars    = Params(gp.pops, nupars1, denspars, deltapars1,\
                                normpars1, nupars2, deltapars2, normpars2)
        gp.parstep = Params(gp.pops, nuparstep1, densparstep, deltaparstep1,\
                                normparstep1, nuparstep2, deltaparstep2, normparstep2)


    gp.safepars.assign(gp.pars)
    gp.safeparstep.assign(gp.parstep)

    gp.parst = Params(0)
    gp.parst.assign(gp.pars)
    print 'mcmc set up'
    return
