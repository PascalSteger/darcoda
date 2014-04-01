#!/usr/bin/env python3

##
# @file
# read in data and store it in an appropriate class

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import pdb
import numpy as np
import numpy.random as npr
import random
import pickle

import gl_physics as phys
import gl_helper as gh
from gl_analytic import *
from gl_project import rho_INT_Rho, Rho_NORM_rho, rho_SUM_Mr

class Datafile:
    def __init__(self):
        ## smallest radius of the data bins, in [pc]
        self.binmin = []
        ## centers of bins
        self.rbin = []
        ## biggest radius of the data bins, in [pc]
        self.binmax = []

        ## keep mass profile
        self.Mdat = []; self.Merr = []
        
        ## keep radial profile of the tracer density, averaged in 2D-rings
        self.nudat = []; self.Nudat = []
        ## keep error of Nudat, in [Msun/pc^2]
        self.nuerr = []; self.Nuerr = []

        ## keep line of sight velocity dispersion profile, in [km/s]
        self.sigdat = []
        ## keep error of sigdat
        self.sigerr = []

        ## keep fourth velocity moment of the LOS velocities
        self.kapdat = []
        ## keep errors of kapdat
        self.kaperr = []
        return
    ## \fn __init__(self)
    # constructor: set all initial values

    def read_nu(self, gp):
        for pop in np.arange(gp.pops+1):
            Nux, binmin, binmax, Nudat, Nuerr = \
                    gh.readcol5(gp.files.nufiles[pop]) # component 1
            # 3*[rscale], [Nu0], [Nu0]
            Nuerr *= gp.nuerrcorr # [rho0]
            # switch to Munit (msun) and pc here
            Nux    = Nux[:]    * gp.Rscale[pop]         # [pc]
            Nudat  = Nudat[:]  * gp.Nu0pc[pop]          # [Msun/pc^2]
            Nuerr  = Nuerr[:]  * gp.Nu0pc[pop]          # [Msun/pc^2]

            # take the overall bins for rbin, binmin, binmax vals
            if pop == 0:
                self.rbin = Nux                                 # [pc]
                self.binmin = binmin * gp.Rscale[pop]           # [pc]
                self.binmax = binmax * gp.Rscale[pop]           # [pc]
                gp.xipol = self.rbin                            # [pc]
                maxr = max(self.rbin)                           #
                gp.xepol = np.hstack([self.rbin, 2*maxr, 4*maxr, 8*maxr]) #
            
            # deproject, # takes [pc], 2* [munit/pc^2], gives [pc], 2* [munit/pc^3],
            # already normalized to same total mass
            if gp.investigate == 'walk' or gp.investigate == 'obs' or gp.investigate == 'gaia':
                dummy, nudat, nuerr = Rho_NORM_rho(self.rbin, Nudat, Nuerr)
            else:
                print('working in disc symmetry')
                nudat, nuerr = Nudat, Nuerr

            self.Nudat.append(Nudat) # [Msun/pc^2]
            self.Nuerr.append(Nuerr) # [Msun/pc^2]
            self.nudat.append(nudat) # [Msun/pc^3]
            self.nuerr.append(nuerr) # [Msun/pc^3]
        return
    ## \fn read_nu(self)
    # read surface density of tracer stars, deproject, renormalize


    def read_sigma(self, gp):
        for pop in np.arange(gp.pops+1):
            # print(gp.files.sigfiles[pop])
            Dummy, Dummy, Dummy, sigdat, sigerr = gh.readcol5(gp.files.sigfiles[pop])
            # 3*[Rscale], [maxvlos], [maxvlos]
            sigerr /= gp.sigerrcorr     # [maxvlos]
            # change to km/s here
            self.sigdat.append(sigdat[:] * gp.maxvlos[pop]) # [km/s]
            self.sigerr.append(sigerr[:] * gp.maxvlos[pop]) # [km/s]
        return
    ## \fn read_sigma(self)
    # read in line of sight velocity dispersion


    def read_kappa(self, gp):
        for pop in np.arange(gp.pops+1):
            Dummy, Dummy, Dummy, kapdat, kaperr = gh.readcol5(gp.files.kappafiles[pop])
            # 3*[Rscale], [1], [1]
            kaperr /= gp.kaperrcorr # [1]
            self.kapdat.append(kapdat) # [1]
            self.kaperr.append(kaperr) # [1]
        return
    ## \fn read_kappa(self)
    # read in line of sight velocity kurtosis

    
    def set_rhalf(self, xhalf):
        self.xhalf = xhalf
        return
    ## \fn set_rhalf(self, xhalf)
    # set half-light radius, [pc]


    def set_Mtottracer(self, Mtottracer):
        self.Mtottracer  = Mtottracer
        return
    ## \fn set_Mtottracer(self, Mtottracer)
    # set total tracer mass


    def save(self, fn):
        fil = open(fn, 'wb')
        pickle.dump(self, fil)
        fil.close()
        return
    ## \fn save(self, fn)
    # pickle dump 3D quantities


    def load(self, fn):
        fil = open(fn,'r')
        obj = pickle.load(fil)
        fil.close()
        self.binmin = obj.binmin;  self.rbin   = obj.rbin;  self.binmax = obj.binmax
        self.nudat  = obj.nudat;   self.nuerr  = obj.nuerr
        self.kapdat = obj.kapdat;  self.kaperr = obj.kaperr

        return self
    ## \fn load(self, fn)
    # reload previously stored quantities
    # @param fn filename
    

    def copyfrom(self, obj):
        self.binmin = obj.binmin;  self.rbin   = obj.rbin;  self.binmax = obj.binmax
        self.nudat  = obj.nudat;   self.nuerr  = obj.nuerr
        self.kapdat = obj.kapdat;  self.kaperr = obj.kaperr

        return self
    ## \fn copyfrom(self, obj)
    # set to given values, copy constructor
    
# \class Datafile
# stores all data from mock/observation files
