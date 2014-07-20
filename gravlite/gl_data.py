#!/usr/bin/env ipython3

##
# @file
# read in data and store it in an appropriate class

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import gl_physics as phys
import gl_helper as gh
from gl_project import Rho_NORM_rho
from scipy.interpolate import splrep, splev


class Datafile:
    def __init__(self):
        ## smallest radius of the data bins, in [pc]
        self.binmin = []
        ## centers of bins
        self.rbin = []
        ## biggest radius of the data bins, in [pc]
        self.binmax = []

        ## keep mass profile
        self.Mrdat = []; self.Mrerr = []
        self.Mhalf = []; self.rhalf = []
        
        ## keep radial profile of the tracer density, averaged in 2D-rings
        self.nu = []; self.Sig = []
        ## keep error of Sigdat, in [Munit/pc^2]
        self.nuerr = []; self.Sigerr = []
        self.nuhalf = []

        ## keep line of sight velocity dispersion profile, in [km/s]
        self.sig = []
        ## keep error of sigdat
        self.sigerr = []

        ## keep fourth velocity moment of the LOS velocities
        self.kap = []
        ## keep errors of kapdat
        self.kaperr = []
        return
    ## \fn __init__(self)
    # constructor: set all initial values


    def read_Sig(self, gp):
        for pop in np.arange(gp.pops+1):
            Sigx, binmin, binmax, Sigdat, Sigerr = \
                    gh.readcol5(gp.files.Sigfiles[pop])
            # 3*[rscale], [Sig0], [Sig0]

            # switch to Munit (msun) and pc here
            Sigx    = Sigx[:]    * gp.Rscale[pop]         # [pc]
            Sigdat  = Sigdat[:]  * gp.Sig0pc[pop]          # [Munit/pc^2]
            Sigerr  = Sigerr[:]  * gp.Sig0pc[pop]          # [Munit/pc^2]

            # take the overall bins for rbin, binmin, binmax vals
            if pop == 0:
                self.rbin = Sigx                                 # [pc]
                self.binmin = binmin * gp.Rscale[pop]           # [pc]
                self.binmax = binmax * gp.Rscale[pop]           # [pc]
                gp.xipol = self.rbin                            # [pc]
                if gp.iscale >= 0:
                    gp.iscale = np.sum(self.rbin<gp.Rscale[0])
                maxr = max(self.rbin)                           #
                gp.xepol = np.hstack([self.rbin, 2*maxr, 4*maxr, 8*maxr]) #
            
            # deproject, # takes [pc], 2* [Munit/pc^2], gives [pc], 2* [Munit/pc^3],
            # already normalized to same total mass
            if gp.geom == 'sphere':
                dummy, nudat, nuerr, Mrdat = Rho_NORM_rho(self.rbin, Sigdat, Sigerr)
                self.Mrdat.append(Mrdat)     # [Munit]
                Mhalf = Mrdat[-1]/2.
                self.Mhalf.append(Mhalf) # [Munit]
                
                # spline interpolation with M as x axis:
                Mtck = splrep(np.log(Mrdat), np.log(self.binmax), s=0.01)
                r_half = np.exp(splev(np.log(Mhalf), Mtck))
                # spline interpolation of nu:
                nutck = splrep(np.log(gp.xipol), np.log(nudat), s=0.01)
                self.rhalf.append(r_half) # [pc]
                self.nuhalf.append(np.exp(splev(np.log(r_half), nutck))) # [Munit/pc^3]
            else:
                print('working in disc symmetry')
                nudat, nuerr = Sigdat, Sigerr

            self.Sig.append(Sigdat) # [Munit/pc^2]
            self.Sigerr.append(Sigerr) # [Munit/pc^2]
            self.nu.append(nudat)   # [Munit/pc^3]
            self.nuerr.append(nuerr)   # [Munit/pc^3]
        return
    ## \fn read_Sig(self, gp)
    # read surface density of tracer stars, deproject, renormalize
    # @param gp global parameters


    def read_sig(self, gp):
        for pop in np.arange(gp.pops+1):
            # print(gp.files.sigfiles[pop])
            Dummy, Dummy, Dummy, sigdat, sigerr = gh.readcol5(gp.files.sigfiles[pop])
            # 3*[Rscale], [maxvlos], [maxvlos]

            # change to km/s here
            self.sig.append(sigdat[:] * gp.maxvlos[pop]) # [km/s]
            self.sigerr.append(sigerr[:] * gp.maxvlos[pop]) # [km/s]
        return
    ## \fn read_sig(self, gp)
    # read in line of sight velocity dispersion
    # @param gp global parameters


    def read_kappa(self, gp):
        for pop in np.arange(gp.pops+1):
            Dummy, Dummy, Dummy, kapdat, kaperr = gh.readcol5(gp.files.kappafiles[pop])
            # 3*[Rscale], [1], [1]

            self.kap.append(kapdat) # [1]
            self.kaperr.append(kaperr) # [1]
        return
    ## \fn read_kappa(self, gp)
    # read in line of sight velocity kurtosis
    # @param gp global parameters


    def read_zeta(self, gp):
        for pop in np.arange(gp.pops+1):
            D,D,D,zetaa,zetaaerr=gh.readcol5(gp.files.zetaafiles[pop])
            self.zetaa.append(zetaa)
            self.zetaaerr.append(zetaaerr)

            D,D,D,zetab,zetaberr=gh.readcol5(gp.files.zetabfiles[pop])
            self.zetab.append(zetab)
            self.zetaberr.append(zetaberr)
        return
    ## \fn read_zeta(self, gp)
    # read zeta profiles
    # @param gp global parameters
    
    def __repr__(self):
        return "Datafile: radii "+self.rbin
    ## \fn __repr__(self)
    # string representation for ipython

