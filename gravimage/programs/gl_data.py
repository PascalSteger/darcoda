#!/usr/bin/env ipython3

##
# @file
# read in data and store it in an appropriate class

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import pdb
import numpy as np
import gl_helper as gh
import gl_project as glp
from scipy.interpolate import splrep, splev

def introduce_points_in_between(r0, gp):
    rmin = np.log10(min(r0))
    rmax = np.log10(max(r0))
    return np.logspace(rmin, rmax, gp.nfine)
## \fn introduce_points_in_between(r0, gp)
# get gp.fine points logarithmically spaced points
# @param r0 [pc] gp.xipol
# @param gp global parameter

class Datafile:
    def __init__(self):
        ## smallest radius of the data bins, in [pc]
        self.binmin = []
        ## centers of bins
        self.rbin = []
        ## biggest radius of the data bins, in [pc]
        self.binmax = []

        ## keep mass profile
        self.Mr = []; self.Mrerr = []
        self.Mhalf = []; self.rhalf = []

        ## keep radial profile of the tracer density, averaged in 2D-rings
        self.nu = [];    self.nu_epol = []
        self.nuerr = []; self.nuerr_epol = []

        self.nuhalf = []

        self.Sig = [];     self.Sigerr = []

        ## keep line of sight velocity dispersion profile, in [km/s]
        self.sigz2 = []
        ## keep error of sigdat
        self.sigz2err = []

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
            Sigx    = Sigx[:]    * gp.Xscale[pop]         # [pc]
            Sigdat  = Sigdat[:]  * gp.Sig0pc[pop]          # [Munit/pc^2]
            Sigerr  = Sigerr[:]  * gp.Sig0pc[pop]          # [Munit/pc^2]

            # take the overall bins for rbin, binmin, binmax vals
            if pop == 0:
                self.rbin = Sigx                                 # [pc]
                self.binmin = binmin * gp.Xscale[pop]           # [pc]
                self.binmax = binmax * gp.Xscale[pop]           # [pc]
                gp.xipol = self.rbin                            # [pc]

                minr = min(self.rbin)                           # [pc]
                maxr = max(self.rbin)                           # [pc]
                gp.xepol = np.hstack([minr/8., minr/4., minr/2.,\
                                      self.rbin, \
                                      2*maxr, 4*maxr, 8*maxr]) # [pc]
                gp.xfine = introduce_points_in_between(gp.xepol, gp)
            # deproject,
            # takes [pc], 2* [Munit/pc^2], gives [pc], 2* [Munit/pc^3],
            # already normalized to same total mass
            if gp.geom == 'sphere':
                Sigdatnu, Sigerrnu = gh.complete_nu(self.rbin, \
                                                    Sigdat, Sigerr, gp.xfine)
                dummyx, nudatnu, nuerrnu, Mrnu = glp.Sig_NORM_rho(gp.xfine, \
                                                                Sigdatnu, Sigerrnu,\
                                                                gp)
                self.nu_epol.append(gh.linipollog(gp.xfine, nudatnu, gp.xepol))
                self.nuerr_epol.append(gh.linipollog(gp.xfine, nuerrnu, gp.xepol))
                nudat = gh.linipollog(gp.xfine, nudatnu, gp.xipol)
                nuerr = gh.linipollog(gp.xfine, nuerrnu, gp.xipol)
                Mr = gh.linipollog(gp.xfine, Mrnu, gp.xipol)
                self.Mr.append(Mr) # [Munit]
                Mhalf = Mr[-1]/2.     # [Munit]
                self.Mhalf.append(Mhalf) # [Munit]

                # spline interpolation with M as x axis, to get half-mass of system:
                splpar_M = splrep(np.log(Mr), np.log(self.binmax), s=0.01)
                r_half = np.exp(splev(np.log(Mhalf), splpar_M)) # [pc]
                self.rhalf.append(r_half) # [pc]

                # spline interpolation of nu at r_half:
                splpar_nu = splrep(np.log(gp.xipol), np.log(nudat), s=0.01)
                nuhalf = np.exp(splev(np.log(r_half), splpar_nu)) # [pc]
                self.nuhalf.append(nuhalf)
                # [Munit/pc^3]
            else:
                gh.LOG(1, 'working in disc symmetry: reading nu directly')
                dummy1, dummy2, dummy3, nudat, nuerr = \
                        gh.readcol5(gp.files.nufiles[pop])
                self.nuhalf.append(nudat[round(len(nudat)/2)]) #HS ToDo: check validity of this

            self.Sig.append(Sigdat)    # [Munit/pc^2]
            self.Sigerr.append(Sigerr) # [Munit/pc^2]
            self.nu.append(nudat)      # [Munit/pc^3]
            self.nuerr.append(nuerr)   # [Munit/pc^3]
        return
    ## \fn read_Sig(self, gp)
    # read surface density of tracer stars, deproject, renormalize
    # @param gp global parameters


    def read_sig(self, gp):
        for pop in np.arange(gp.pops+1):
            print('read sig')
            pdb.set_trace()
            Dummy1, Dummy2, Dummy3, sigdat, sigerr = gh.readcol5(gp.files.sigfiles[pop])
            # 3*[Xscale], [maxsiglos], [maxsiglos]

            # change to km/s here
            self.sig.append(sigdat[:] * gp.maxsiglos[pop]) # [km/s]
            self.sigerr.append(sigerr[:] * gp.maxsiglos[pop]) # [km/s]
        return
    ## \fn read_sig(self, gp)
    # read in line of sight velocity dispersion
    # @param gp global parameters


    def read_sigz2(self, gp):
        dummy, dummy, dummy, sigz2dat, sigz2err = gh.readcol5(gp.files.sigfiles[0])
        self.sigz2.append(sigz2dat[:]) # [km/s]
        self.sigz2err.append(sigz2err[:]) # [km/s]

        return
    ## \fn read_sig(self, gp)
    # read in line of sight velocity dispersion
    # @param gp global parameters
    # H Silverwood 20/11/14

    def read_nu(self, gp):
        bincenters, binmins, binmaxs, nudat, nuerr = gh.readcol5(gp.files.nufiles[0])
        self.nu.append(nudat[:]) # [#stars/kpc^3]
        self.nuerr.append(nuerr[:])
        gp.z_bincenters = bincenters # [kpc]
        gp.z_binmins = binmins
        gp.z_binmaxs = binmaxs
        #gp.z_all_pts = np.append(np.append([0.0], bincenters), [binmaxs[-1]]) #[kpc]
        gp.z_all_pts = np.append([0.0], bincenters)  #[kpc]


    def read_kappa(self, gp):
        for pop in np.arange(gp.pops+1):
            Dummy1, Dummy2, Dummy3, kapdat, kaperr = gh.readcol5(gp.files.kappafiles[pop])
            # 3*[Xscale], [1], [1]

            self.kap.append(kapdat) # [1]
            self.kaperr.append(kaperr) # [1]
        return
    ## \fn read_kappa(self, gp)
    # read in line of sight velocity kurtosis
    # @param gp global parameters


    def read_zeta(self, gp):
        for pop in np.arange(gp.pops+1):
            D1,D2,D3,zetaa,zetaaerr=gh.readcol5(gp.files.zetaafiles[pop])
            self.zetaa.append(zetaa)
            self.zetaaerr.append(zetaaerr)

            D1,D2,D3,zetab,zetaberr=gh.readcol5(gp.files.zetabfiles[pop])
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
