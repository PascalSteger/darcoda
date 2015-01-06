#!/usr/bin/env ipython3

##
# @file
# read in data and store it in an appropriate class

# (c) GPL v3 2014 Pascal Steger, psteger@phys.ethz.ch

import pdb
import numpy as np
import gi_helper as gh
import gi_project as glp
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
        self.nu = []
        self.nuerr = []
        self.nu_epol = []
        self.nuerr_epol = []
        ## to keep the nr parameters
        self.nuhalf = []
        self.nrnu = []
        self.nrnuerr = []
        self.Sig = []
        self.Sigerr = []
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

                # calculate n(r) parameters as used in gi_physics from the nu(r) profile
                rleft = gp.xfine[gp.xfine <= r_half]
                rleft = rleft[::-1]
                rright= gp.xfine[gp.xfine > r_half]
                nuleft = nudatnu[gp.xfine <= r_half]
                nuleft = nuleft[::-1]
                nuright = nudatnu[gp.xfine > r_half]

                rlast = 1.*r_half
                nulast = 1.*nuhalf
                sloperight = []
                for r0 in rright:
                    i = np.argmin(np.abs(rright - r0))
                    Deltanu = -(np.log(nuright[i]) - np.log(nulast))
                    Deltar  = np.log(rright[i]) - np.log(rlast)
                    sloperight.append(Deltanu/Deltar)
                    nulast = nuright[i]
                    rlast = rright[i]
                rlast = 1.*r_half
                nulast = 1.*nuhalf
                slopeleft = []
                # work through the array from the left, from r_half
                for r0 in rleft:
                    i = np.argmin(np.abs(rleft - r0))
                    Deltanu = np.log(nuleft[i]) - np.log(nulast)
                    Deltar  = np.log(rlast) - np.log(rleft[i])
                    slopeleft.append(Deltanu/Deltar)
                    nulast = nuleft[i]
                    rlast = rleft[i]
                # inverse order of slopeleft to have it sorted according increasing r
                slopeleft = slopeleft[::-1]
                slopes = np.hstack([slopeleft, sloperight])
                nrpar = 1.*slopes[:-1]
                #Deltalogr = np.log(gp.xfine[1:]) - np.log(gp.xfine[:-1])
                #nrpar = (slopes[1:]-slopes[:-1])/Deltalogr
                spl_nrpar = splrep(gp.xfine[:-1], nrpar, k=1)
                nre = splev( gp.xipol, spl_nrpar)
                extleft = splev(gp.xepol[0:3], spl_nrpar) # np.array([nrpar[0], nrpar[0], nrpar[0]])
                extright = splev(gp.xepol[-3:], spl_nrpar) #[nrpar[-1], nrpar[-1], nrpar[-1]])
                maxnre = max(nre)
                self.nrnu.append(np.hstack([nuhalf,
                                            nre[0],
                                            extleft,
                                            nre,
                                            extright,
                                            nre[-1]]))
                errnre = np.ones(1+len(extleft)+len(nre)+len(extright)+1)*maxnre/10.
                self.nrnuerr.append(np.hstack([nuhalf/3.,
                                               errnre]))

                # import gi_physics as phys
                # from pylab import loglog, axvline, axhline, plot, xscale, clf
                # loglog(gp.xepol, self.nu_epol[0], 'b.-', lw=1)
                # axvline(r_half)
                # axhline(nuhalf)
                # rh = phys.rho(gp.xepol, self.nrnu, 0, gp)
                # rhmin = phys.rho(gp.xepol, self.nrnu - self.nrnuerr, 0, gp)
                # rhmax = phys.rho(gp.xepol, self.nrnu + self.nrnuerr, 0, gp)
                # loglog(gp.xepol, rh, 'r.-', linewidth=2)
                # loglog(gp.xepol, rhmin, 'g.-')
                # loglog(gp.xepol, rhmax, 'g--')
                # pdb.set_trace()
                # clf()
                # plot(gp.xfine[:-1], nrpar, '.-')
                # plot(gp.xipol, nre, '.-')
                # xscale('log')
                # axvline(r_half)
                # pdb.set_trace()
                # print(nrpar)
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
        # only ever use the tracer particles for velocity information
        self.sig.append(np.zeros(gp.nipol))
        self.sigerr.append(np.zeros(gp.nipol))
        for pop in np.arange(gp.pops):
            print('read_sig on file ', gp.files.sigfiles[pop+1])
            Dummy1, Dummy2, Dummy3, sigdat, sigerr = gh.readcol5(gp.files.sigfiles[pop+1])
            # 3*[Xscale], [maxsiglos], [maxsiglos]
            # change to km/s here
            self.sig.append(sigdat[:] * gp.maxsiglos[pop+1]) # [km/s]
            self.sigerr.append(sigerr[:] * gp.maxsiglos[pop+1]) # [km/s]
        return
    ## \fn read_sig(self, gp)
    # read in line of sight velocity dispersion for tracer particles
    # @param gp global parameters

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
        return "Datafile class, call gp.dat.rbin, .Sig, .sig"
    ## \fn __repr__(self)
    # string representation for ipython
