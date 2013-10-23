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

import gl_params as gp
if gp.geom == 'sphere':
    import physics_sphere as phys
else:
    import physics_disc as phys

import gl_helper as gh
# import gl_plot as gpl
from gl_analytic import *
from gl_project import rho_INT_Rho, Rho_NORM_rho, rho_SUM_Mr


## stores all data from mock/observation files
class Datafile:

    ## constructor: set all initial values
    def __init__(self):
        self.xhalf = 1000.0        ##< half-light radius in [pc]
        ## total mass of tracers in [Msun]
        self.Mtottracer = 0.0
        ## TODO: unknown quantity
        self.DMpredict = np.array([])
        ## smallest radius of the data bins, in [pc]
        self.binmin = np.array([])
        ## biggest radius of the data bisn, in [pc]
        self.binmax = np.array([])
        ## keep center of radial bins, in [pc]
        self.Mx = np.array([])
        ## 3D summed overall mass, in [Msun]
        self.Mdat = np.array([])
        ## error of 3D summed overall mass, in [Msun]
        self.Merr = np.array([])
        ## keep center of radial bins, in [pc]
        self.Mx_2D=np.array([])
        ## 2D summed mass, in [Msun]
        self.Mdat_2D=np.array([])
        ## error of 2D summed mass, in [Msun]
        self.Merr_2D=np.array([])
        ## keep center of radial bins
        self.nux1 = np.array([])
        ## keep 3D model of the tracer density falloff, in [Msun/pc^3]
        self.nudat1=np.array([])
        ## keep errors of nudat1, in [Msun/pc^3]
        self.nuerr1=np.array([])
        ## keep center of radial bins
        self.nux1_2D=np.array([]);
        ## keep radial profile of the tracer density, averaged in 2D-rings
        self.nudat1_2D=np.array([]);
        ## keep error of nudat1_2D, in [Msun/pc^2]
        self.nuerr1_2D=np.array([])
        ## keep center of radial bins
        self.densx = np.array([]);
        ## keep overall 3D tracer density, 
        self.densdat = np.array([])
        ## keep error of densdat
        self.denserr = np.array([])
        ## keep center of radial bins, in [pc]
        self.densx_2D=np.array([])
        ## keep 2D density of overall density
        self.densdat_2D = np.array([]);
        ## keep error of overall density
        self.denserr_2D = np.array([])
        ## keep center of radial bins
        self.sigx1 = np.array([]);
        ## keep line of sight velocity dispersion profile, in [km/s]
        self.sigdat1=np.array([]);
        ## keep error of sigdat1
        self.sigerr1=np.array([])
        ## keep center of radial bins
        self.kapx1 = np.array([]);
        ## keep fourth velocity moment of the LOS velocities, pop1
        self.kapdat1=np.array([]);
        ## keep errors of kapdat1
        self.kaperr1=np.array([])
        ## keep number of tracers in pop 1
        self.tracers1 = np.array([])

        if(gp.pops==2):
            ## keep center of radial bins
            self.nux2 = np.array([]);
            ## keep 3D model of the tracer density falloff
            self.nudat2=np.array([]);
            ## keep errors of nudat2
            self.nuerr2=np.array([])
            ## keep center of radial bins
            self.nux2_2D=np.array([]);
            ## keep radial profile of the tracer density, averaged in 2D-rings
            self.nudat2_2D=np.array([]);
            ## keep error of nudat2
            self.nuerr2_2D=np.array([])
            ## keep center of radial bins
            self.sigx2 = np.array([]);
            ## keep second velocity moment, the velocity dispersion, as a function of radius
            self.sigdat2=np.array([]);
            ## keep error of sigdata2
            self.sigerr2=np.array([])
            ## keep center of radial bins
            self.kapx2 = np.array([]);
            ## keep fourth LOS velocity moment as function of radius
            self.kapdat2=np.array([]);
            ## keep kappa error as function of radius
            self.kaperr2=np.array([])
            ## keep number of tracers in pop 2
            self.tracers2 = np.array([])



        


    ## read enclosed baryonic 2D mass, convert it to surface density, convert that to 3D mass density
    def read_mass(self):
        self.Mx_2D, self.binmin, self.binmax, self.Mdat_2D,self.Merr_2D = gh.readcol5(gp.files.massfile)
        # [rcore], [totmass], [totmass]

        # switch to Munit (msun) and pc here
        self.Mx_2D      = self.Mx_2D[:]   * gp.rcore_2D[0]     # [pc]
        self.binmin     = self.binmin[:]  * gp.rcore_2D[0]
        self.binmax     = self.binmax[:]  * gp.rcore_2D[0]
        self.Mdat_2D    = self.Mdat_2D[:] * gp.totmass[0]      # [munit/pc**2]
        self.Merr_2D    = self.Merr_2D[:] * gp.totmass[0]      # [munit/pc**2]
        
        # calculate surface density 
        self.densx_2D   = self.Mx_2D       # [pc]
        # use binmax here, as mass is enclosed withing max radius of bin!
        self.densdat_2D = phys.calculate_surfdens(self.binmax, self.Mdat_2D) # [munit/pc**2]
        self.denserr_2D = self.densdat_2D * self.Merr_2D/self.Mdat_2D # [munit/pc**2]

        # deproject, already normalized to same total mass
        # possibility: read density file, is not possible as this is only 2D density
        self.densx, self.densdat, self.denserr = Rho_NORM_rho(self.densx_2D[:],self.densdat_2D[:],self.denserr_2D[:])
        # takes [pc], 2* [munit/pc**2], gives [pc], 2*[munit/pc**3]
        
        # needed for baryonic prior
        self.Mx   = self.densx[:]                       # [pc,3D]
        self.Mdat = rho_SUM_Mr(self.densx, self.densdat) # [Munit], 3D
        self.Merr = self.Merr_2D * self.Mdat/self.Mdat_2D # [Munit], 3D
        return






    ## read surface density of tracer stars, deproject, renormalize
    def read_nu(self):
        self.nux1_2D, dummy, dummy, self.nudat1_2D, self.nuerr1_2D = \
          gh.readcol5(gp.files.nufiles[1]) # component 1 in Walker
        # [rcore], [dens0], [dens0]
        self.nuerr1_2D *= gp.nuerrcorr #[dens0]

        # switch to Munit (msun) and pc here
        self.nux1_2D    = self.nux1_2D[:]    * gp.rcore_2D[1] # component 1
        self.nudat1_2D  = self.nudat1_2D[:]  * gp.dens0pc_2D[1]
        self.nuerr1_2D  = self.nuerr1_2D[:]  * gp.dens0pc_2D[1]
        
        # deproject, # takes [pc], 2* [munit/pc^2], gives [pc], 2* [munit/pc^3],
        # already normalized to same total mass
        if gp.geom=='sphere':
            self.nux1, self.nudat1, self.nuerr1 = Rho_NORM_rho(self.nux1_2D, self.nudat1_2D, self.nuerr1_2D)
        else:
            self.nux1, self.nudat1, self.nuerr1 = self.nux1_2D, self.nudat1_2D, self.nuerr1_2D


        if gp.pops == 2:
            self.nux2_2D,dummy,dummy, self.nudat2_2D, self.nuerr2_2D = gh.readcol5(gp.files.nufiles[2])
            # [rcore], [dens0], [dens0]
            self.nuerr2_2D *= gp.nuerrcorr #[dens0]
            
            # switch to Munit (msun) and pc here
            self.nux2_2D    = self.nux2_2D[:]    * gp.rcore_2D[2]
            self.nudat2_2D  = self.nudat2_2D[:]  * gp.dens0pc_2D[2]
            self.nuerr2_2D  = self.nuerr2_2D[:]  * gp.dens0pc_2D[2]
        
            # deproject, # takes [pc], 2* [munit/pc^2], gives [pc], 2* [munit/pc^3],
            # already normalized to same total mass
            # TODO: change nux2 from nux2_2D
            if gp.geom == 'sphere':
                self.nux2,self.nudat2,self.nuerr2=Rho_NORM_rho(self.nux2_2D,self.nudat2_2D,self.nuerr2_2D)
            else:
                self.nux2   = self.nux2_2D
                self.nudat2 = self.nudat2_2D
                self.nuerr2 = self.nuerr2_2D

        return



    ## read in line of sight velocity dispersion        
    def read_sigma(self):
        print(gp.files.sigfiles[1])
        self.sigx1,Dummy,Dummy,self.sigdat1,self.sigerr1 = gh.readcol5(gp.files.sigfiles[1])
        # 3*[Rcore], [maxvlos], [maxvlos]

        # change to pc, km/s here
        self.sigx1   = self.sigx1[:]   * gp.rcore_2D[1]     # [pc]
        self.sigdat1 = self.sigdat1[:] * gp.maxvlos[1]   # [km/s]
        self.sigerr1 = self.sigerr1[:] * gp.maxvlos[1]   # [km/s]
        
        self.sigerr1 /= gp.sigerrcorr #[km/s]
            
        if gp.pops==2:
            self.sigx2,dummy,dummy,self.sigdat2,self.sigerr2 = gh.readcol5(gp.files.sigfiles[2])
            # 3*[Rcore], [maxvlos], [maxvlos]

            # change to pc, km/s here
            self.sigx2   = self.sigx2[:]   * gp.rcore_2D[2]     # [pc]
            self.sigdat2 = self.sigdat2[:] * gp.maxvlos[2]   # [km/s]
            self.sigerr2 = self.sigerr2[:] * gp.maxvlos[2]   # [km/s]

            
            self.sigerr2 /= gp.sigerrcorr #[maxvlos]

    ## read in line of sight velocity kurtosis
    def read_kappa(self):
        self.kapx1,Dummy,Dummy,self.kapdat1,self.kaperr1 = gh.readcol5(gp.files.kappafiles[1])
        # 3*[Rcore], [maxvlos], [maxvlos]

        # change to pc, km/s here
        self.kapx1   = self.kapx1[:]   * gp.rcore_2D[1]     # [pc]
        self.kaperr1 /= gp.kaperrcorr # [1]
            
        if gp.pops==2:
            self.kapx2,dummy,dummy,self.kapdat2,self.kaperr2 = gh.readcol5(gp.files.kappafiles[2])
            # 3*[Rcore], [maxvlos], [maxvlos]

            # change to pc, km/s here
            self.kapx2   = self.kapx2[:]   * gp.rcore_2D[2]     # [pc]
            self.kaperr2 /= gp.kaperrcorr # [1]
        return

    ## set half-light radius
    def set_rhalf(self, xhalf):
        self.xhalf = xhalf
        return

    ## set total tracer mass
    def set_Mtottracer(self, Mtottracer):
        self.Mtottracer  = Mtottracer
        return

    ## interpolate 3D quantities to new regularly-spaced radial bins
    def interpol(self, dat):
        if gp.consttr:
            print('really? have const tracer number per bin, so keep same radii!')
            pdb.set_trace()

        gp.xmin  = min(dat.nux1); gp.xmax = max(dat.nux1) # [pc]
        # TODO: use radii from constant tracers per bin
        if gp.lograd:
            gp.xipol = np.logspace(np.log10(gp.xmin),np.log10(gp.xmax),gp.nipol) # [pc]
        else:
            gp.xipol = np.linspace(gp.xmin,gp.xmax,gp.nipol)    # [pc]

        if gp.geom == 'sphere':
            self.nux1_2D     = gp.xipol                                # [pc]
            self.nudat1_2D   = gh.ipollog(dat.nux1_2D, dat.nudat1_2D, self.nux1_2D) # [munit/pc^3]
            self.nuerr1_2D   = gh.ipollog(dat.nux1_2D, dat.nuerr1_2D, self.nux1_2D) # [munit/pc^3]


        self.nux1     = gp.xipol                                # [pc]
        self.nudat1   = gh.ipollog(dat.nux1, dat.nudat1, self.nux1) # [munit/pc^3]
        self.nuerr1   = gh.ipollog(dat.nux1, dat.nuerr1, self.nux1) # [munit/pc^3]

        self.binmin   = dat.binmin
        self.binmax   = dat.binmax
        self.Mx       = gp.xipol                                # [pc]
        self.Mx_2D    = gp.xipol
        # do not use ipollog here to avoid slightly wiggling output
        self.Mdat     = gh.ipol(dat.Mx, dat.Mdat, self.Mx)      # [munit]
        self.Merr     = gh.ipol(dat.Mx, dat.Merr, self.Mx)  # [munit]

        self.densx     = gp.xipol #[pc]
        self.densdat   = gh.ipol(dat.densx, dat.densdat, self.densx) # [munit/pc^3]
        # use density calculation only if real density wanted. above, we take Kz parameter
        # self.densdat   = phys.calculate_dens(self.Mx,self.Mdat) #[totmass/rcore**3]
        
        self.denserr   = self.Merr/self.Mdat*self.densdat # [munit/pc^3]

        self.sigx1    = gp.xipol                                # [pc]
        self.sigdat1  = gh.ipollog(dat.sigx1, dat.sigdat1, self.sigx1) # [km/s]
        self.sigerr1  = gh.ipollog(dat.sigx1, dat.sigerr1, self.sigx1) # [km/s]

        self.kapx1    = gp.xipol                                # [pc]
        self.kapdat1  = gh.ipollog(dat.kapx1, dat.kapdat1, self.kapx1) # [km/s]
        self.kaperr1  = gh.ipollog(dat.kapx1, dat.kaperr1, self.kapx1) # [km/s]

        if gp.pops==2:
            if gp.geom == 'sphere':
                self.nux2_2D   = gp.xipol                                # [pc]
                self.nudat2_2D = gh.ipollog(dat.nux2_2D, dat.nudat2_2D, self.nux2_2D) # [munit/pc^3]
                self.nuerr2_2D = gh.ipollog(dat.nux2_2D, dat.nuerr2_2D, self.nux2_2D) # [munit/pc^3]

            self.nux2    = gp.xipol # [pc]
            self.nudat2  = gh.ipollog(dat.nux2, dat.nudat2, self.nux2) #[munit/pc^3]
            self.nuerr2  = gh.ipollog(dat.nux2, dat.nuerr2, self.nux2) #[munit/pc^3]

            self.sigx2   = gp.xipol # [pc]
            self.sigdat2 = gh.ipollog(dat.sigx2, dat.sigdat2, self.sigx2) # [km/s]
            self.sigerr2 = gh.ipollog(dat.sigx2, dat.sigerr2, self.sigx2) # [km/s]

            self.kapx2   = gp.xipol # [pc]
            self.kapdat2 = gh.ipollog(dat.kapx2, dat.kapdat2, self.kapx2) # [km/s]
            self.kaperr2 = gh.ipollog(dat.kapx2, dat.kaperr2, self.kapx2) # [km/s]

        return self
        # dens0 = 0.40 #Msun/pc**3
        # dens0 = dens0/self.Mtottracer*self.rhalf**3 # in internal units
        # rDM  = gp.r_DM/self.rhalf #pc to internal units
        # denspredict = dens0*(self.Mr/rDM)**(-gp.gamma_DM)*\
        #                  (1+(self.Mr/rDM)**(gp.alpha_DM))**((gp.gamma_DM-gp.beta_DM)*1.0/gp.alpha_DM)
        # vol    = 4.*np.pi/3.*(self.Mr/rDM)**3
        # Mshell = vol*denspredict
        # self.DMpredict = np.cumsum(Mshell)
        # radii for samplings of M, nu, beta in nipol bins



    ## quick output of all 3D quantities
    def output(self):
        print('Mx = ',   gh.pretty(self.Mx), ' pc')
        print('Mdat  = ',gh.pretty(self.Mdat), ' munit')
        print('Merr = ', gh.pretty(self.Merr), ' munit')

        print('nux1  = ',gh.pretty(self.nux1), ' pc')
        print('nudat1 =',gh.pretty(self.nudat1), ' munit/pc^3')
        print('nuerr1=', gh.pretty(self.nuerr1), ' munit/pc^3')
        
        print('sigx1 = ',  gh.pretty(self.sigx1), ' pc')
        print('sigdat1 = ',gh.pretty(self.sigdat1), 'km/s')
        print('sigerr1 = ',gh.pretty(self.sigerr1), 'km/s')

        print('kapx1 = ',  gh.pretty(self.kapx1), ' pc')
        print('kapdat1 = ',gh.pretty(self.kapdat1), 'km/s')
        print('kaperr1 = ',gh.pretty(self.kaperr1), 'km/s')


        if gp.pops==2:
            print('nudat2 =',gh.pretty(self.nudat2), ' munit/pc^3')
            print('nuerr2=', gh.pretty(self.nuerr2), ' munit/pc^3')
            
            print('sigdat2 = ',gh.pretty(self.sigdat2), ' km/s')
            print('sigerr2 = ',gh.pretty(self.sigerr2), ' km/s')

            print('kapdat2 = ',gh.pretty(self.kapdat2), ' km/s')
            print('kaperr2 = ',gh.pretty(self.kaperr2), ' km/s')




    ## pickle dump 3D quantities
    def save(self,fn):
        fil = open(fn,'w')
        pickle.dump(self, fil)
        fil.close()
        return

    ## reload previously stored quantities
    def load(self,fn):
        fil = open(fn,'r')
        obj = pickle.load(fil)
        fil.close()
        
        self.Mx = obj.Mx;   self.Mdat = obj.Mdat;    self.Merr = obj.Merr
        self.densx = obj.densx;   self.densdat = obj.densdat;     self.denserr = obj.denserr
        self.nux1 = obj.nux1;     self.nudat1 = obj.nudat1;       self.nuerr1 = obj.nuerr1
        self.sigx1 = obj.sigx1;   self.sigdat1 = obj.sigdat1;     self.sigerr1 = obj.sigerr1
        self.kapx1 = obj.kapx1;   self.kapdat1 = obj.kapdat1;     self.kaperr1 = obj.kaperr1

        if gp.pops == 2:
            self.nux2 = obj.nux2;     self.nudat2 = obj.nudat2;      self.nuerr2 = obj.nuerr2
            self.sigx2 = obj.sigx2;   self.sigdat2 = obj.sigdat2;    self.sigerr2 = obj.sigerr2
            self.kapx2 = obj.kapx2;   self.kapdat2 = obj.kapdat2;    self.kaperr2 = obj.kaperr2

        return self


    ## set to given values, copy constructor
    def copyfrom(self,obj):
        gp.xipol = obj.densx
        self.Mx = obj.Mx;   self.Mdat = obj.Mdat;    self.Merr = obj.Merr
        self.binmin = obj.binmin;  self.binmax = obj.binmax
        self.Mx_2D = obj.Mx_2D; self.Mdat_2D = obj.Mdat_2D; self.Merr_2D = obj.Merr_2D
        self.densx = obj.densx;   self.densdat = obj.densdat;     self.denserr = obj.denserr
        self.densx_2D=obj.densx_2D; self.densdat_2D=obj.densdat_2D; self.denserr_2D=obj.denserr_2D
        self.nux1 = obj.nux1;     self.nudat1 = obj.nudat1;       self.nuerr1 = obj.nuerr1
        self.nux1_2D=obj.nux1_2D; self.nudat1_2D=obj.nudat1_2D;   self.nuerr1_2D=obj.nuerr1_2D
        self.sigx1 = obj.sigx1;   self.sigdat1 = obj.sigdat1;     self.sigerr1 = obj.sigerr1
        self.kapx1 = obj.kapx1;   self.kapdat1 = obj.kapdat1;     self.kaperr1 = obj.kaperr1
        

        if gp.pops == 2:
            self.nux2 = obj.nux2;     self.nudat2 = obj.nudat2;      self.nuerr2 = obj.nuerr2
            self.nux2_2D=obj.nux2_2D; self.nudat2_2D=obj.nudat2_2D;  self.nuerr2_2D=obj.nuerr2_2D
            self.sigx2 = obj.sigx2;   self.sigdat2 = obj.sigdat2;    self.sigerr2 = obj.sigerr2
            self.kapx2 = obj.kapx2;   self.kapdat2 = obj.kapdat2;    self.kaperr2 = obj.kaperr2

        return self
