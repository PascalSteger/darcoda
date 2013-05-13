#!/usr/bin/python
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
'''read in data and store it in appropriate class'''

import numpy as np
import numpy.random as npr
import random
import pickle
import gl_params as gp
import gl_helper as gh
import gl_plot as gpl
import pdb
from gl_analytic import *
if gp.geom == 'sphere':
    import physics_sphere as phys
elif gp.geom == 'disc':
    import physics_disc as phys







class Datafile:
    '''store all data from mock/observation files'''



    def __init__(self):

        self.xhalf = 1000.0; self.Mtottracer = 0.0
        self.DMpredict = np.array([])
        self.Mx = np.array([]);self.Mdat = np.array([]); self.Merr = np.array([])
        self.Mx_2D=np.array([]);self.Mdat_2D=np.array([]);self.Merr_2D=np.array([])
        
        self.nux1 = np.array([]); self.nudat1=np.array([]); self.nuerr1=np.array([])
        self.nux1_2D=np.array([]);self.nudat1_2D=np.array([]);self.nuerr1_2D=np.array([])
        self.densx = np.array([]); self.densdat = np.array([]); self.denserr = np.array([])
        self.densx_2D=np.array([]); self.densdat_2D = np.array([]); self.denserr_2D = np.array([])
        self.sigx1 = np.array([]); self.sigdat1=np.array([]); self.sigerr1=np.array([])
        self.tracers1 = np.array([])

        if(gp.pops==2):
            self.nux2 = np.array([]); self.nudat2=np.array([]); self.nuerr2=np.array([])
            self.nux2_2D=np.array([]);self.nudat2_2D=np.array([]);self.nuerr2_2D=np.array([])
            self.sigx2 = np.array([]); self.sigdat2=np.array([]); self.sigerr2=np.array([])
            self.tracers2 = np.array([])



        



    def read_mass(self):
        'read encl bary 2D mass, convert it to surf dens, convert that to 3D mass density'
        self.Mx_2D,self.Mdat_2D,self.Merr_2D = gh.readcol(gp.files.massfile)
        # [rcore], [totmass], [totmass]

        # switch to Munit (msun) and pc here
        self.Mx_2D      = self.Mx_2D[:]   * gp.rcore_2D[0]     # [pc]
        self.Mdat_2D    = self.Mdat_2D[:] * gp.totmass[0]      # [munit/pc**2]
        self.Merr_2D    = self.Merr_2D[:] * gp.totmass[0]      # [munit/pc**2]
        
        # calculate surface density 
        self.densx_2D   = self.Mx_2D       # [pc]
        self.densdat_2D = phys.calculate_surfdens(self.Mx_2D, self.Mdat_2D) # [munit/pc**2]
        self.denserr_2D = self.densdat_2D * self.Merr_2D/self.Mdat_2D # [munit/pc**2]

        # deproject, already normalized to same total mass
        self.densx, self.densdat, self.denserr = phys.deproject(self.densx_2D[:],\
                                                                self.densdat_2D[:],\
                                                                self.denserr_2D[:])
        # takes [pc], 2x [munit/pc**2], gives [pc], 2*[munit/pc**3]
        


        self.Mx   = self.densx[:]                       # [pc,3D]
        self.Mdat = phys.Mr3D(self.densx, self.densdat) # [Munit], 3D
        self.Merr = self.Merr_2D * self.Mdat/self.Mdat_2D # [Munit], 3D


        # normalize is deprecated, try to run MCMC with pc, Msun, km/s instead
        # assume rcore is the same for 2D and 3D.
        # TODO: check, how much half-mass radius really is in 3D case
        # gp.rcore.append(gp.rcore_2D[0])
        # gp.dens0pc.append(self.densdat[0])

        # self.densx   = self.densx[:]/gp.rcore[0]     # [rcore], 3D
        # self.densdat = self.densdat[:]/gp.dens0pc[0]   # [dens0], 3D
        # self.denserr = self.denserr[:]/gp.dens0pc[0]   # [dens0], 3D

        # self.densx_2D   = self.densx_2D[:]/gp.rcore[0]      # [rcore], 2D
        # self.densdat_2D = self.densdat_2D[:]/gp.dens0pc_2D[0] # [dens0], 2D
        # self.denserr_2D = self.denserr_2D[:]/gp.dens0pc_2D[0] # [dens0], 2D

        


        # gpl.clf(); gpl.yscale('log');
        # TODO (low priority): check difference between both plots underneath. compare surface densities first
        # # ultimate check: plot data (green) and calculated surface density (red) from 3D data (as done in gl_plot)
        # lbound = (self.Mdat_2D - gp.dat.Merr_2D)*gp.totmass[0] #[munit]
        # ubound = (self.Mdat_2D + gp.dat.Merr_2D)*gp.totmass[0] #[munit]
        # gpl.fill_between(self.Mx_2D*gp.rcore[0], lbound, ubound, alpha=0.5, color='g')
        
        # # densdat => surfden => M2D
        # from gl_int import *
        # M2D = int_project(self.densx*gp.rcore[0], self.densdat*gp.dens0pc[0])
        # gpl.plot(self.densx*gp.rcore[0], M2D)
        







    def read_nu(self):
        'read surface density of tracer stars, deproject, renormalize'

        self.nux1_2D, self.nudat1_2D, self.nuerr1_2D = gh.readcol(gp.files.nufiles[1]) # component 1 in Walker
        # [rcore], [dens0], [dens0]
        self.nuerr1_2D *= gp.nuerrcorr #[dens0]

        # switch to Munit (msun) and pc here
        self.nux1_2D    = self.nux1_2D[:]    * gp.rcore_2D[1] # component 1
        self.nudat1_2D  = self.nudat1_2D[:]  * gp.dens0pc_2D[1]
        self.nuerr1_2D  = self.nuerr1_2D[:]  * gp.dens0pc_2D[1]
        
        # deproject, # takes [pc], 2x [munit/pc^2], gives [pc], 2x [munit/pc^3],
        # already normalized to same total mass
        if gp.geom=='sphere':
            self.nux1, self.nudat1, self.nuerr1 = phys.deproject(self.nux1_2D,\
                                                                 self.nudat1_2D,\
                                                                 self.nuerr1_2D)
        else:
            self.nux1, self.nudat1, self.nuerr1 = self.nux1_2D, self.nudat1_2D, self.nuerr1_2D
        # check mass is the same
        # totmass_2D = phys.Mr2D(self.nux1_2D*gp.rcore_2D[1], self.nudat1_2D*gp.dens0pc_2D[1])
        # totmass_3D = phys.Mr3D(self.nux1, self.nudat1)  #[totmass], 3D

        # normalize, is deprecated, work in Msun, pc, km/s from now on
        # assume rcore is the same for 2D and 3D.
        # TODO: check, how much half-mass radius really is in 3D case
        # gp.rcore.append(gp.rcore_2D[1])
        # gp.dens0pc.append(self.nudat1[0])
        # self.nux1   /= gp.rcore[1]     #[rcore], 3D
        # self.nudat1 /= gp.dens0pc[1]   #[dens0], 3D
        # self.nuerr1 /= gp.dens0pc[1]   #[dens0], 3D

        # TODO: low priority, tweak integration
        # from gl_int import *
        # Sig2D = int_surfden(self.nux1*gp.rcore[1], self.nudat1*gp.dens0pc[1])


        if gp.pops == 2:
            self.nux2_2D, self.nudat2_2D, self.nuerr2_2D = gh.readcol(gp.files.nufiles[2])
            # [rcore], [dens0], [dens0]
            self.nuerr2_2D *= gp.nuerrcorr #[dens0]
            
            # switch to Munit (msun) and pc here
            self.nux2_2D    = self.nux2_2D[:]    * gp.rcore_2D[2]
            self.nudat2_2D  = self.nudat2_2D[:]  * gp.dens0pc_2D[2]
            self.nuerr2_2D  = self.nuerr2_2D[:]  * gp.dens0pc_2D[2]
        
            # deproject, # takes [pc], 2x [munit/pc^2], gives [pc], 2x [munit/pc^3],
            # already normalized to same total mass
            # TODO: change nux2 from nux2_2D
            if gp.geom == 'sphere':
                self.nux2, self.nudat2, self.nuerr2 = phys.deproject(self.nux2_2D,\
                                                                     self.nudat2_2D,\
                                                                     self.nuerr2_2D)
            else:
                self.nux2   = self.nux2_2D
                self.nudat2 = self.nudat2_2D
                self.nuerr2 = self.nuerr2_2D

            # check mass is the same
            # totmass_2D = phys.Mr2D(self.nux2_2D*gp.rcore_2D[2], self.nudat2_2D*gp.dens0pc_2D[2])
            # totmass_3D = phys.Mr3D(self.nux2, self.nudat2)  #[totmass], 3D
            
            # normalization, deprecated
            # # assume rcore is the same for 2D and 3D.
            # # TODO: check, how much half-mass radius really is in 3D case
            # gp.rcore.append(gp.rcore_2D[2])
            # gp.dens0pc.append(self.nudat2[0])
            # self.nux2   /= gp.rcore[2]     #[rcore], 3D
            # self.nudat2 /= gp.dens0pc[2]   #[dens0], 3D
            # self.nuerr2 /= gp.dens0pc[2]   #[dens0], 3D
            
            # # to check right integration. TODO: low priority: tweak
            # # from gl_int import *
            # # Sig2D = int_surfden(self.nux2*gp.rcore[2], self.nudat2*gp.dens0pc[2])






        
    def read_sigma(self):
        'read in line of sight velocity dispersion'
        self.sigx1,self.sigdat1,self.sigerr1 = gh.readcol(gp.files.sigfiles[1])
        #[rcore], [maxvlos], [maxvlos]

        # change to pc, km/s here
        self.sigx1   = self.sigx1[:]   * gp.rcore_2D[1]     # [pc]
        self.sigdat1 = self.sigdat1[:] * gp.maxvlos[1]   # [km/s]
        self.sigerr1 = self.sigerr1[:] * gp.maxvlos[1]   # [km/s]
        
        self.sigerr1 /= gp.sigerrcorr #[km/s]
            
        if gp.pops==2:
            self.sigx2,self.sigdat2,self.sigerr2 = gh.readcol(gp.files.sigfiles[2])
            # [rcore], [maxvlos], [maxvlos]

            # change to pc, km/s here
            self.sigx2   = self.sigx2[:]   * gp.rcore_2D[2]     # [pc]
            self.sigdat2 = self.sigdat2[:] * gp.maxvlos[2]   # [km/s]
            self.sigerr2 = self.sigerr2[:] * gp.maxvlos[2]   # [km/s]

            
            self.sigerr2 /= gp.sigerrcorr #[maxvlos]
            

    def set_rhalf(self, xhalf):
        self.xhalf = xhalf

    def set_Mtottracer(self, Mtottracer):
        self.Mtottracer  = Mtottracer


    def interpol(self, dat):
        'interpolate 3D quantities to new regularly-spaced radial bins'
        gp.xmin  = min(dat.nux1); gp.xmax = max(dat.nux1) # [pc]
        if gp.lograd:
            gp.xipol = np.logspace(np.log10(gp.xmin),np.log10(gp.xmax),gp.nipol) # [pc]
        else:
            gp.xipol = np.linspace(gp.xmin,gp.xmax,gp.nipol)    # [pc]


        self.nux1_2D     = gp.xipol                                # [pc]
        self.nudat1_2D   = gh.ipollog(dat.nux1_2D, dat.nudat1_2D, self.nux1_2D) # [munit/pc^3]
        self.nuerr1_2D   = gh.ipollog(dat.nux1_2D, dat.nuerr1_2D, self.nux1_2D) # [munit/pc^3]


        self.nux1     = gp.xipol                                # [pc]
        self.nudat1   = gh.ipollog(dat.nux1, dat.nudat1, self.nux1) # [munit/pc^3]
        self.nuerr1   = gh.ipollog(dat.nux1, dat.nuerr1, self.nux1) # [munit/pc^3]

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

        if gp.pops==2:
            self.nux2_2D     = gp.xipol                                # [pc]
            self.nudat2_2D   = gh.ipollog(dat.nux2_2D, dat.nudat2_2D, self.nux2_2D) # [munit/pc^3]
            self.nuerr2_2D   = gh.ipollog(dat.nux2_2D, dat.nuerr2_2D, self.nux2_2D) # [munit/pc^3]

            self.nux2    = gp.xipol # [pc]
            self.nudat2  = gh.ipollog(dat.nux2, dat.nudat2, self.nux2) #[munit/pc^3]
            self.nuerr2  = gh.ipollog(dat.nux2, dat.nuerr2, self.nux2) #[munit/pc^3]

            self.sigx2   = gp.xipol # [pc]
            self.sigdat2 = gh.ipollog(dat.sigx2, dat.sigdat2, self.sigx2) # [km/s]
            self.sigerr2 = gh.ipollog(dat.sigx2, dat.sigerr2, self.sigx2) # [km/s]

        # dens0 = 0.40 #Msun/pc**3
        # dens0 = dens0/self.Mtottracer*self.rhalf**3 # in internal units
        # rDM  = gp.r_DM/self.rhalf #pc to internal units
        # denspredict = dens0*(self.Mr/rDM)**(-gp.gamma_DM)*\
        #                  (1+(self.Mr/rDM)**(gp.alpha_DM))**((gp.gamma_DM-gp.beta_DM)*1.0/gp.alpha_DM)
        # vol    = 4.*np.pi/3.*(self.Mr/rDM)**3
        # Mshell = vol*denspredict
        # self.DMpredict = np.cumsum(Mshell)
        # radii for samplings of M, nu, beta in nipol bins



    def output(self):
        'quick output of all 3D quantities'
        print >> 'Mx = ',   gh.pretty(self.Mx), ' pc'
        print >> 'Mdat  = ',gh.pretty(self.Mdat), ' munit'
        print >> 'Merr = ', gh.pretty(self.Merr), ' munit'

        print >> 'nux1  = ',gh.pretty(self.nux1), ' pc'
        print >> 'nudat1 =',gh.pretty(self.nudat1), ' munit/pc^3'
        print >> 'nuerr1=', gh.pretty(self.nuerr1), ' munit/pc^3'
        
        print >> 'sigx1 = ',  gh.pretty(self.sigx1), ' pc'
        print >> 'sigdat1 = ',gh.pretty(self.sigdat1), 'km/s'
        print >> 'sigerr1 = ',gh.pretty(self.sigerr1), 'km/s'
        if gp.pops==2:
            print >> 'nudat2 =',gh.pretty(self.nudat2), ' munit/pc^3'
            print >> 'nuerr2=', gh.pretty(self.nuerr2), ' munit/pc^3'
            
            print >> 'sigdat2 = ',gh.pretty(self.sigdat2), ' km/s'
            print >> 'sigerr2 = ',gh.pretty(self.sigerr2), ' km/s'




    def save(self,fn):
        'pickle dump 3D quantities'
        fil = open(fn,'w')
        pickle.dump(self, fil)
        fil.close()



    def load(self,fn):
        'reload previously stored quantities'
        fil = open(fn,'r')
        obj = pickle.load(fil)
        fil.close()
        
        self.Mx = obj.Mx;   self.Mdat = obj.Mdat;    self.Merr = obj.Merr
        self.densx = obj.densx;   self.densdat = obj.densdat;     self.denserr = obj.denserr
        self.nux1 = obj.nux1;     self.nudat1 = obj.nudat1;       self.nuerr1 = obj.nuerr1
        self.sigx1 = obj.sigx1;   self.sigdat1 = obj.sigdat1;     self.sigerr1 = obj.sigerr1

        if gp.pops == 2:
            self.nux2 = obj.nux2;     self.nudat2 = obj.nudat2;      self.nuerr2 = obj.nuerr2
            self.sigx2 = obj.sigx2;   self.sigdat2 = obj.sigdat2;    self.sigerr2 = obj.sigerr2

        return self
