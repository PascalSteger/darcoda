#!/usr/bin/python
import numpy as np
import numpy.random as npr
import random
import pickle
import gl_params as gp
import gl_helper as gh
import pdb
if gp.geom == 'sphere':
    import physics_sphere as phys
elif gp.geom == 'disc':
    import physics_disc as phys

class Datafile:
    '''store all data from mock/observation files'''
    def __init__(self):
        self.Mx = np.array([])
        self.xhalf = 1000.0; self.Mtottracer = 0.0
        self.DMpredict = np.array([])
        self.Mdat = np.array([]); self.Merr = np.array([])
        
        self.nux1 = np.array([]); self.nudat1=np.array([]); self.nuerr1=np.array([])
        self.densx = np.array([]); self.densdat = np.array([]); self.denserr = np.array([])
        self.sigx1 = np.array([]); self.sigdat1=np.array([]); self.sigerr1=np.array([])
        self.tracers1 = np.array([])

        if(gp.pops==2):
            self.nux2 = np.array([]); self.nudat2=np.array([]); self.nuerr2=np.array([])
            self.sigx2 = np.array([]); self.sigdat2=np.array([]); self.sigerr2=np.array([])
            self.tracers2 = np.array([])
        
    def read_nu(self):
        self.nux1,self.nudat1,self.nuerr1 = gh.readcol(gp.files.nufiles[0])
        # [rcore], [dens0], [dens0]
        self.nuerr1 *= gp.nuerrcorr #[dens0]
        if gp.pops == 2:
            self.nux2,self.nudat2,self.nuerr2 = gh.readcol(gp.files.nufiles[1])
            # [rcore], [dens0], [dens0]
            self.nuerr2 *= gp.nuerrcorr

    def read_mass(self):
        self.Mx,self.Mdat,self.Merr = gh.readcol(gp.files.massfile)
        #[rcore], [totmass], [totmass]
        self.densx   = self.Mx #[rcore]
        self.densdat = phys.calculate_dens(self.Mdat,self.Mx) #[totmass/rcore**3]
        self.denserr = self.Merr/self.Mdat*self.densdat #[totmass/rcore**3]

    def read_sigma(self):
        self.sigx1,self.sigdat1,self.sigerr1 = gh.readcol(gp.files.sigfiles[0])
        #[rcore], [maxvlos], [maxvlos]
        self.sigerr1 /= gp.sigerrcorr #[maxvlos]
        if gp.pops==2:
            self.sigx2,self.sigdat2,self.sigerr2 = gh.readcol(gp.files.sigfiles[1])
            # [rcore], [maxvlos], [maxvlos]
            self.sigerr2 /= gp.sigerrcorr #[maxvlos]
            
    def set_rhalf(self, xhalf):
        self.xhalf = xhalf

    def set_Mtottracer(self, Mtottracer):
        self.Mtottracer  = Mtottracer

    def interpol(self, dat):
        gp.xmin  = min(dat.nux1); gp.xmax = max(dat.nux1) #[rcore]
        binlength = (gp.xmax-gp.xmin)/(gp.nipol-1.)
        # TODO: change s.t. binning is proportional to tracer density,
        # i.e. more bins where more tracers
        gp.xipol = np.arange(gp.nipol)*binlength + gp.xmin #[rcore]
        
        self.nux1     = gp.xipol #[rcore]
        self.nudat1   = gh.ipollog(dat.nux1, dat.nudat1, self.nux1) #[dens0]
        self.nuerr1   = gh.ipollog(dat.nux1, dat.nuerr1, self.nux1) #[dens0]

        self.Mx       = gp.xipol #[rcore]
        # do not use ipollog here to avoid slightly wiggling output
        if gp.scalemass:
            self.Mdat = gh.ipol(dat.Mx,dat.Mdat/max(dat.Mdat),self.Mx) #[maxM]
        else:
            self.Mdat = gh.ipol(dat.Mx, dat.Mdat, self.Mx) #[totmass]
        self.Merr     = gh.ipol(dat.Mx, dat.Merr, self.Mx) #[totmass]

        self.densx     = gp.xipol #[rcore]
        self.densdat   = gh.ipol(dat.densx, dat.densdat, self.densx)
        # use density calculation only if real density wanted. above, we take Kz parameter
        # self.densdat   = phys.calculate_dens(self.Mdat,self.Mx) #[totmass/rcore**3]
        
        self.denserr   = self.Merr/self.Mdat*self.densdat #[totmass/rcore**3]

        self.sigx1    = gp.xipol #[rcore]
        self.sigdat1  = gh.ipollog(dat.sigx1, dat.sigdat1, self.sigx1) #[maxvlos]
        self.sigerr1  = gh.ipollog(dat.sigx1, dat.sigerr1, self.sigx1) #[maxvlos]

        if gp.pops==2:
            self.nux2    = gp.xipol #[rcore]
            self.nudat2  = gh.ipollog(dat.nux2, dat.nudat2, self.nux2) #[dens0]
            self.nuerr2  = gh.ipollog(dat.nux2, dat.nuerr2, self.nux2) #[dens0]

            self.sigx2   = gp.xipol #[rcore]
            self.sigdat2 = gh.ipollog(dat.sigx2, dat.sigdat2, self.sigx2) #[maxvlos]
            self.sigerr2 = gh.ipollog(dat.sigx2, dat.sigerr2, self.sigx2) #[maxvlos]

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
        print 'Mx = ',   gh.pretty(self.Mx)
        print 'Mdat  = ',gh.pretty(self.Mdat)
        print 'Merr = ', gh.pretty(self.Merr)

        print 'nux1  = ',gh.pretty(self.nux1)
        print 'nudat1 =',gh.pretty(self.nudat1)
        print 'nuerr1=', gh.pretty(self.nuerr1)
        
        print 'sigx1 = ',  gh.pretty(self.sigx1)
        print 'sigdat1 = ',gh.pretty(self.sigdat1)
        print 'sigerr1 = ',gh.pretty(self.sigerr1)
        if gp.pops==2:
            print 'nudat2 =',gh.pretty(self.nudat2)
            print 'nuerr2=', gh.pretty(self.nuerr2)
            
            print 'sigdat2 = ',gh.pretty(self.sigdat2)
            print 'sigerr2 = ',gh.pretty(self.sigerr2)

    def save(self,fn):
        fil = open(fn,'w')
        pickle.dump(self, fil)
        fil.close()

    def load(self,fn):
        fil = open(fn,'r')
        obj = pickle.load(fil)
        fil.close()
        
        self.Mx = obj.Mx
        self.Mdat = obj.Mdat
        self.Merr = obj.Merr

        self.densx = obj.densx
        self.densdat = obj.densdat
        self.denserr = obj.denserr

        self.nux1 = obj.nux1
        self.nudat1 = obj.nudat1
        self.nuerr1 = obj.nuerr1

        self.sigx1 = obj.sigx1
        self.sigdat1 = obj.sigdat1
        self.sigerr1 = obj.sigerr1

        if gp.pops == 2:
            self.nux2 = obj.nux2
            self.nudat2 = obj.nudat2
            self.nuerr2 = obj.nuerr2

            self.sigx2 = obj.sigx2
            self.sigdat2 = obj.sigdat2
            self.sigerr2 = obj.sigerr2

        return self
