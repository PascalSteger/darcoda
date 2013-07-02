#!/usr/bin/python
# (c) 2013 Pascal S.P. Steger
'''calculate approximative center of mass, assuming constant stellar mass'''

import numpy as np
import sys
import pdb

from pylab import *
import gl_params as gp
import gr_params as gpr
from gl_helper import expDtofloat
from gl_class_files import *

def run():
    print 'input:'
    print gpr.fil
    x0,y0,z0,vz0,vb0,Mg0,PM0,comp0=np.genfromtxt(gpr.fil,skiprows=0,unpack=True,\
                                                 usecols=(0,1,2,11,12,13,19,20),\
                                                 dtype="d17",\
                                                 converters={0:expDtofloat,  # x0  in pc \
                                                             1:expDtofloat,  # y0  in pc \
                                                             2:expDtofloat,  # z0  in pc \
                                                             11:expDtofloat, # vz0 in km/s\
                                                             12:expDtofloat, # vb0(LOS due binary), km/s\
                                                             13:expDtofloat, # Mg0 in Angstrom\
                                                             19:expDtofloat, # PM0 [1]\
                                                             20:expDtofloat}) # comp0 1,2,3(background)
    
    
    
    # TODO: use component 6-1 instead of 12-1 for z velocity, to include observational errors
    
    # only use stars which are members of the dwarf: exclude pop3 by construction
    pm = (PM0 >= gpr.pmsplit) * (comp0<3) # exclude outliers .  TODO: do that without cheating
    pm1 = pm*(comp0==1)                   # are updated below if gp.metalpop is set to True
    pm2 = pm*(comp0==2)                   # same same
    pm3 = pm*(comp0==3)

    
    if gp.metalpop:
        # drawing of populations based on metallicity
        # get parameters from function in pymcmetal.py
        import pymcmetal as pmc
        p,mu1,sig1,mu2,sig2 = pmc.bimodal_gauss(Mg0,pm)
        pm1, pm2 = pmc.assign_pop(Mg0,pm,p,mu1,sig1,mu2,sig2)
        # output: changed pm1, pm2

    # TODO: cutting pm_i to a maximum of ntracers particles:
    # from random import shuffle
    # pm1 = shuffle(pm1)[:ntracers1]
    # pm2 = shuffle(pm2)[:ntracers2]
    
    x  = x0[pm]; y=y0[pm]; z=z0[pm]; vz = vz0[pm]
    
    # center of mass
    com_x = np.sum(x)/(1.*len(x)) # [pc]
    com_y = np.sum(y)/(1.*len(y)) # [pc]
    com_z = np.sum(z)/(1.*len(z)) # [pc]
    com_vz = np.sum(vz)/(1.*len(vz)) # [km/s]
    print 'COM [pc]: ', com_x, com_y, com_z
    print 'VOM [km/s]', com_vz
    
    x0 -= com_x; y0 -= com_y; z0 -= com_z # [pc]
    
    
    # test scalelength of stars: as given in name, 100/10 pc, or as given in parameter file, 1.pc?
    # x0 = x0[pm1]; y0=y0[pm1]; z0=z0[pm1]
    # r0 = np.sqrt(x0**2+y0**2+z0**2)
    # plt.ion(); plt.subplot(111)
    # from gl_analytic import rhohern
    # rho0 = 100;
    # ma = 2000. #[pc]
    # rplot= np.arange(100)*ma/100.
    # rh = rhohern(rplot, 1000.,rho0, 2., 5., 1.)
    # plt.hist(r0[(r0<ma)],bins=20)
    
    
    
    
##########################################################################
    # calculate v_{LOS} after subtracting bulk motion of dwarf
    # overall line-of-sight velocity of the whole dwarf galaxy, in [km/s]
    vz0 -= com_vz #[km/s]
    
    r0 = np.sqrt(x0**2+y0**2) # [pc]
    rc = r0[pm] # [pc]
    rc.sort() # [pc]
    for i in range(len(rc)-1):
        if rc[i]>rc[i+1]: #[pc]
            print 'sorting error!'
            exit(1)
    rhalf = rc[floor(len(rc)/2)] # [pc]
    rcore = rhalf # or gpr.r_DM # [pc]
    print 'rcore = ',rcore,' pc'
    print 'max(r) = ',max(rc),' pc'
    print 'last element of r : ',rc[-1],' pc'
    print 'total number of stars: ',len(rc)
    
    x0 = x0/rcore; y0 = y0/rcore # [r_core]
    
    i = -1
    for pmn in [pm,pm1,pm2,pm3]:
        pmr = (r0<(gpr.rprior*rcore)) # TODO: read from gl_class_file
        pmn = pmn*pmr # [1]
        print "fraction of members = ",1.0*sum(pmn)/len(pmn)
        i = i+1
        x=x0[pmn]; y=y0[pmn]; z=z0[pmn]; vz=vz0[pmn]; vb=vb0[pmn]; #[1], [km/s]
        Mg=Mg0[pmn]; comp=comp0[pmn]; PMN=PM0[pmn] # [ang], [1], [1]
        m = np.ones(len(pmn))
        
        r = np.sqrt(x*x+y*y) #[r_core]
        
        # print "x y z" on first line, to interprete data later on
        crcore = open(gpr.get_params_file(i),'w')
        print >> crcore, '# rcore in [pc], surfdens_central (=dens0) in [munit/rcore**2], and in [munit/pc**2], and totmass [munit], and max(v_LOS) in [km/s]'
        print >> crcore, rcore
        crcore.close()

        print 'output: ',gpr.get_com_file(i)
        c = open(gpr.get_com_file(i),'w')
        print >> c,'# x [rcore],','y [rcore],','vLOS [km/s],','rcore = ',rcore,' pc'
        for k in range(len(x)):
            print >> c,x[k],y[k],vz[k] #[rcore], [rcore], [km/s]
        c.close()
        
        
        if not gp.testplot_read: continue
        
        ion(); subplot(111)
        res = (abs(x)<3)*(abs(y)<3)
        x = x[res]; y = y[res] #[rcore]
        en = len(x)
        if en == 0: continue
        scatter(x[:en], y[:en], c=pmn[:en], s=35, vmin=0.95, vmax=1.0, lw=0.0, alpha=0.5)
        # xscale('log'); yscale('log')
        if i == 0: colorbar()
        circ_HL=Circle((0,0), radius=rcore/rcore, fc='None', ec='g', lw=3)
        circ_DM=Circle((0,0), radius=gpr.r_DM/rcore, fc='None', ec='r', lw=3)
        gca().add_patch(circ_HL); gca().add_patch(circ_DM)
        
        # visible region
        maxr = max(np.abs(x));  mayr = max(np.abs(y)) #[rcore]
        width2 = max([maxr,mayr]) #[rcore]
        xlim([-width2,width2]); ylim([-width2,width2])
        axes().set_aspect('equal')
    
        xlabel(r'$x [R_s]$'); ylabel(r'$y [R_s]$')
        # title(gpr.fil)
        savefig(gpr.get_com_png(i))
        if gpr.showplots:
            ioff();show();clf()
    
