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


def com_mean(x,y,pm):
    '''mean COM, weighted by probability of membership'''
    com_x = 1.*np.sum(x*pm)/np.sum(pm) # [pc]
    com_y = 1.*np.sum(y*pm)/np.sum(pm) # [pc]
    return com_x, com_y


def com_shrinkcircle(x,y,pm):
    print 'shrinking sphere'
    eps = 1e-6
    com_x = 1.*np.sum(x*pm)/np.sum(pm);    com_y = 1.*np.sum(y*pm)/np.sum(pm)
    bucom_x = 0.+com_x; bucom_y = 0.+com_y
    x -= com_x; y -= com_y
    dr = np.sqrt(com_x**2+com_y**2)
    r0 = np.sqrt(x**2+y**2)

    nit = 0; minlen = len(x)*0.666666666
    while nit < 200 and len(x) > minlen:
        nit += 1
        print 'iteration ',nit,' with ',len(x), ' particles has overall COM of: ',bucom_x,bucom_y,' with remaining offset ',dr

        # shrink sphere:
        # 1) calc radius
        r0 = np.sqrt(x**2+y**2)
        # 2) sort remaining particles
        order = np.argsort(r0)
        r0 = np.array(r0)[order]; x = np.array(x)[order]; y = np.array(y)[order]; pm = np.array(pm)[order]

        # 3) cut x,y,z,pm after 1-10%
        end = len(r0)*0.95
        r0 = r0[:end]; x = x[:end]; y = y[:end]; pm = pm[:end]

        # calculate new COM
        com_x = 1.*np.sum(x*pm)/np.sum(pm);    com_y = 1.*np.sum(y*pm)/np.sum(pm)
        dr = np.sqrt(com_x**2+com_y**2)

        # add to bucom
        bucom_x += com_x; bucom_y += com_y

        # recenter particles
        x -= com_x; y -= com_y

    return bucom_x, bucom_y


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
    pm = (PM0 >= gpr.pmsplit) # exclude foreground contamination, outliers
    PM0 = PM0[pm]
    comp0 = comp0[pm]; x0=x0[pm]; y0=y0[pm]; z0=z0[pm]; vz0=vz0[pm]; vb0=vb0[pm]; Mg0=Mg0[pm]
    pm1 = (comp0 == 1) # will be overwritten below if gp.metalpop
    pm2 = (comp0 == 2) # same same
    
    
    if gp.metalpop:
        # drawing of populations based on metallicity
        # get parameters from function in pymcmetal.py
        import pymcmetal2 as pmc
        p,mu1,sig1,mu2,sig2, M = pmc.bimodal_gauss(Mg0)
        pm1, pm2 = pmc.assign_pop(Mg0,p,mu1,sig1,mu2,sig2)
        # output: changed pm1, pm2

    # cutting pm_i to a maximum of ntracers particles:
    from random import shuffle
    ind = np.arange(len(x0))
    np.random.shuffle(ind)
    ind = ind[:gp.files.ntracer]
    x0=x0[ind]; y0=y0[ind]; z0=z0[ind]; comp0=comp0[ind]; vz0=vz0[ind]; vb0=vb0[ind]; Mg0=Mg0[ind]
    PM0 = PM0[ind]; pm1 = pm1[ind]; pm2 = pm2[ind]; pm = pm1+pm2
    
    # get center of mass with means
    #com_x, com_y = com_mean(x0,y0,PM0) # [TODO]
    
    # get COM with shrinking sphere method
    com_x, com_y = com_shrinkcircle(x0,y0,PM0)
    print 'COM [pc]: ', com_x, com_y


    com_vz = np.sum(vz0*PM0)/np.sum(PM0) # [km/s]
    print 'VOM [km/s]', com_vz


    
    x0 -= com_x; y0 -= com_y # [pc]
    vz0 -= com_vz #[km/s]
    
    r0 = np.sqrt(x0**2+y0**2) # [pc]
    rc = r0 # [pc]
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
    for pmn in [pm,pm1,pm2]:
        pmr = (r0<(gpr.rprior*rcore)) # TODO: read from gl_class_file
        pmn = pmn*pmr # [1]
        print "fraction of members = ",1.0*sum(pmn)/len(pmn)
        i = i+1
        x=x0[pmn]; y=y0[pmn]; vz=vz0[pmn]; vb=vb0[pmn]; #[1], [km/s]
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
        
        
        if not gp.showplot_readout: continue
        
        ion(); subplot(111)
        res = (abs(x)<3)*(abs(y)<3)
        x = x[res]; y = y[res] #[rcore]
        en = len(x)
        if en == 0: continue
        scatter(x[:en], y[:en], c=pmn[:en], s=35, vmin=0.95, vmax=1.0, lw=0.0, alpha=0.2)
        # xscale('log'); yscale('log')
        if i == 0: colorbar()
        circ_HL=Circle((0,0), radius=rcore/rcore, fc='None', ec='g', lw=1)
        circ_DM=Circle((0,0), radius=gpr.r_DM/rcore, fc='None', ec='r', lw=1)
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
    
if __name__=='__main__':
    # gp.showplot_readout = True
    run()
