#!/usr/bin/env ipython3

##
# @file
# calculate approximative center of mass, assuming constant stellar mass
# 3D version, see grw_COM for 2D

# (c) 2013 Pascal S.P. Steger

import numpy as np
import sys, pdb

from pylab import *
ion()
import gr_params as gpr
from gl_helper import expDtofloat
from gl_class_files import *
from gl_centering import *

def run(gp):
    ## read input measurements
    print('input: ', gpr.fil)
    x0,y0,z0,vb0,vz0,Mg0,PM0,comp0=np.genfromtxt(gpr.fil,skiprows=0,unpack=True,\
                                                 usecols=(0, 1, 2, 5, 11, 13, 19, 20),\
                                                 dtype="d17",\
                                                 converters={0:expDtofloat,  # x0  in pc \
                                                             1:expDtofloat,  # y0  in pc \
                                                             2:expDtofloat,  # z0  in pc \
                                                             5:expDtofloat, # vz0 in km/s\
                                                             12:expDtofloat, # vb0(LOS due binary), km/s\
                                                             13:expDtofloat, # Mg0 in Angstrom\
                                                             19:expDtofloat, # PM0 [1]\
                                                             20:expDtofloat}) # comp0 1,2,3(background)
    # use component 12-1 instead of 6-1 for z velocity, to exclude observational errors

    # only use stars which are members of the dwarf: exclude pop3 by construction
    pm  = (PM0 >= gpr.pmsplit) # exclude foreground contamination, outliers
    PM0 = PM0[pm]
    comp0 = comp0[pm]
    x0 = x0[pm]
    y0 = y0[pm]
    z0 = z0[pm]
    vz0 = vz0[pm]; vb0 = vb0[pm]; Mg0 = Mg0[pm]

    pm1 = (comp0 == 1) # will be overwritten below if gp.metalpop
    pm2 = (comp0 == 2) # same same
    pm3 = (comp0 == 3)


    if gp.metalpop:
        # drawing of populations based on metallicity
        # get parameters from function in pymcmetal.py
        import pickle
        fi = open('metalsplit.dat', 'rb')
        DATA = pickle.load(fi)
        fi.close()
        p, mu1, sig1, mu2, sig2, M, pm1, pm2 = DATA

        # TODO: do this in python2 once before invoking run before
        #import pymcmetal as pmc
        #p, mu1, sig1, mu2, sig2, M = pmc.bimodal_gauss(Mg0)
        #pm1, pm2 = pmc.assign_pop(Mg0, p, mu1, sig1, mu2, sig2)
        #DATA = [p, mu1, sig1, mu2, sig2, M, pm1, pm2]
        #import pickle
        #fi = open('metalsplit.dat', 'wb')
        #pickle.dump(DATA, fi)
        #fi.close()

    # cutting pm_i to a maximum of ntracers particles:
    from random import shuffle
    ind = np.arange(len(x0))
    np.random.shuffle(ind)
    ind = ind[:np.sum(gp.ntracer)]

    x0 = x0[ind];   y0 = y0[ind]; z0 = z0[ind]; comp0 = comp0[ind]
    vz0 = vz0[ind]; vb0=vb0[ind]; Mg0 = Mg0[ind]
    PM0 = PM0[ind]; pm1 = pm1[ind]; pm2 = pm2[ind]; pm3 = pm3[ind];
    pm = pm1+pm2+pm3

    # get COM with shrinking sphere method
    com_x, com_y, com_z = com_shrinkcircle(x0,y0,z0,PM0)
    print('COM [pc]: ', com_x, com_y, com_z)


    com_vz = np.sum(vz0*PM0)/np.sum(PM0) # [km/s]
    print('VOM [km/s]', com_vz)

    # from now on, continue to work with 3D data. store to different files

    x0 -= com_x; y0 -= com_y; z0 -= com_z # [pc]
    vz0 -= com_vz #[km/s]

    # but still get the same radii as from 2D method, to get comparison of integration routines right
    r0 = np.sqrt(x0*x0+y0*y0+z0*z0) # [pc]
    rhalf = np.median(r0) # [pc]
    rscale = rhalf                       # or gpr.r_DM # [pc]

    print('rscale = ', rscale,  ' pc')
    print('max(R) = ', max(r0) ,' pc')
    print('last element of R : ',r0[-1],' pc')
    print('total number of stars: ',len(r0))

    pop = -1
    for pmn in [pm, pm1, pm2]:
        pmr = (r0<(gp.maxR*rscale)) # [1] based on [pc]
        pmn = pmn*pmr                  # [1]
        print("fraction of members = ", 1.0*sum(pmn)/len(pmn))
        pop = pop + 1
        x  = x0[pmn];  y = y0[pmn]; z = z0[pmn]; vz = vz0[pmn]; vb = vb0[pmn];  # [pc], [km/s]
        Mg = Mg0[pmn]; comp = comp0[pmn]; PMN = PM0[pmn]   # [ang], [1], [1]
        m = np.ones(len(pmn))

        rscalei = np.median(np.sqrt(x*x+y*y+z*z))

        # print("x y z" on first line, to interprete data later on)
        crscale = open(gp.files.get_scale_file(pop)+'_3D','w')
        print('# rscale in [pc], surfdens_central (=dens0) in [Munit/rscale0^2], and in [Munit/pc^2], and totmass_tracers [Munit], and max(sigma_LOS) in [km/s]', file=crscale)
        print(rscalei, file=crscale) # use 3 different half-light radii
        crscale.close()

        # store recentered positions and velocity
        print('output: ',gp.files.get_com_file(pop)+'_3D')
        c = open(gp.files.get_com_file(pop)+'_3D','w')
        print('# x [rscale],','y [rscale],', 'z [rscale]','vLOS [km/s],','rscale = ',rscalei,' pc', file=c)
        for k in range(len(x)):
            print(x[k]/rscalei, y[k]/rscalei, z[k]/rscalei, vz[k], file=c) # 3* [pc], [km/s]
        c.close()

        if gpr.showplots and False:
            from mpl_toolkits.mplot3d import Axes3D
            import matplotlib.pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            #res = (abs(x)<3*rscalei)*(abs(y)<3*rscalei)
            #x = x[res]; y = y[res]; z = z[res]
            en = len(x)

            ax.scatter3D(x[:en], y[:en], z[:en], c=pmn[:en], s=35, \
                    vmin=0.95, vmax=1.0, lw=0.0, alpha=0.2)

            #circ_HL=Circle((0,0), radius=rscalei, fc='None', ec='b', lw=1)
            #gca().add_patch(circ_HL)
            #circ_DM=Circle((0,0), radius=gpr.r_DM, fc='None', ec='r', lw=1)
            #gca().add_patch(circ_DM)
            pdb.set_trace()
            gpr.show_part_pos(x, y, pmn, rscalei)


if __name__=='__main__':
    gpr.showplots = True
    import gl_params
    gp = gl_params.Params()
    run(gp)
