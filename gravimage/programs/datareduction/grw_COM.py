#!/usr/bin/env ipython3

##
# @file
# calculate approximative center of mass, assuming constant
# stellar mass 2D version, see grw_com for 3D

# (c) 2013 Pascal S.P. Steger

import numpy as np
import sys, pdb
from pylab import *
ion()
from random import shuffle

import gl_helper as gh
import gl_file as gf
from gl_helper import expDtofloat
from gl_centering import com_shrinkcircle_v


def concat_pops(x1, x2, y1, y2, z1, z2, vz1, vz2, gp):
    x0 = np.hstack([x1, x2])
    y0 = np.hstack([y1, y2])
    z0 = np.hstack([z1, z2])
    vz0 = np.hstack([vz1, vz2])
    N1  = min(len(x1), gp.ntracer[1-1])
    N2  = min(len(x2), gp.ntracer[2-1])
    pm1 = np.hstack([np.ones(N1, dtype=bool), np.zeros(N2, dtype=bool)])
    pm2 = np.hstack([np.zeros(N1, dtype=bool), np.ones(N2,  dtype=bool)])
    pm = pm1 + pm2
    return x0, y0, z0, vz0, pm1, pm2, pm
## \fn concat_pops(x1, x2, y1, y2, z1, z2, vz1, vz2)
# concatenate all arrays for two populations
# @param x1
# @param x2
# @param y1
# @param y2
# @param z1
# @param z2
# @param vz1
# @param vz2


def select_pm(x, y, z, comp, vz, vb, Mg, PM, pm):
    return x[pm], y[pm], z[pm], comp[pm], vz[pm], vb[pm], Mg[pm], PM[pm]
## \fn select_pm(x, y, z, comp, vz, vb, Mg, PM, pm)
# extract only parts of the arrays given
# @param x
# @param y
# @param z
# @param comp
# @param vz
# @param vb
# @param Mg
# @param PM
# @param pm


def read_data(filename):
    # x0 in pc y0 in pc z0 in pc vz0 in km/s vb0(LOS due binary), km/s
    # Mg0 in Angstrom PM0 [1] comp0 1,2,3(background) use component
    # 12-1 instead of 6-1 for z velocity, to exclude observational
    # errors

    x0,y0,z0,vb0,vz0,Mg0,PM0,comp0 = np.genfromtxt(filename, skiprows = 0, unpack = True,\
                                                   usecols=(0, 1, 2, 5, 11, 13, 19, 20),\
                                                   dtype="d17",\
                                                   converters={0:expDtofloat,\
                                                               1:expDtofloat,\
                                                               2:expDtofloat,\
                                                               5:expDtofloat,\
                                                               11:expDtofloat,\
                                                               13:expDtofloat,\
                                                               19:expDtofloat,\
                                                               20:expDtofloat})
    return x0, y0, z0, vb0, vz0, Mg0, PM0, comp0
## \fn read_data(filename)
# read and convert data
# @param filename string


def run(gp):
    import gr_params
    gpr = gr_params.Params(gp)
    print('input: ', gpr.fil)
    x0,y0,z0,vb0,vz0,Mg0,PM0,comp0 = read_data(gpr.fil)
    # [pc], [km/s], [1]

    # only use stars which are members of the dwarf: exclude pop3 by
    # construction
    pm = (PM0 >= gpr.pmsplit) # exclude foreground contamination,
                              #outliers

    x0, y0, z0, comp0, vb0, vz0, Mg0, PM0 = select_pm(x0, y0, z0, comp0, vb0, vz0, Mg0, PM0, pm)

    # assign population
    if gp.pops==2:
        pm1 = (comp0 == 1) # will be overwritten below if gp.metalpop
        pm2 = (comp0 == 2) # same same
    elif gp.pops==1:
        pm1 = (comp0 < 3)
        pm2 = (comp0 == -1) # assign none, but of same length as comp0

    if gp.metalpop:
        # drawing of populations based on metallicity get parameters
        # from function in pymcmetal.py

        import pickle
        fi = open('metalsplit.dat', 'rb')
        DATA = pickle.load(fi)
        fi.close()
        p, mu1, sig1, mu2, sig2, M, pm1, pm2 = DATA

    x1, y1, z1, comp1, vb1, vz1, Mg1, PM1 = select_pm(x0, y0, z0, comp0, vb0, vz0, Mg0, PM0, pm1)
    x2, y2, z2, comp2, vb2, vz2, Mg2, PM2 = select_pm(x0, y0, z0, comp0, vb0, vz0, Mg0, PM0, pm2)

    # cut to subsets
    ind1 = gh.draw_random_subset(x1, gp.ntracer[1-1])
    x1, y1, z1, comp1, vb1, vz1, Mg1, PM1 = select_pm(x1, y1, z1, comp1, vb1, vz1, Mg1, PM1, ind1)

    ind2 = gh.draw_random_subset(x2, gp.ntracer[2-1])
    x2, y2, z2, comp2, vb2, vz2, Mg2, PM2 = select_pm(x2, y2, z2, comp2, vb2, vz2, Mg2, PM2, ind2)

    # use vz for no contamination, or vb for with contamination
    x0, y0, z0, vz0, pm1, pm2, pm = concat_pops(x1, x2, y1, y2, z1, z2, vz1, vz2, gp)
    com_x, com_y, com_z, com_vz = com_shrinkcircle_v(x0, y0, z0, vz0, pm) # [pc]
    print('COM [pc]: ', com_x, com_y, com_z)   # [pc]
    print('VOM [km/s]', com_vz)                # [km/s]

    # from now on, work with 2D data only; z0 was only used to get
    # center in (x,y) better
    x0 -= com_x # [pc]
    y0 -= com_y # [pc]
    vz0 -= com_vz # [km/s]

    R0 = np.sqrt(x0**2+y0**2) # [pc]
    Rhalf = np.median(R0) # [pc]
    Rscale = Rhalf        # [pc] from all tracer points

    pop = -1
    for pmn in [pm, pm1, pm2]:
        pop = pop + 1                    # population number
        pmr = ( R0 < (gp.maxR*Rscale) )  # read max extension for data
                                         #(rprior*Rscale) from
                                         #gl_params
        pmn = pmn*pmr                    # [1]
        print("fraction of members = ", 1.0*sum(pmn)/len(pmn))

        x, y, z, comp, vz, vb, Mg, PMN = select_pm(x0, y0, z0, comp0, vz0, vb0, Mg0, PM0, pmn)
        R = np.sqrt(x*x+y*y)             # [pc]
        Rscalei = np.median(R)
        gf.write_Xscale(gp.files.get_scale_file(pop), Rscalei)
        gf.write_data_output(gp.files.get_com_file(pop), x/Rscalei, y/Rscalei, vz, Rscalei)

        if gpr.showplots:
            gpr.show_part_pos(x, y, pmn, Rscale)

if __name__=='__main__':
    gpr.showplots = True
    import gl_params
    gp = gl_params.Params()
    run(gp)
