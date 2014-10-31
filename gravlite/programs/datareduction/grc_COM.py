#!/usr/bin/env ipython3

##
# @file
# calculate approximative center of mass, assuming constant
# stellar mass 2D version

# (c) 2014 Pascal S.P. Steger

import numpy as np
import pdb
from pylab import *
ion()

import gl_helper as gh
import gl_file as gf
from gl_helper import expDtofloat, floattoint
from gl_centering import com_shrinkcircle_v

def concat_pops(M1, M2, x1, x2, y1, y2, z1, z2, vx1, vx2, vy1, vy2, vz1, vz2, gp):
    M0 = np.hstack([M1, M2])
    x0 = np.hstack([x1, x2])
    y0 = np.hstack([y1, y2])
    z0 = np.hstack([z1, z2])
    vx0 = np.hstack([vx1, vx2])
    vy0 = np.hstack([vy1, vy2])
    vz0 = np.hstack([vz1, vz2])
    N1  = min(len(x1), gp.ntracer[1-1])
    N2  = min(len(x2), gp.ntracer[2-1])
    pm1 = np.hstack([np.ones(N1, dtype=bool), np.zeros(N2, dtype=bool)])
    pm2 = np.hstack([np.zeros(N1, dtype=bool), np.ones(N2,  dtype=bool)])
    pm = pm1 + pm2
    return M0, x0, y0, z0, vx0, vy0, vz0, pm1, pm2, pm
## \fn concat_pops(M1, M2, x1, x2, y1, y2, z1, z2, vx1, vx2, vy1, vy2, vz1, vz2, gp)
# concatenate all arrays for two populations
# @param M1
# @param M2
# @param x1
# @param x2
# @param y1
# @param y2
# @param z1
# @param z2
# @param vx1
# @param vx2
# @param vy1
# @param vy2
# @param vz1
# @param vz2
# @param gp

def select_pm(M, x, y, z, vx, vy, vz, comp, pm):
    return M[pm], x[pm], y[pm], z[pm], vx[pm], vy[pm], vz[pm], comp[pm]
## \fn select_pm(M, x, y, z, vx, vy, vz, comp)
# extract only parts of the arrays given
# @param M
# @param x
# @param y
# @param z
# @param vx
# @param vy
# @param vz
# @param comp
# @param pm

def read_data(filename):
    # M0 in Msun, x0, y0, z0 in pc, vx,vy,vz in km/s
    # comp0 1,2 component (float in file, but really only int)
    M0, x0,y0,z0,vx0,vy0,vz0,comp0 = np.genfromtxt(filename, skiprows = 0, unpack = True,\
                                                   usecols=(0, 1, 2, 3, 4, 5, 6, 7),\
                                                   dtype="d17",\
                                                   converters={0:expDtofloat,\
                                                               1:expDtofloat,\
                                                               2:expDtofloat,\
                                                               3:expDtofloat,\
                                                               4:expDtofloat,\
                                                               5:expDtofloat,\
                                                               6:expDtofloat,\
                                                               7:floattoint})
    return M0, x0, y0, z0, vx0, vy0, vz0, comp0
## \fn read_data(filename)
# read and convert data
# @param filename string

def run(gp):
    import gr_params
    gpr = gr_params.Params(gp)
    print('input: ', gpr.fil)
    M0,x0,y0,z0,vx0, vy0, vz0, comp0 = read_data(gpr.fil)
    # [Msun], 3*[pc], 3*[km/s], [1]

    # assign population
    if gp.pops==2:
        pm1 = (comp0 == 1) # will be overwritten below if gp.metalpop
        pm2 = (comp0 == 2) # same same
    elif gp.pops==1:
        pm1 = (comp0 < 3)
        pm2 = (comp0 == -1) # assign none, but of same length as comp0

    # cut to subsets
    ind1 = gh.draw_random_subset(x1, gp.ntracer[1-1])
    M1, x1, y1, z1, vx1, vy1, vz1, comp1 = select_pm(M1, x1, y1, z1, vx1, vy1, vz1, comp1, ind1)

    ind2 = gh.draw_random_subset(x2, gp.ntracer[2-1])
    M2, x2, y2, z2, vx2, vy2, vz2, comp2 = select_pm(M2, x2, y2, z2, vx2, vy2, vz2, comp2, ind2)

    # use vz for no contamination, or vb for with contamination
    M0, x0, y0, z0, vx0, vy0, vz0 = concat_pops(M1, M2, x1, x2, y1, y2, z1, z2, vx1, vx2, vy1, vy2, vz1, vz2, gp)
    com_x, com_y, com_z, com_vz = com_shrinkcircle_v(x0, y0, z0, vz0, pm0) # [pc]
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
