#!/usr/bin/env ipython3

##
# @file
# calculate approximative center of mass, assuming constant stellar mass
# 3D version, see grw_COM for 2D

# (c) 2013 ETHZ Pascal S.P. Steger, psteger@phys.ethz.ch

import numpy as np
import sys, pdb


def select_pm(x, y, z, vz, pm):
    return x[pm], y[pm], z[pm], vz[pm]
## \fn select_pm(x, y, z, vz, pm)
# extract only parts of the arrays given
# @param x
# @param y
# @param z
# @param vz
# @param pm


def run(gp):
    import gr_params
    gpr = gr_params.Params(gp)

    print('input:',gpr.fil)
    x0, y0, z0, vx, vy, vz = np.transpose(np.loadtxt(gpr.fil))
    # for purely tangential beta=-0.5 models, have units of kpc instead of pc
    if gp.case == 9 or gp.case == 10:
        x0 *= 1000. # [pc]
        y0 *= 1000. # [pc]
        z0 *= 1000. # [pc]

    # cutting pm_i to a maximum of ntracers particles:
    from random import shuffle
    import gl_helper as gh
    ind1 = gh.draw_random_subset(x0, gp.ntracer[1-1])
    x0, y0, z0, vz0 = select_pm(x0, y0, z0, vz, ind1)

    PM = np.ones(len(x0)) # assign all particles the full probability of membership
    import gl_centering as glc
    com_x, com_y, com_z, com_vz = glc.com_shrinkcircle_v(x0, y0, z0, vz, PM)

    # from now on, work with 2D data only;
    # z0 was only used to get center in (x,y) better

    x0 -= com_x  # [pc]
    y0 -= com_y  # [pc]
    vz -= com_vz # [km/s]
    R0 = np.sqrt(x0*x0+y0*y0) # [pc]
    Rscale = np.median(R0) # [pc]

    import gl_file as gf
    for pop in range(gp.pops+1):      # gp.pops +1 for all components together
        pmr = (R0<(gp.maxR*Rscale))
        m = np.ones(len(R0))
        x = x0[pmr] # [pc]
        y = y0[pmr] # [pc]
        R = np.sqrt(x*x+y*y) # [pc]
        Rscalei = np.median(R)
        # print("x y z" on first line, to interprete data later on)
        gf.write_Xscale(gp.files.get_scale_file(pop), Rscalei)
        gf.write_data_output(gp.files.get_com_file(pop), x/Rscalei, y/Rscalei, vz, Rscalei)

        if gpr.showplots:
            gpr.show_part_pos(x, y, pmn, Rscale)

if __name__=='__main__':
    gpr.showplots = True
    import gl_params
    gp = gl_params.Params()
    run(gp)
