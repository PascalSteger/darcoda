#!/usr/bin/env ipython3

##
# @file
# read simulation files (x,y,z, vx,vy,vz), center, cut to N particles, write 2D x,y,vlos

# (c) 2013 ETHZ, Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import random, pdb

import gr_params as gpr
import gl_units as gu
# TODO hernquist gp.G1 replacement in gu
import gl_helper as gh
from gl_centering import com_shrinkcircle_v

def run(gp):
    gu.G1__pcMsun_1km2s_2 = 1.  # as per definition
    gp.anM = 1. #
    gp.ana = 1. #

    print('grh_com: input: ', gpr.simpos)
    xall, yall, zall = np.loadtxt(gpr.simpos, skiprows=1, unpack=True) # 3*[gp.ana]
    vxall,vyall,vzall= np.loadtxt(gpr.simvel, skiprows=1, unpack=True) # 3*[gp.ana]
    nall = len(xall)                                                 # [1]

    # shuffle and restrict to ntracer random points
    ndm = int(min(gp.ntracer[0], nall-1))
    trace = random.sample(range(nall), nall)
    if gp.pops > 1:
        gh.LOG(1, 'implement more than 2 pops for hern')
        pdb.set_trace()

    PM = [1. for i in trace] # [1]=const, no prob. of membership info in dataset
    x  = [ xall[i]    for i in trace ] # [gp.ana]
    y  = [ yall[i]    for i in trace ] # [gp.ana]
    z  = [ zall[i]    for i in trace ] # [gp.ana]
    vz = [ vzall[i]   for i in trace ] # [km/s]
    PM = np.array(PM); x=np.array(x); y=np.array(y); z=np.array(z); vz=np.array(vz)

    com_x, com_y, com_z, com_vz = com_shrinkcircle_v(x,y,z,vz,PM) # 3*[gp.ana], [velocity]
    print('COM [gp.ana]: ', com_x, com_y, com_z, com_vz)

    xnew = (x-com_x) #*gp.ana      # [pc]
    ynew = (y-com_y) #*gp.ana      # [pc]
    znew = (z-com_z) # *gp.ana      # [pc]
    vznew = (vz-com_vz) #*1e3*np.sqrt(gu.G1__pcMsun_1km2s_2*gp.anM/gp.ana) # [km/s], from conversion from system with L=G=M=1

    R0 = np.sqrt(xnew**2+ynew**2)   # [pc]
    Rhalf = np.median(R0)           # [pc]
    Rscale = Rhalf                  # or gpr.r_DM # [pc]

    print('Rscale/pc = ', Rscale)

    # only for 0 (all) and 1 (first and only population)
    for pop in range(gp.pops+1):
        crscale = open(gp.files.get_scale_file(pop),'w')
        print('# Rscale in [pc],',' surfdens_central (=dens0) in [Munit/rscale**2],',\
              ' and totmass_tracers [Munit],',\
              ' and max(sigma_LOS) in [km/s]', file=crscale)
        print(Rscale, file=crscale)
        crscale.close()

        gh.LOG(2, 'grh_com: output: ', gp.files.get_com_file(pop))
        filepos = open(gp.files.get_com_file(pop), 'w')
        print('# x [Rscale]','y [Rscale]','vLOS [km/s]', file=filepos)
        for k in range(ndm):
            print(xnew[k]/Rscale, ynew[k]/Rscale, vznew[k], file=filepos)
        filepos.close()
        gh.LOG(2, '')

        # disabled as too many particles are to be plotted
        # use random subset, if at all
        #if gpr.showplots:
        #    gpr.show_part_pos(x, y, PM, Rscale)

if __name__=='__main__':
    gpr.showplots = True
    import gl_params
    gp = gl_params.Params()
    run(gp)
