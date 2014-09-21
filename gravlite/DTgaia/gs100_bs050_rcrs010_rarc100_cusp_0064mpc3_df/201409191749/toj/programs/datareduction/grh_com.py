#!/usr/bin/env ipython3

##
# @file
# read simulation files (x,y,z, vx,vy,vz), center, cut to N particles, write 2D x,y,vlos

# (c) 2013 ETHZ, Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import random, pdb
import gr_params as gpr

Mscale = 1.
ascale = 1.

def run(gp):
    print('grh_com: input: ', gpr.simpos)
    xall, yall, zall = np.loadtxt(gpr.simpos, skiprows=1, unpack=True) # 3*[ascale]
    vxall,vyall,vzall= np.loadtxt(gpr.simvel, skiprows=1, unpack=True) # 3*[ascale]
    nall = len(xall)                                                 # [1]

    # shuffle and restrict to ntracer random points
    if(gp.ntracer[0]>0):
        ndm = min(gp.ntracer[0], nall-1)
        trace = random.sample(range(nall), nall)
    else:
        ndm = min(gp.ntracer[1], nall-1)
        trace = random.sample(range(nall), nall)
    if(gp.ntracer[0]+gp.ntracer[1] == 0):
        ndm = nall
        trace = np.arange(nall)

    PM = [1. for i in trace] # [1]=const, no prob. of membership info in dataset
    x  = [ xall[i]    for i in trace ] # [ascale]
    y  = [ yall[i]    for i in trace ] # [ascale]
    z  = [ zall[i]    for i in trace ] # [ascale]
    vz = [ vzall[i]   for i in trace ] # [km/s]
    PM = np.array(PM); x=np.array(x); y=np.array(y); z=np.array(z); vz=np.array(vz)

    from gl_centering import com_shrinkcircle_v
    com_x, com_y, com_z, com_vz = com_shrinkcircle_v(x,y,z,vz,PM) # 3*[ascale], [velocity]
    print('COM [ascale]: ', com_x, com_y, com_z, com_vz)

    xnew = (x-com_x)*ascale      # [pc]
    ynew = (y-com_y)*ascale      # [pc]
    znew = (z-com_z)*ascale      # [pc]
    vznew = (vz-com_vz)*np.sqrt(gp.G1*Mscale/ascale) # [km/s], from conversion from system with L=G=M=1

    R0 = np.sqrt(xnew**2+ynew**2)   # [pc]
    Rhalf = np.median(R0) # [pc]
    Rscale = Rhalf                  # or gpr.r_DM # [pc]

    print('Rscale = ', Rscale)

    # only for 0 (all) and 1 (first and only population)
    for pop in range(gpr.pops):
        crscale = open(gp.files.get_scale_file(pop),'w')
        print('# Rscale in [pc],',' surfdens_central (=dens0) in [Munit/rscale**2],',\
              ' and in [Munit/pc**2],',' and totmass [Munit],',\
              ' and max(v_LOS) in [km/s]', file=crscale)
        print(Rscale, file=crscale)
        crscale.close()

        print('grh_com: output: ', gpr.fileposcenter[pop])
        filepos = open(gpr.fileposcenter[pop],'w')
        print('# x [Rscale]','y [Rscale]','vLOS [km/s]', file=filepos)
        for k in range(ndm):
            print(xnew[k]/Rscale, ynew[k]/Rscale, vznew[k], file=filepos)
        filepos.close()
        print('')

        
if __name__=='__main__':
    gpr.showplots = True
    import gl_params
    gp = gl_params.Params()
    run(gp)
