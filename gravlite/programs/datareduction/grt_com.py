#!/usr/bin/env ipython3

##
# @file
# calculate approximative center of mass, assuming constant stellar mass, for triax. system

# (c) 2013 Pascal S.P. Steger, psteger@phys.ethz.ch


import numpy as np
import sys, pdb
from pylab import *
ion()

import gr_params as gpr
from gl_helper import expDtofloat
from gl_class_files import *
from gl_centering import *
import gl_file as gf

def run(gp):
    print('input: ', gpr.fil)
    x0,y0,vlos = np.genfromtxt(gpr.fil, skiprows=0, unpack =  True,
                               usecols = (0,1,5))

    # use only 3000 random particles:
    ind = np.arange(len(x0))
    np.random.shuffle(ind)
    ind = ind[:3000]
    x0 = x0[ind];    y0 = y0[ind];    vlos = vlos[ind]

    x0 *= 1000.                         # [pc]
    y0 *= 1000.                         # [pc]

    # shrinking sphere method
    pm = np.ones(len(x0))
    com_x, com_y, com_vz = com_shrinkcircle_v_2D(x0, y0, vlos, pm)

    x0 -= com_x # [pc]
    y0 -= com_y # [pc]
    vlos -= com_vz #[km/s]

    for pop in range(2):
        Rc = np.sqrt(x0**2+y0**2) # [pc]
        Rhalf = np.median(Rc) # [pc]
        Rscale = Rhalf # or gpr.r_DM # [pc]
        gp.Xscale.append(Rscale) # [pc]

        print('Rscale = ', Rscale,' pc')
        print('max(R) = ', max(Rc),' pc')
        print('total number of stars: ', len(Rc))

        R0 = np.sqrt(x0**2+y0**2)/Rscale
        sel = (R0 < gp.maxR)
        x = x0[sel]/Rscale; y = y0[sel]/Rscale # [Rscale]
        vz = vlos[sel] # [km/s]
        m = np.ones(len(x))
        R = np.sqrt(x*x+y*y)*Rscale # [pc]

        gf.write_Xscale(gp.files.get_scale_file(pop), np.median(R))

        c = open(gpr.get_com_file(pop), 'w')
        print('# x [Xscale],','y [Xscale],','vLOS [km/s],','Xscale = ', \
              Rscale, ' pc', file=c)
        for k in range(len(x)):
            print(x[k], y[k], vz[k], file=c) #[rscale], [rscale], [km/s]
        c.close()

        if gpr.showplots:
            gpr.show_part_pos(x, y, np.ones(len(x)), Rscale)

if __name__=='__main__':
    gpr.showplots = True
    import gl_params
    gp = gl_params.Params()
    run(gp)
