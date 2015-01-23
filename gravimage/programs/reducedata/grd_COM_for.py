#!/usr/bin/env ipython3

##
# @file
# calculate approximative center of mass, assuming constant stellar mass
# 2D version for Fornax only, based on deBoer+2013 data

# (c) GPL v3 2015 Pascal S.P. Steger

import numpy as np
import pdb
import gi_helper as gh
import gi_file as gf
from gi_centering import com_shrinkcircle_2D

def concat_pops(x1, x2, y1, y2, vz1, vz2, gp):
    x0 = np.hstack([x1, x2])
    y0 = np.hstack([y1, y2])
    vz0 = np.hstack([vz1, vz2])
    N1  = min(len(x1), gp.ntracer[1-1])
    N2  = min(len(x2), gp.ntracer[2-1])
    pm1 = np.hstack([np.ones(N1, dtype=bool), np.zeros(N2, dtype=bool)])
    pm2 = np.hstack([np.zeros(N1, dtype=bool), np.ones(N2,  dtype=bool)])
    pm = pm1 + pm2
    return x0, y0, vz0, pm1, pm2, pm
## \fn concat_pops(x1, x2, y1, y2, z1, z2, vz1, vz2)
# concatenate all arrays for two populations
# @param x1
# @param x2
# @param y1
# @param y2
# @param vz1
# @param vz2

def select_pm(x, y, vz, Mg, PM, pm):
    return x[pm], y[pm], vz[pm], Mg[pm], PM[pm]
## \fn select_pm(x, y, comp, vz, Mg, PM, pm)
# extract only parts of the arrays given
# @param x
# @param y
# @param vz
# @param Mg
# @param PM
# @param pm

def run(gp):
    import gr_params
    gpr = gr_params.grParams(gp)
    gpr.fil = gpr.dir+"/deBoer/table1.dat"
    ALL = np.loadtxt(gpr.fil)
    RAh = ALL[:,0]
    RAm = ALL[:,1]
    RAs = ALL[:,2]
    DEd = ALL[:,3]
    DEm = ALL[:,4]
    DEs = ALL[:,5]
    # that's all we read in for now. Crude assumptions: each star belongs to Fornax, and has mass 1Msun

    # only use stars which are members of the dwarf
    sig = abs(RAh[0])/RAh[0]
    RAh = RAh/sig
    xs = 15*(RAh*3600+RAm*60+RAs)*sig       # [arcsec/15]
    sig = abs(DEd[0])/DEd[0]
    DEd = DEd/sig
    ys = (DEd*3600+DEm*60+DEs)*sig          # [arcsec]
    arcsec = 2.*np.pi/(360.*60.*60) # [pc]
    kpc = 1000 # [pc]
    DL = {0: lambda x: x * (138),#+/- 8 for Fornax
          1: lambda x: x * (101),#+/- 5 for Carina
          2: lambda x: x * (79),  #+/- 4 for Sculptor
          3: lambda x: x * (86) #+/- 4 for Sextans
      }[gp.case](kpc)

    xs *= (arcsec*DL) # [pc]
    ys *= (arcsec*DL) # [pc]
    x0 = np.copy(xs)
    y0 = np.copy(ys) # [pc]
    com_x, com_y = com_shrinkcircle_2D(x0, y0) # [pc], [km/s]
    # from now on, work with 2D data only; z0 was only used to get center in (x,y) better
    # x0 -= com_x; y0 -= com_y # [pc]
    # vz0 -= com_vz #[km/s]
    R0 = np.sqrt(x0**2+y0**2) # [pc]
    Rhalf = np.median(R0) # [pc]
    Rscale = Rhalf # [pc] overall
    pop = 0
    pmr = (R0<(gp.maxR*Rscale)) # read max extension for data (rprior*Rscale) from gi_params
    x=1.*x0[pmr]
    y=1.*y0[pmr]
    R = np.sqrt(x*x+y*y)            # [pc]
    Rscalei = np.median(R)          # [pc]
    gf.write_Xscale(gp.files.get_scale_file(pop), Rscalei) # [pc]
    gf.write_data_output(gp.files.get_com_file(pop), x/Rscalei, y/Rscalei, np.zeros(len(x)), Rscalei) # [pc]
## \fn run(gp)
# main functionality: get center of mass of deBoer data
# @param gp global parameters






if __name__=='__main__':
    # for debugging input issues here:
    import gi_params
    gp = gi_params.Params()
    # instead of run(gp), use the code directly, such that after execution, the variables are available
