#!/usr/bin/env ipython3

##
# @file
# calculate approximative center of mass, assuming constant stellar mass
# 2D version for Fornax, Sculptor, Sextans, .. observed dwarfs

# (c) 2013 Pascal S.P. Steger

import numpy as np
import sys, pdb
import gr_params as gpr
import gl_helper as gh
import gl_file as gf
from gl_centering import com_shrinkcircle_v_2D

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
    gpr.fil = gpr.dir+"/table_merged.bin"
    delim = [0,22,3,3,6,4,3,5,6,6,7,5,6,5,6,5,6]
    ID = np.genfromtxt(gpr.fil, skiprows=29, unpack=True,\
                       usecols=(0,1),delimiter=delim)

    RAh,RAm,RAs,DEd,DEm,DEs,Vmag,VI,\
      VHel,e_VHel,SigFe,e_SigFe,\
      SigMg,e_SigMg,PM = np.genfromtxt(gpr.fil, skiprows=29, unpack=True, \
                                       usecols=tuple(range(2,17)), delimiter=delim, filling_values=-1)

    # only use stars which are members of the dwarf
    pm = (PM>=0.95)
    print("fraction of members = ",1.*sum(pm)/len(pm))
    ID=ID[1][pm]; RAh=RAh[pm]; RAm=RAm[pm]; RAs=RAs[pm]; DEd=DEd[pm]; DEm=DEm[pm]; DEs=DEs[pm];
    Vmag = Vmag[pm]; VI=VI[pm]; VHel=VHel[pm]; e_VHel=e_VHel[pm];
    SigFe=SigFe[pm]; e_SigFe=e_SigFe[pm]; SigMg=SigMg[pm]; e_SigMg=e_SigMg[pm];PM=PM[pm]

    Mg0 = SigMg
    sig = abs(RAh[0])/RAh[0]
    print('RAh: signum = ',sig)
    RAh = RAh/sig
    xs = 15*(RAh*3600+RAm*60+RAs)*sig       # [arcsec/15]

    sig = abs(DEd[0])/DEd[0]
    print('DEd: signum = ',sig)
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

    PM0 = np.copy(PM); x0 = np.copy(xs); y0 = np.copy(ys) # [pc]
    vz0 = np.copy(VHel) # [km/s]

    # only use stars which are members of the dwarf: exclude pop3 by construction
    pm = (PM0 >= gpr.pmsplit) # exclude foreground contamination, outliers
    x0, y0, vz0, Mg0, PM0 = select_pm(x0, y0, vz0, Mg0, PM0, pm)

    # assign population
    if gp.pops == 2:
        # drawing of populations based on metallicity
        # get parameters from function in pymcmetal.py
        [p, mu1, sig1, mu2, sig2] = np.loadtxt(gp.files.dir+'metalsplit.dat')
        [pm1, pm2] = np.loadtxt(gp.files.dir+'metalsplit_assignment.dat')
        pm1 = (pm1>0)
        pm2 = (pm2>0)

    elif gp.pops == 1:
        pm1 = (PM >= 0)
        pm2 = (PM <  0) # assign none, but of same length as xs

    x1, y1, vz1, Mg1, PM1 = select_pm(x0, y0, vz0, Mg0, PM0, pm1)
    x2, y2, vz2, Mg2, PM2 = select_pm(x0, y0, vz0, Mg0, PM0, pm2)

    # cutting pm_i to a maximum of ntracers_i particles each:
    from random import shuffle

    ind1 = np.arange(len(x1))
    np.random.shuffle(ind1)     # random.shuffle already changes ind
    ind1 = ind1[:gp.ntracer[1-1]]

    ind2 = np.arange(len(x2))
    np.random.shuffle(ind2)     # random.shuffle already changes ind
    ind2 = ind2[:gp.ntracer[2-1]]

    x1, y1, vz1, Mg1, PMS1 = select_pm(x1, y1, vz1, Mg1, PM1, ind1)
    x2, y2, vz2, Mg2, PMS2 = select_pm(x2, y2, vz2, Mg2, PM2, ind2)

    x0, y0, vz0, pm1, pm2, pm = concat_pops(x1, x2, y1, y2, vz1, vz2, gp)

    # optimum: get 3D center of mass with means
    # com_x, com_y, com_z = com_mean(x0,y0,z0,PM0) # 3*[pc],  z component included if available

    com_x, com_y, com_vz = com_shrinkcircle_v_2D(x0, y0, vz0, pm) # [pc], [km/s]

    # from now on, work with 2D data only; z0 was only used to get center in (x,y) better
    # x0 -= com_x; y0 -= com_y # [pc]
    # vz0 -= com_vz #[km/s]

    R0 = np.sqrt(x0**2+y0**2) # [pc]
    Rhalf = np.median(R0) # [pc]
    Rscale = Rhalf # [pc] overall

    pop = -1
    for pmn in [pm, pm1, pm2]:
        pop = pop+1
        pmr = (R0<(gp.maxR*Rscale)) # read max extension for data
                                    # (rprior*Rscale) from gl_params
        pmn = pmn*pmr                   # [1]
        print("fraction of members = ", 1.0*sum(pmn)/len(pmn))

        x, y, vz, Mg, PMN = select_pm(x0, y0, vz0, Mg0, PM0, pmn)

        m = np.ones(len(pmn))
        R = np.sqrt(x*x+y*y)            # [pc]
        Rscalei = np.median(R)          # [pc]
        gf.write_Xscale(gp.files.get_scale_file(pop), Rscalei) # [pc]
        gf.write_data_output(gpr.get_com_file(pop), x/Rscalei, y/Rscalei, vz, Rscalei) # [pc]

        if gpr.showplots:
            gpr.show_part_pos(x, y, pmn, Rscale)

if __name__=='__main__':
    # for debugging input issues here:
    gpr.showplots = True
    import gl_params
    gp = gl_params.Params()
    # instead of run(gp), use the code directly, such that after execution, the variables are available

    gpr.fil = gpr.dir+"/table_merged.bin"
    delim = [0,22,3,3,6,4,3,5,6,6,7,5,6,5,6,5,6]
    ID = np.genfromtxt(gpr.fil, skiprows=29, unpack=True,\
                       usecols=(0,1),delimiter=delim)

    RAh,RAm,RAs,DEd,DEm,DEs,Vmag,VI,\
      VHel,e_VHel,SigFe,e_SigFe,\
      SigMg,e_SigMg,PM = np.genfromtxt(gpr.fil, skiprows=29, unpack=True, \
                                       usecols=tuple(range(2,17)), delimiter=delim, filling_values=-1)

    # only use stars which are members of the dwarf
    pm = (PM>=0.95)
    print("fraction of members = ", 1.*sum(pm)/len(pm))
    ID=ID[1][pm]; RAh=RAh[pm]; RAm=RAm[pm]; RAs=RAs[pm]; DEd=DEd[pm]; DEm=DEm[pm]; DEs=DEs[pm];
    Vmag = Vmag[pm]; VI=VI[pm]; VHel=VHel[pm]; e_VHel=e_VHel[pm];
    SigFe=SigFe[pm]; e_SigFe=e_SigFe[pm]; SigMg=SigMg[pm]; e_SigMg=e_SigMg[pm];PM=PM[pm]

    Mg0 = SigMg
    sig = abs(RAh[0])/RAh[0]
    print('RAh: signum = ',sig)
    RAh = RAh/sig
    xs = 15*(RAh*3600+RAm*60+RAs)*sig       # [arcsec] (use 360 deg/12 hrs)

    sig = abs(DEd[0])/DEd[0]
    print('DEd: signum = ',sig)
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

    PM0 = np.copy(PM); x0 = np.copy(xs); y0 = np.copy(ys)
    vz0 = np.copy(VHel)

    # only use stars which are members of the dwarf: exclude pop3 by construction
    pm = (PM0 >= gpr.pmsplit) # exclude foreground contamination, outliers
    x0, y0, vz0, Mg0, PM0 = select_pm(x0, y0, vz0, Mg0, PM0, pm)

    # assign population (OLD, new way is to run grd_split after grd_COM in gl_file
    # if gp.pops == 2:
    #     import pymcmetal as pmc
    #     p, mu1, sig1, mu2, sig2, M = pmc.bimodal_gauss(Mg0)
    #     pm1, pm2 = pmc.assign_pop(Mg0, p, mu1, sig1, mu2, sig2)
    #     fi = open(gp.files.dir+'metalsplit.dat', 'w')
    #     fi.write(str(p)+'\n')
    #     fi.write(str(mu1)+'\n')
    #     fi.write(str(sig1)+'\n')
    #     fi.write(str(mu2)+'\n')
    #     fi.write(str(sig2)+'\n')
    #     fi.close()
    #     np.savetxt(gp.files.dir+'metalsplit_assignment.dat', np.array([pm1, pm2]))
