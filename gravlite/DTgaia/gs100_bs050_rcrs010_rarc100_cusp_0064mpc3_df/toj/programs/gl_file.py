#!/usr/bin/env ipython3

##
# @file
# all file related functions

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import sys, pdb
import numpy as np
import gl_physics as phys
from gl_data import Datafile

def bin_data(gp):
    if gp.investigate == 'hern':
        import grh_com, grh_Pos, grh_MCMCbin
        grh_com.run(gp)
        grh_Pos.run()
        grh_MCMCbin.run(gp)
        grh_MCMCbin.run(gp)
    elif gp.investigate == 'gaia':
        import grg_COM, grg_MCMCbin
        grg_COM.run(gp)
        grg_MCMCbin.run(gp)
    elif gp.investigate == 'walk':
        import grw_COM, grw_MCMCbin # inside there, split by metallicity
        grw_COM.run(gp)
        grw_MCMCbin.run(gp)
        # run for 3D models as well if model is set (needed in rhotot_walk)
        if gp.walker3D:
            import grw_com, grw_mcmcbin
            grw_com.run()
            grw_mcmcbin.run()
    elif gp.investigate == 'triax':
        import grt_com
        grt_com.run(gp)
        import grt_dens
        grt_dens.run(gp)
        import grt_siglos
        grt_siglos.run(gp)
    elif gp.investigate == 'obs':
        import grd_COM, grd_MCMCbin
        grd_COM.run(gp)
        grd_MCMCbin.run(gp)
    elif gp.investigate == 'discmock':
        import grdm_write
        grdm_write.run(gp)
    elif gp.investigate == 'discsim':
        import grds_write
        grds_write.run(gp)
    return
## \fn bin_data(gp)
# get data, bin it anew (e.g. if gp.nbin changed)
# @param gp global parameter


def get_data(gp):
    gp.dat = Datafile()
    for pop in range(gp.pops+1):
        A = np.loadtxt(gp.files.get_scale_file(pop), unpack=False, skiprows=1)
        gp.Xscale.append(A[0])
        gp.Sig0pc.append(A[1])
        gp.totmass.append(A[2])
        gp.maxsiglos.append(A[3])
    #if gp.investigate == "walk":
    #    for pop in range(gp.pops):
    #        gp.ntracer[pop] = gp.totmass[pop+1]

    gp.dat.read_Sig(gp)    # set gp.xipol in here
    gp.dat.read_sig(gp)
    if gp.usekappa:
        gp.dat.read_kappa(gp)
    return gp.dat
## \fn get_data(gp)
# read in data, store in a gl_data class
# @param gp global parameters


def get_rhohalfs(gp):
    if gp.geom == 'sphere':
        # Wolf, Walker method for M_half, r_half
        # assuming isotropic Plummer profile:
        r_half = gp.dat.rhalf[0] # [pc] from overall rho*
        sigv = np.median(gp.dat.sig[0]) # [km/s] from overall rho*
        M_half = 5.*r_half*sigv**2/(2.*gp.G1) # [Munit] Walker Penarrubia 2011
        # other estimate: Wolf+2010,
        # M(4/3 r_half = R_half) = 4 r_half sigv**2/gp.G
        gp.rhohalf = M_half/(4.*np.pi*r_half**3/3.)
    elif gp.geom == 'disc':
        return -1.08
## \fn get_rhohalfs(gp)
# get informed priors on 3D densities at half-light radius
# via deprojection for nu
# and M_half = M(<r_half) for rho
# @param gp global parameters


def arraydump(fname, arrays, app='a', narr=1):
    fn=open(fname,app)
    if narr == 1:
        print(" ".join(map(str,arrays)), file=fn)
    else:
        anew = np.transpose(arrays)
        for line in anew:
            if(isinstance(line,list)):
                print(" ".join(map(str, line)), file=fn)
            else:
                print(line, file=fn)
    fn.close()
    return 0
## \fn arraydump(fname, arrays, app, narr)
# This routine takes a number, narr, of equal length arrays
# and appends/writes them to a specified file in columnated data format.
# @param fname filename, string
# @param arrays =[[arr1], [arr2]...] 
# @param app  ='a' appending?
# @param narr =1 number of arrays


def bufcount(filename):
    f = open(filename)
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read # loop optimization
    
    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)
        
    return lines
## \fn bufcount(filename)
# count lines of a file
# @param filename string


def write_headers_2D(gp, pop):
    de = open(gp.files.Sigfiles[pop], 'w')
    print('Rbin [Xscale];','Binmin [Xscale];','Binmax [Xscale];',\
          'Sig(R) [Munit/pc^2];','error [1]', file=de)
    
    em = open(gp.files.massfiles[pop],'w')
    print('R [Xscale];','Binmin [Xscale];','Binmax [Xscale];',\
          'M(<Binmax) [Munit];','error [Munit]', file=em)
    
    sigfil = open(gp.files.sigfiles[pop],'w')
    print('R [Xscale];','Binmin [Xscale];','Binmax [Xscale];',\
          'sigma_p(R) [maxsiglos];','error [km/s]', file=sigfil)

    kappafil = open(gp.files.kappafiles[pop],'w')
    print('R [Xscale];','Binmin [Xscale];','Binmax [Xscale];',\
          'kappa_los(R) [1];','error [1]', file=kappafil)
    return de, em, sigfil, kappafil
## \fn write_headers_2D(gp, pop)
# write headers for datareduction output files, and return file handlers
# @param gp global parameters
# @param pop component (0: all, 1,2,...)


def write_headers_3D(gp, pop):
    de = open(gp.files.Sigfiles[pop]+'_3D', 'w')
    print('rbin [rscale];','binmin [rscale];','binmax [rscale];','nu(r) [nu(0)];', 'error', file=de)
    
    em = open(gp.files.massfiles[pop]+'_3D','w')
    print('rbin [rscale];','binmin [rscale];','binmax [rscale];','M(<r) [Munit];','error', file=em)
    return de, em
## \fn write_headers_3D(gp, pop)
# write headers for the 3D quantity output
# @param gp global parameters
# @param pop population int


def empty(filename):
    if bufcount(filename)<2:
        return True
    return False
## \fn empty(filename)
# determine if file is empty or not
# @param filename string


def read_Xscale(filename):
    crscale = open(filename, 'r')
    Xscale = np.loadtxt(crscale, comments='#', skiprows=1, unpack=False)
    crscale.close()
    return Xscale
## \fn read_Xscale(filename)
# read scale radius from file
# @param filename string


def write_tracer_file(filename, totmass):
    tr = open(filename, 'w')
    print(totmass, file=tr)
    tr.close()
## \fn write_tracer_file(filename, totmass)
# write tracer file
# @param filename
# @param totmass
    

def write_density(filename, Dens0pc, totmass):
    cdens = open(filename, 'a')
    print(Dens0pc, file=cdens)                      # [Munit/pc^2]
    print(totmass, file=cdens)                      # [Munit]
    cdens.close()
## \fn write_density(filename, Dens0pc, totmass)
# output density
# @param filename
# @param Dens0pc central density [Munit/pc^2]
# @param totmass
    

def write_data_output(filename, x, y, vz, Xscale):
    print('output: ', filename)
    c = open(filename,'w')
    print('# x [Xscale], y [Xscale], vLOS [km/s], Xscale = ',Xscale,' pc', file=c)
    for k in range(len(x)):
        print(x[k], y[k], vz[k], file=c)      # [Xscale], [Xscale], [km/s]
    c.close()
    return
## \fn write_data_output(filename, x, y, vz, Xscale)
# write x,y,vz to files
# @param filename
# @param x
# @param y
# @param vz
# @param Xscale


def write_Xscale(filename, Xscale):
    crscale = open(filename, 'w')
    print('# Xscale in [pc], surfdens_central (=dens0) in [Munit/Xscale^2], and in [Munit/pc^2], and totmass [Munit], and max(v_LOS) in [km/s]', file=crscale)
    print(Xscale, file=crscale)
    crscale.close()
    return
## \fn write_Xscale(filename, Xscale)
# store central values in scale_ file
# @param filename string
# @param Xscale float, [pc]

