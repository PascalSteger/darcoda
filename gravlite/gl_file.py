#!/usr/bin/env python3

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
        grg_COM.run()
        grg_MCMCbin.run()
    elif gp.investigate == 'walk':
        # TODO: call main again after first iteration, if gp.metalpop set
        import grw_COM, grw_MCMCbin # inside there, split by metallicity
        grw_COM.run(gp)
        grw_MCMCbin.run(gp)
        # run for 3D models as well if model is set (needed in rhowalktot)
        if gp.model:
            import grw_com, grw_mcmcbin
            grw_com.run()
            grw_mcmcbin.run()
    elif gp.investigate == 'triax':
        import grt_com
        grt_com.run()
        import grt_dens
        grt_dens.run()
        import grt_siglos
        grt_siglos.run()
    elif gp.investigate == 'obs':
        import grd_COM, grd_MCMCbin
        grd_COM.run(gp)
        grd_MCMCbin.run(gp)
    elif gp.investigate == 'discsim':
        import grs_com_align # centering, if not aligned yet
        import grs_rho
        import grs_siglos
    return
## \fn bin_data(gp)
# get data, bin it anew (e.g. if gp.nbin changed)


def get_data(gp):
    gp.dat = Datafile()
    if gp.investigate == 'discmock':
        import gl_disc_mock as gs
        gs.disc_mock(gp)
    elif gp.investigate == 'discsim':
        import gl_disc_sim as gs
        gs.disc_sim(gp)
    else: # for all dwarfs, read from files
        for i in range(gp.pops+1):
            A = np.loadtxt(gp.files.get_scale_file(i), unpack=False, skiprows=1)
            gp.Rscale.append(A[0]) # TODO: error in Gaia case:
                                   # sometimes, we miss
                                   # entries here
            gp.Nu0rscale.append(A[1])
            gp.Nu0pc.append(A[2])
            gp.totmass.append(A[3])
            gp.maxvlos.append(A[4])

        gp.rstarhalf = gp.Rscale[0]
        gp.dat.read_nu(gp)    # set gp.xipol in here
        gp.dat.read_sigma(gp)
        if gp.usekappa:
            gp.dat.read_kappa(gp)
    return gp.dat
## \fn get_data()
# read in data, store in a gl_data class


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

def write_headers(gp, comp):
    de = open(gp.files.nufiles[comp], 'w')
    print('Rbin [Rscale]','Binmin [Rscale]','Binmax [Rscale]',\
          'Nu(R)/Nu(0) [1]','error [1]', file=de)
    
    em = open(gp.files.massfiles[comp],'w')
    print('R [Rscale]','Binmin [Rscale]','Binmax [Rscale]',\
          'M(<Binmax) [Msun]','error [Msun]', file=em)
    
    sigfil = open(gp.files.sigfiles[comp],'w')
    print('R [Rscale]','Binmin [Rscale]','Binmax [Rscale]',\
          'sigma_r(R) [km/s]','error [km/s]', file=sigfil)

    kappafil = open(gp.files.kappafiles[comp],'w')
    print('R [Rscale]','Binmin [Rscale]','Binmax [Rscale]',\
          'kappa_los(R) [1]','error [1]', file=kappafil)
    return de, em, sigfil, kappafil
## \fn write_headers(gp, comp)
# write headers for datareduction output files, and return file handlers
# @param gp global parameters
# @param comp component (0: all, 1,2,...)
