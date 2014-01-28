#!/usr/bin/env python3

##
# @file
# all file related functions

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import pdb
import sys
import numpy as np
import gl_params as gp
import gl_physics as phys

from gl_data import Datafile

def bin_data():
    if gp.investigate == 'hern':
        import grh_com
        import grh_Pos
        import grh_MCMCbin
        grh_MCMCbin.run()
    elif gp.investigate == 'gaia':
        import grg_COM, grg_MCMCbin
        grg_COM.run()
        grg_MCMCbin.run()
    elif gp.investigate == 'walk':
        # TODO: call main again after first iteration, if gp.metalpop set
        import grw_COM, grw_MCMCbin # inside there, split by metallicity
        grw_COM.run()
        grw_MCMCbin.run()
        # run for 3D models as well if model is set (needed in rhowalktot)
        if gp.model:
            import grw_com, grw_mcmcbin
            grw_com.run()
            grw_mcmcbin.run()
    elif gp.investigate == 'discsim':
        import grs_com_align # centering, if not aligned yet
        import grs_rho
        import grs_siglos
    return
## \fn bin_data()
# get data, bin it anew (e.g. if gp.nbin changed)


def get_data():
    gp.dat = Datafile()
    if gp.investigate == 'discmock':
        import gl_discmock as gs
        gs.discmock()
    elif gp.investigate == 'discsim':
        import gl_discsim as gs
        gs.discsim()
    else: # for all dwarfs, read from files
        for i in range(gp.pops+1):
            A = np.loadtxt(gp.files.get_scale_file(i), unpack=False, skiprows=1)
            gp.Rscale.append(A[0]) # TODO: error in Gaia case: sometimes, we do not have entries here
            gp.Nu0rscale.append(A[1])
            gp.Nu0pc.append(A[2])
            gp.totmass.append(A[3])
            gp.maxvlos.append(A[4])

        gp.rstarhalf = gp.Rscale[0]
        gp.dat.read_nu()    # set gp.xipol in here
        gp.dat.read_sigma()
        gp.dat.read_kappa()
    return gp.dat
## \fn get_data()
# read in data, store in a gl_data class


def arraydump(fname,arrays,app='a',narr=1):
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
## \fn arraydump(fname, arrays, app='a', narr=1)
# This routine takes a number, narr, of equal length arrays
# and appends/writes them to a specified file in columnated data format.
# @param fname filename, string
# @param arrays =[[arr1], [arr2]...] 
# @param app='a' appending?
# @param narr=1 number of arrays


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
