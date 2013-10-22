#!/usr/bin/env ipython-python3.2

##
# @file
# all file related functions
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch


import pdb
import sys
import gl_params as gp
if gp.geom == 'sphere':
    import physics_sphere as phys
else:
    import physics_disc as phys

from gl_data import *


## get data, bin it anew (e.g. if gp.nbin changed)
def bin_data():
    if gp.investigate == 'hernquist':
        import grh_com
        import grh_Pos
        import grh_MCMCbin
        grh_MCMCbin.run()
    elif gp.investigate == 'gaia':
        import grg_COM, grg_MCMCbin
        grg_COM.run()
        grg_MCMCbin.run()
    elif gp.investigate == 'walker':
        # TODO: call main again after first iteration, if gp.metalpop set
        import grw_COM, grw_MCMCbin # inside there, split by metallicity
        grw_COM.run()
        grw_MCMCbin.run()
        # run for 3D models as well if model is set (needed in rhowalkertot)
        if gp.model:
            import grw_com, grw_mcmcbin
            grw_com.run()
            grw_mcmcbin.run()
    elif gp.investigate == 'sim':
        import grs_com_align # centering, if not aligned yet
        import grs_dens
        import grs_siglos


## read in data, store in a gl_data class
def get_data():
    gp.dat = Datafile()
    if gp.investigate == 'simple':
        import gl_disc_simple as gs
        gs.disc_simple()
    elif gp.investigate == 'sim':
        import gl_disc_sim as gs
        gs.disc_sim()
    else: # for all dwarfs, read from files
        if gp.investigate == 'walker' or gp.investigate == 'triaxial' or gp.investigate == 'gaia':
            for i in range(gp.pops+1):
                A = np.loadtxt(gp.files.get_scale_file(i), unpack=False, skiprows=1)
                gp.rcore_2D.append(A[0])
                gp.dens0rcore_2D.append(A[1])
                gp.dens0pc_2D.append(A[2])
                gp.totmass.append(A[3])
                gp.maxvlos.append(A[4])

        gp.dat.read_mass()
        gp.dat.read_nu()
        gp.dat.read_sigma()
        gp.dat.read_kappa()

    if gp.bprior:
        gp.blow = gp.dat.Mdat - gp.dat.Merr


    # Binning in z:
    if (gp.xpmin<0) : gp.xpmin = min(gp.dat.Mx)
    if (gp.xpmax<0) : gp.xpmax = max(gp.dat.Mx)
    return gp.dat

## interpolate all data to nipol bins with same range of r (called ripol)
# @return interpolated values in gl_data class
def ipol_data():
    gp.ipol = Datafile()
    if gp.consttr:                      # if set const tracer number, do NOT interpolate
        gp.ipol.copyfrom(gp.dat)
    else:
        gp.ipol.interpol(gp.dat)
    return gp.ipol


## write MCMC characteristics to .txt file
def write_key_data_parameters():
    twelve=open(gp.files.get_outtxt(),'w')
    print('Number of terms [M,nu,delta]  & iterations:', file=twelve)
    print(gp.nipol, gp.niter, file=twelve)
    print('Run parameters [gprior, cprior, bprior]:', file=twelve)
    print(gp.gprior, gp.cprior, gp.bprior, file=twelve)
    twelve.close()
    return 0

## write all MCMC parameters at the end of a run
def adump():
    write_key_data_parameters()
    arraydump(gp.files.get_outdat(), gp.xipol, 'w')

    profM, profdens, profnus, profdeltas, profsigs, profkaps = gp.files.get_outprofs()
    arraydump(profM, gp.xipol, 'w')
    arraydump(profdens, gp.xipol, 'w')
    arraydump(profnus[0], gp.xipol, 'w')
    arraydump(profdeltas[0], gp.xipol, 'w')
    arraydump(profsigs[0], gp.xipol, 'w')
    arraydump(profkaps[0], gp.xipol, 'w')
    if gp.pops==2:
        arraydump(profnus[1], gp.xipol, 'w')
        arraydump(profdeltas[1], gp.xipol, 'w')
        arraydump(profsigs[1], gp.xipol, 'w')
        arraydump(profkaps[1], gp.xipol, 'w')
    return 0

## write profiles to output files in directory
def write_outfile():
    profM, profdens, profnus, profdeltas, profsigs, profkaps = gp.files.get_outprofs()
    arraydump(profM,         gp.M_x)
    arraydump(profdens,      gp.dens_x) # [Msun/pc^3] in spherical case
    arraydump(profnus[0],    gp.nu1_x)  # [Msun/pc^3]
    arraydump(profdeltas[0], gp.d1_x)   # [1]
    arraydump(profsigs[0],   gp.sig1_x) # [km/s]
    arraydump(profkaps[0],   gp.kap1_x) # [km/s]
    if gp.pops == 2:
        arraydump(profnus[1],    gp.nu2_x)  # [Msun/pc^3]
        arraydump(profdeltas[1], gp.d2_x)   # [1]
        arraydump(profsigs[1],   gp.sig2_x) # [km/s]
        arraydump(profkaps[1],   gp.kap2_x) # [km/s]
    return 0

## This routine takes a number, narr, of equal length arrays:
#  arrays=[[arr1], [arr2]...] and appends/writes them to a specified
#  file (fname) in columnated data format.
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

## store parameters into an array of 'good' parameters
def store_old_params(pars,chi2):
    gp.run_configs.append((pars,chi2))
    return gp.run_configs

## write working parameters to file .txt if finished with init phase
def store_working_pars(n,pars,chi2,parstep):
    gp.init_configs.append([pars,chi2,parstep])
    if not gp.initphase:
        twelve = open( gp.files.get_outtxt(), 'a')
        print(n, chi2, file=twelve)
        twelve.close()
    return gp.init_configs
    

## restore 'good' parameters if any exception happened
# @return parameters in class
def get_working_pars(scale=False):
    if len(gp.init_configs)<1:
        gp.pars.assign(gp.safepars); gp.parstep.assign(gp.safeparstep)
        gp.chi2 = gp.safechi2
        return gp.pars
    else:
        gp.pars, gp.chi2, gp.parstep = gp.init_configs.pop()
        gp.parst.assign(gp.pars)
        if scale: gp.parstep.adaptworst(gp.stepafterrunaway)
    # TODO: jump to end of main loop, independet of where call originated
    return gp.pars


## count lines of a file
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
