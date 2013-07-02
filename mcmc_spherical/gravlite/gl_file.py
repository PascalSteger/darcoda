#!/usr/bin/python
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
'''all file related functions'''

import gl_params as gp
import pdb
import sys
import logging
logging.basicConfig(stream=sys.stderr, level=logging.INFO)
LOG = logging.getLogger(__name__)
from gl_data import *
if gp.geom == 'sphere':
    import physics_sphere as phys
elif gp.geom == 'disc':
    import physics_disc as phys




def bin_data():
    if gp.investigate == 'hernquist':
        import grh_com
        import grh_pos2D
        import grh_dens2D
        import grh_siglos2D
    elif gp.investigate == 'walker':
        # TODO: call main again after first iteration, if gp.metalpop set
        import grw_com                  # inside there, split by metallicity
        grw_com.run()
        import grw_dens
        grw_dens.run()
        import grw_siglos
        grw_siglos.run()
    elif gp.investigate == 'sim':
        import grs_com_align # centering, if not aligned yet
        import grs_dens
        import grs_siglos
        




def get_data():

    gp.dat = Datafile()
    if gp.investigate == 'simple':
        import gl_disc_simple as gs
        gs.disc_simple()
    elif gp.investigate == 'sim':
        import gl_disc_sim as gs
        gs.disc_sim()
    else: # for all dwarfs, read from files
        if gp.investigate == 'walker':
            for i in range(3):
                A = np.loadtxt(gp.files.get_scale_file(i), unpack=False, skiprows=1)
                gp.rcore_2D.append(A[0])
                gp.dens0rcore_2D.append(A[1])
                gp.dens0pc_2D.append(A[2])
                gp.totmass.append(A[3])
                gp.maxvlos.append(A[4])

        gp.dat.read_mass()
        gp.dat.read_nu()
        gp.dat.read_sigma()


    if gp.bprior:
        gp.blow = gp.dat.Mdat - gp.dat.Merr


    # Binning in z:
    if (gp.xpmin<0) : gp.xpmin = min(gp.dat.Mx)
    if (gp.xpmax<0) : gp.xpmax = max(gp.dat.Mx)
    return gp.dat





def ipol_data():
    '''interpolate all data to nipol bins with same range of r (called ripol)'''
    gp.ipol = Datafile()
    gp.ipol.interpol(gp.dat)
    return gp.ipol






def write_key_data_parameters():
    twelve=open(gp.files.get_outtxt(),'w')
    print>>twelve,'Number of terms [M,nu,delta]  & iterations:'
    print>>twelve, gp.nipol, gp.niter
    print>>twelve,'Run parameters [gprior, cprior, bprior]:'
    print>>twelve, gp.gprior, gp.cprior, gp.bprior
    twelve.close()
    return 0





def adump():
    write_key_data_parameters()
    arraydump(gp.files.get_outdat(), gp.xipol, 'w')

    profM, profdens, profnus, profdeltas, profsigs = gp.files.get_outprofs()
    arraydump(profM, gp.xipol, 'w')
    arraydump(profdens, gp.xipol, 'w')
    arraydump(profnus[0], gp.xipol, 'w')
    arraydump(profdeltas[0], gp.xipol, 'w')
    arraydump(profsigs[0], gp.xipol, 'w')
    if gp.pops==2:
        arraydump(profnus[1], gp.xipol, 'w')
        arraydump(profdeltas[1], gp.xipol, 'w')
        arraydump(profsigs[1], gp.xipol, 'w')
    return 0






def write_outfile():
    '''write profiles to output files in directory'''
    if gp.geom == 'sphere':
        M = phys.Mzdefault(gp.pars.dens)
    else:
        M = gp.M_x  # TODO: check meaning, possible inclusion in physics_disc.Mzdefault
        
    profM, profdens, profnus, profdeltas, profsigs = gp.files.get_outprofs()
    arraydump(profM, M)
    arraydump(profdens, gp.dens_x) # [Msun/pc^3] in spherical case
    arraydump(profnus[0],   phys.nu(gp.pars.nu1))  # [Msun/pc^3]
    arraydump(profdeltas[0], gp.pars.delta1)       # [1]
    arraydump(profsigs[0],  gp.sig1_x)             # [km/s]
    if gp.pops == 2:
        arraydump(profnus[1],   phys.nu(gp.pars.nu2)) # [Msun/pc^3]
        arraydump(profdeltas[1], gp.pars.delta2)      # [1]
        arraydump(profsigs[1],  gp.sig2_x)            # [km/s]
    return 0






def arraydump(fname,arrays,app='a',narr=1):
    '''This routine takes a number, narr, of equal length arrays:
    arrays=[[arr1], [arr2]...] and appends/writes them to a specified
    file (fname) in columnated data format.'''
    fn=open(fname,app)
    if narr == 1:
        print >> fn," ".join(map(str,arrays))
    else:
        anew = np.transpose(arrays)
        for line in anew:
            if(isinstance(line,list)):
                print>>fn," ".join(map(str, line))
            else:
                print>>fn,line
    fn.close()
    return 0










    

def store_old_params(pars,chi2):
    gp.run_configs.append((pars,chi2))
    return gp.run_configs






def store_working_pars(n,pars,chi2,parstep):
    gp.init_configs.append([pars,chi2,parstep])
    twelve = open( gp.files.get_outtxt(), 'a')
    print>>twelve, n, chi2
    twelve.close()
    return gp.init_configs





    
def get_working_pars(scale=True):
    if len(gp.init_configs)<1:
        gp.pars.assign(gp.safepars); gp.parstep.assign(gp.safeparstep)
        gp.chi2 = gp.safechi2
        return
    else:
        gp.pars, gp.chi2, gp.parstep = gp.init_configs.pop()
        gp.parst.assign(gp.pars)
        if scale: gp.parstep.adaptworst(gp.stepafterrunaway)
    return gp.pars




def bufcount(filename):
    '''count lines of a file'''
    f = open(filename)
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read # loop optimization
    
    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)
        
    return lines
