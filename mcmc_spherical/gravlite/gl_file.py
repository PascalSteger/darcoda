#!/usr/bin/python2.7
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


    if gp.bprior and gp.investigate != 'simple':
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






def write_outfile():
    '''write profiles to output files in directory'''
    M = phys.Mzdefault(gp.pars.dens)
    profM, profnus, profdeltas, profsigs = gp.files.get_outprofs()
    arraydump(profM, M)
    arraydump(profnus[0],   gp.pars.nu1)
    arraydump(profdeltas[0], gp.pars.delta1)
    arraydump(profsigs[0],  gp.sig1_x)
    if gp.pops == 2:
        arraydump(profnus[1],   gp.pars.nu2)
        arraydump(profdeltas[1], gp.pars.delta2)
        arraydump(profsigs[1],  gp.sig2_x)
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






def adump():
    write_key_data_parameters()
    arraydump(gp.files.get_outdat(), gp.xipol, 'w')

    profM, profnus, profdeltas, profsigs = gp.files.get_outprofs()
    arraydump(profM, gp.xipol, 'w')
    arraydump(profnus[0], gp.xipol, 'w')
    arraydump(profdeltas[0], gp.xipol, 'w')
    arraydump(profsigs[0], gp.xipol, 'w')
    if gp.pops==2:
        arraydump(profnus[1], gp.xipol, 'w')
        arraydump(profdeltas[1], gp.xipol, 'w')
        arraydump(profsigs[1], gp.xipol, 'w')
    return 0




    

def store_old_params(pars,chisq):
    gp.run_configs.append((pars,chisq))
    return gp.run_configs






def store_working_pars(n,pars,chisq,parstep):
    gp.init_configs.append([pars,chisq,parstep])
    twelve = open( gp.files.get_outtxt(), 'a')
    print>>twelve, n, chisq
    twelve.close()
    return gp.init_configs





    
def get_working_pars():
    if len(gp.init_configs)<1:
        gp.parst  = gp.pars
        gp.chisqt = gp.chisq
        return
    else:
        gp.parst,gp.chisqt,gp.parstep = gp.init_configs.pop()
        gp.parstep.adaptworst(1./gp.stepafterrunaway)
    return gp.parst
