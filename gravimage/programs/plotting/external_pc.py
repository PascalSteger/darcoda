#!/usr/bin/env python3

# include runtime gi_params, but other files all from current directory

import os
import sys
import pdb
import pickle
import numpy as np

basename = sys.argv[1]
tt = sys.argv[2]
minlinelen = int(sys.argv[3])
basedir = basename+tt+'/'

def prepare_output_folder(outdir):
    os.system('mkdir -p '+ outdir + 'output/data/')
    os.system('mkdir -p '+ outdir + 'output/ascii/')
    os.system('mkdir -p '+ outdir + 'output/pdf/')
    os.system('mkdir -p '+ outdir + 'output/analytic/')
    return 0
## \fn prepare_output_folder(outdir)
# create directory structure in output folder
# @param outdir string, data folder

def read_scale(gp):
    gp.Xscale = []
    gp.Sig0pc = []
    gp.maxsiglos = []
    for pop in range(gp.pops+1):
        A = np.loadtxt(gp.files.get_scale_file(pop), unpack=False, skiprows=1)
        gp.Xscale.append(A[0])  # [pc]
        gp.Sig0pc.append(A[1])  # [Munit/pc^2]
        # totmass_tracers is A[2] # [Munit]
        gp.nu0pc.append(A[3])   # [Munit/pc^3]
        gp.maxsiglos.append(A[4]) # [km/s]
## \fn read_scale(gp)
# read scale file, store into gp.*scale
# @param gp

def pcload_single_entries(bn, gp):
    import gi_collection as glc
    pc = glc.ProfileCollection(gp.pops, gp.nepol)
    import re
    tmp = re.split('/DT', bn)[-1]
    path = str.join('/', re.split('/', tmp)[:-1])
    # get number of iterations from run_info file
    import gi_base as gb
    bp = gb.get_basepath()
    fil = open(bp+"/run_info", "r")
    ln = 0
    for line in fil:
        ln = ln +1
        if ln > minlinelen:
            break
        if re.search(path, line):
            line2 = re.sub(r'\n', '', line)
            if not re.search("File ", line2):
                runparams = line2
    fil.close()
    numofmodels = int(re.split('\t', runparams)[2])
    current = 0
    with open(bn+'pc2.save', 'rb') as fi:
        dum = pickle.load(fi) # dummy variable, was used to create file
        try:
            while current<minlinelen: # if needing to cut to max no. iterations
                current += 1
                if current%100 == 0:
                    gh.progressbar((1.0*current)/numofmodels)
                MODEL = pickle.load(fi)
                pc.add(MODEL)
        except EOFError:
            pass
    print("")
    return pc
## \fn pcload_single_entries(bn, gp)
# load all data into [chi^2, profiles] pairs, and add to a profile collection
# @param bn string
# @param gp global parameters

import import_path as ip
ip.insert_sys_path(basedir+'programs/')
ip.insert_sys_path(basedir+'programs/sphere')
import gi_params as ngip
ngp = ngip.Params(tt)
print(ngp)
print('ngp.rinfty = ',ngp.rinfty)
import select_run as sr
ngp.pops = sr.get_pops(basedir)
print('working with ', ngp.pops, ' populations')
prepare_output_folder(basedir)

# check whether we need to read in ngp.dat, or whether we are plotting from inside gravimage main program
if len(ngp.dat.Sig) == 0:
    import gi_file as glf
    ngp.dat = glf.get_binned_data(ngp)
read_scale(ngp) # store half-light radii in  gp.Xscale
import gi_helper as gh
Radii, Binmin, Binmax, Sigdat1, Sigerr1 = gh.readcol5(ngp.files.Sigfiles[0])
# [Xscale0], [Munit/Xscale0^2]
# verified that indeed the stored files in the run directory are used
ngp.xipol = Radii * ngp.Xscale[0]       # [pc]
maxR = max(Radii)                     # [pc]
minR = min(Radii)                     # [pc]
Radii = np.hstack([minR/8, minR/4, minR/2, Radii, 2*maxR, 4*maxR, 8*maxR])
ngp.xepol = Radii * ngp.Xscale[0]       # [pc]

pc = pcload_single_entries(basedir, ngp)

with open(basedir+'pc', 'wb') as fn:
    pickle.dump(pc, fn)
