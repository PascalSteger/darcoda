#!/usr/bin/env ipython3

##
# @file
# select a completed or still running MCMC run for plotting

# (c) 2013 ETHZ, psteger@phys.ethz.ch

import os, sys, time, glob, shutil
import numpy as np
import pdb

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
# determine no. lines optimally
# @param filename filename


def removeDir(path):
    if os.path.isdir(path):
        shutil.rmtree(path)
# \fn removeDir(path)
# delete directory recursively, if it exists
# @param path path to directory


def remove_empty_folders(fdl):
    for x in fdl:
        try:
            with open(x[0]+"/ev.dat"):
                # delete any folder that has empty ev.dat
                if bufcount(x[0]+'/ev.dat')<=0:
                    print('empty '+x[0]+'/ev.dat, delete containing directory')
                    removeDir(x[0])
                continue
        except IOError:
            print(x[0]+'/ev.dat does not exist, delete containing directory')
            removeDir(x[0])
    return
# \fn remove_empty_folders(fdl)
# remove all folders in list fdl without any ev.dat file in it = where no iterations of the MultiNest algorithm were stored
# @param fdl [(name, datestamp)] of all dirs in current mode


def list_files(basedir):
    from stat import S_ISREG, ST_CTIME, ST_MODE
    dirs = list(filter(os.path.isdir, glob.glob(basedir + "201*")))
    
    from datetime import datetime
    # delete all folders without or with empty ev.dat inside!
    fdl = [(x, datetime.strptime(x[x.find('201'):x.find('201')+14], '%Y%m%d%H%M')) for x in dirs]
    remove_empty_folders(fdl)
    dirs = list(filter(os.path.isdir, glob.glob(basedir + "201*")))
    fdl = [(x,\
            datetime.strptime(x[x.find('201'):x.find('201')+14], '%Y%m%d%H%M'),\
            bufcount(x+'/ev.dat')) for x in dirs]

    fdl.sort(key=lambda x: x[1])

    for i in range(len(fdl)):
        fil = open(fdl[i][0]+'/gl_params.py','r')
        pops = 0                          # default: 0 populations, error
        for line in fil:
            if ('self.pops' in line.split()): # new version only (old
                                              #would have to search
                                              #for 'pops' only)
                pops = int(line.split()[2])
        print("%2d"%(i+1),': ',\
              datetime.strftime(fdl[i][1],'%Y-%m-%d %H:%M'),\
              ' with ', "%7d"%fdl[i][2], ' iterations', \
              ' and ', pops, ' populations')
    return np.transpose(np.array(fdl))[:][0]
## \fn list_files(basedir)
# return all working or completed MCMC run directories
# @param basedir string of directory


def get_investigate():
    default = 'walk'
    invalid=True
    while(invalid):
        try:
            user_input = input('Investigate: (default: '+str(default)+", gaia, triax, hern, obs, discsim, discmock): ")
            if not user_input:
                user_input = str(default)
            sel = user_input
        except ValueError:
            print("error in input")
        invalid = True
        if sel == 'walk' or sel == 'gaia' or sel == 'triax' or sel == 'hern'\
          or sel == 'obs' or sel == 'discsim' or sel == 'discmock':
            invalid = False
    return sel
## \fn get_investigate(default)
# ask user for choice on investigation (walk, gaia, triax, fornax, hern, discsim, discmock)
# @return string


def choose_obs(sel):
    if sel == '0':
        return 'for'
    if sel == '1':
        return 'car'
    if sel == '2':
        return 'scl'
    if sel == '3':
        return 'sex'
    return 'for'
## \fn choose_obs(sel)
# mapping of selection (int) to directory name (string of length 3)
# @param sel selection (int)


def get_case(investigate):
    default = 4
    invalid = True
    while(invalid):
        try:
            user_input = input('Case: (default '+str(default)+"): ")
            if not user_input:
                user_input = str(default)
            sel = int(user_input)
        except ValueError:
            print("error in input")
        if 0<=sel and ((sel <= 5 and investigate == 'walk') or \
                      (sel <= 8 and investigate == 'gaia') or \
                       (sel <=4 and investigate == 'obs')):
            invalid = False
        # assign string if working with observed dwarfs
        if investigate == 'obs':
            sel = choose_obs(sel)
    return sel
## \fn get_case()
# ask user for choice on case number (see gl_params and gl_class_files for definitions)
# @return sel integer for case number


def get_run(default):
    # interactive input
    invalid = True
    while(invalid):
        try:
            user_input = input('Input: (default '+str(default)+"): ")
            if not user_input:
                user_input = str(default)
            selection = int(user_input)
        except ValueError:
            print("error in input")
        if 0 < selection and selection <= default:
            invalid = False
    return selection - 1
## \fn get_run(default)
# ask user for choice on run number
# @param default default value: last one
# @return selection-1 integer 


def get_action():
    # interactive input for action: p - print, k - kill and delete
    default = 'p'
    invalid=True
    while(invalid):
        try:
            user_input = input('Action: (p - print (default), k - kill): ')
            if not user_input:
                user_input = default
            selection = user_input
        except ValueError:
            print("error in input")
        if selection == 'p' or selection == 'k':
            invalid = False
    return selection
## \fn get_action
# ask user for choice on action to take: print or delete


def get_pop():
    default = '1'
    invalid = True
    while(invalid):
        try:
            user_input = input('population: 1 (default), 2 (if available): ')
            if not user_input:
                user_input = default
            sel = user_input
        except ValueError:
            print('error in input')
        if sel=='1' or (sel=='2' and gp.pops>1):
            invalid = False
    return int(sel)
## \fn get_pop
# @return int of population (0: all, 1: first, ...)


def get_prof():
    default = 'rho'
    invalid=True
    while(invalid):
        try:
            user_input = input('profile: rho (default), nr, M, beta, betastar, Sig, sig: ')
            if not user_input:
                user_input = default
            sel = user_input
        except ValueError:
            print("error in input")
        if sel=='rho' or sel=='M' or sel=='beta' or sel=='betastar' or sel=='Sig' or sel=='sig' or sel=='nr':
            invalid = False
    return sel
## \fn get_prof
# interactive input for prof: rho, M, beta, betastar, nu, sig
# @return string




def get_pops(basename):
    pops = 1
    with open(basename+'/gl_params.py', 'r') as f:
        for line in f:
            if 'pops      = ' in line:
                import re
                pops = int(re.split('[=\n]', line)[-2])
    return pops
## \fn get_pops(basename)
# return number of populations from gl_params stored in output directory
# @param basename string of output dir
# @return integer number of populations


def run():
    action = 'k'
    investigate = get_investigate()
    case = get_case(investigate)
    if os.path.exists('/home/ast/read'):
        machine = '/home/ast/read/dark/gravlite/'
    else:
        machine = '/home/psteger/sci/gravlite/'
    basedir = os.path.abspath(machine+'/DT'+investigate+'/'+str(case)+'/')+'/'
    
    # import import_path as ip
    # ip.import_path(basedir+'/gl_params.py')
    while(action == 'k'):
        fdl = list_files(basedir)
        sel = get_run(len(fdl))
        action = get_action()
        if action == 'k':
            import shutil
            shutil.rmtree(fdl[sel])
        
    prof = get_prof()
    basename = fdl[sel] # full directory path, without '/'
    timestamp = basename.split('/')[-1] # extract last bit
    if prof == 'rho' or prof=='nr':
        return timestamp, basename+'/', prof, 0
    pop  = get_pop()
    return timestamp, basename+'/', prof, pop
## \fn run()
# display possible runs of the current investigation method, select one
# @return basename, prof (string)

if __name__ == '__main__':
    # default: take overall gl_params for parameters
    run()
