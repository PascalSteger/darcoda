#!/usr/bin/env ipython3

##
# @file
# select a completed or still running MCMC run for plotting

# (c) 2015 ETHZ, Pascal Steger, pascal@steger.aero

import os
import glob
import shutil
import pdb
import re

import numpy as np
import gi_base as gb

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
# not used anymore
# @param filename filename

def removeDir(path):
    if os.path.isdir(path):
        shutil.rmtree(path)
## \fn removeDir(path)
# delete directory recursively, if it exists
# @param path path to directory

def remove_empty_folders(fdl):
    for x in fdl:
        try:
            with open(x[0]+"ev.dat"):
                # delete any folder that has empty ev.dat
                if bufcount(x[0]+'ev.dat') <= 0:
                    #removeDir(x[0])
                    print('empty '+x[0]+'ev.dat, remove dir?')
                continue
        except IOError:
            #removeDir(x[0])
            print(x[0]+'ev.dat does not exist, remove empty directory '+x[0]+'?')
    return
## \fn remove_empty_folders(fdl)
# remove all folders in list fdl without any ev.dat file in it = where
# no iterations of the MultiNest algorithm were stored
# @param fdl [(name, datestamp)] of all dirs in current mode

def list_files(basedir):
    dirs = list(filter(os.path.isdir, glob.glob(basedir + "201*")))
    from datetime import datetime
    fdl = [(x, datetime.strptime(x[x.find('201'):x.find('201')+14],\
                                 '%Y%m%d%H%M')) for x in dirs]
    #remove_empty_folders(fdl)
    dirs = list(filter(os.path.isdir, glob.glob(basedir + "201*")))
    fdl = []
    for x in dirs:
        timestamp = datetime.strptime(x[x.find('201'):x.find('201')+14], '%Y%m%d%H%M')
        try:
            co = bufcount(x+'ev.dat')
        except:
            print('file not found')
            co = 0
        fdl.append((x, timestamp, co ))
    fdl.sort(key=lambda x: x[1])
    for i in range(len(fdl)):
        try:
            fil = open(fdl[i][0]+'programs/gi_params.py','r')
        except:
            print('missing '+fdl[i][0]+'programs/gi_params.py')
            continue
        pops = 0                # default: 0 populations, error
        nipol = 0
        nbeta = 0
        for line in fil:
            lsp = line.split()
            if len(lsp) < 1:
                continue
            if lsp[0] == 'self.pops':
                pops = int(lsp[2])
            elif lsp[0] == 'self.nipol':
                nipol = int(lsp[2])
            elif lsp[0] == 'self.nbeta':
                nbeta = int(lsp[2])
            elif lsp[0] == 'self.monotonic':
                mono = lsp[2]
        print("%2d"%(i+1),': ',\
              datetime.strftime(fdl[i][1],'%Y-%m-%d %H:%M'),\
              "%5d"%fdl[i][2], 'its,', \
              pops, 'pops,',\
              nipol, 'bins,', 'nbeta=', nbeta, 'mono=', mono)
    out = np.transpose(np.array(fdl))[:][0]
    return out
## \fn list_files(basedir)
# return all working or completed MCMC run directories
# NOT USED ANYMORE
# @param basedir string of directory

def list_files_readout(basedir, investigate, case):
    bp = gb.get_basepath()
    fil = open(bp+"/run_info", "r")
    fdl=[]
    for line in fil:
        if re.search("DT"+str(investigate)+"/"+str(case)+"/", line):
            line2 = re.sub(r'\n', '', line)
            print(line2)
            #exclude "File not found" errors
            if not re.search("File ", line2):
                fdl.append(line2)
    fil.close()
    return fdl
## \fn list_files_readout(basedir, investigate, case)
# return list of MCMC run directories
# @param basedir string of base directory, e.g. DTgaia/3/
# @param investigate str investigation, gaia, walk, ..
# @param case int for case to look at

def get_investigate():
    default = 'gaia'
    invalid = True
    while(invalid):
        try:
            user_input = input('Investigate: (default: '+str(default)+\
                               ", walk, gaia, triax, hern, obs, discsim, discmock, simplenu): ")
            if not user_input:
                user_input = str(default)
            sel = user_input
        except ValueError:
            print("error in input")
        invalid = True
        if sel == 'walk' or sel == 'gaia' or sel == 'triax' or sel == 'hern'\
          or sel == 'obs' or sel == 'discsim' or sel == 'discmock' or sel == 'simplenu':
            invalid = False
    return sel
## \fn get_investigate(default)
# ask user for choice on investigation
# (walk, gaia, triax, fornax, hern, discsim, discmock)
# @return string

def choose_obs(sel):
    if sel == '1':
        return 'for'
    if sel == '2':
        return 'car'
    if sel == '3':
        return 'scl'
    if sel == '4':
        return 'sex'
    if sel == '5':
        return 'dra'
    return 'for'
## \fn choose_obs(sel)
# mapping of selection (int) to directory name (string of length 3)
# @param sel selection (int)

def get_case(investigate):
    default = 1
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
                       (sel <= 10 and investigate == 'gaia') or \
                       (sel <= 5 and investigate == 'obs') or \
                       (sel <= 8 and investigate == 'triax') or \
                       (sel == 0 and investigate == 'discmock') or\
                       (sel == 0 and investigate == 'simplenu')):
            invalid = False
        # assign string if working with observed dwarfs
        #if investigate == 'obs':
        #    sel = choose_obs(sel)
    return sel
## \fn get_case()
# ask user for choice on case number (see gi_params and gi_class_files
# for definitions)
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
# @return selection integer


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
        if sel=='1' or sel=='2':
            invalid = False
    return int(sel)
## \fn get_pop()
# get population integer
# @return int of population (0: all, 1: first, ...)


def get_prof():
    default = 'all'
    invalid=True
    while(invalid):
        try:
            user_input = input('profile: all, rho, nr, beta, betastar, Sig, sig, chi2: ')
            if not user_input:
                user_input = default
            sel = user_input
        except ValueError:
            print("error in input")
        if sel=='rho' or sel=='M' or sel=='beta' or sel=='betastar' or \
          sel=='Sig' or sel=='sig' or sel=='nr' or sel=='all' or sel=='chi2':
            invalid = False
    return sel
## \fn get_prof()
# interactive input for prof: rho, M, beta, betastar, nu, sig
# @return string


def get_pops(basename):
    pops = 1
    with open(basename+'programs/gi_params.py', 'r') as f:
        for line in f:
            lss = line.split()
            if len(lss) == 0:
                continue
            if lss[0]=='self.pops':
                pops = int(lss[2])
                # import re
                # pops = int(re.split('[=\n]', line)[-2])
    return pops
## \fn get_pops(basename)
# return number of populations from gi_params stored in output directory
# @param basename string of output dir
# @return integer number of populations


def get_nipol(basename):
    nipol = 0
    with open(basename+'programs/gi_params.py', 'r') as f:
        for line in f:
            if 'nipol      = ' in line:
                import re
                nipol = int(re.split('[=\n]', line)[-2])
    return nipol
## \fn get_nipol(basename)
# return number bins from gi_params stored in output directory
# @param basename string of output dir
# @return integer number of bins


def get_nbeta(basename):
    nbeta = 0
    with open(basename+'programs/gi_params.py', 'r') as f:
        for line in f:
            if 'nbeta      = ' in line:
                import re
                nbeta = int(re.split('[=\n]', line)[-2])
    return nbeta
## \fn get_nbeta(basename)
# return number of beta parameters  from gi_params stored in output directory
# @param basename string of output dir
# @return integer number of beta* parameters


def run(investigate="", case=-1, latest=False):
    if investigate == "":
        investigate = get_investigate()
    if case == -1:
        case = get_case(investigate)
    basepath = gb.get_basepath()
    basedir = os.path.abspath(basepath+'DT'+investigate+'/'+str(case)+'/')

    print(' - searching directory ', basedir, ' for output files')
    if latest:
        fdl = list_files_readout(basedir, investigate, case)
        sel = -1
    else:
        action = 'k'
        while(action == 'k'):
            fdl = list_files_readout(basedir, investigate, case)
            sel = get_run(len(fdl))
            action = get_action()
            if action == 'k':
                import shutil
                basename = re.split('\t', fdl[sel])[1]
                shutil.rmtree(basepath+basename)

    line = fdl[sel] # full directory path, without '/'
    basename = re.split('\t', line)[1]
    info = re.split('/', basename)
    timestamp = info[2]
    basename = info[0]+'/'+info[1]
    return timestamp, basepath+basename+'/'
## \fn run(investigate, case, latest)
# display possible runs of the current investigation method, select one
# @param investigate string of investigation case, hern, gaia, walk, discmock
# @param case int for case
# @param latest boolean if looking at latest only
# @return basename, prof (string)

if __name__ == '__main__':
    # default: take overall gi_params for parameters
    run()
