#!/usr/bin/env python3

##
# @file
# select a completed or still running MCMC run for plotting

# (c) 2013 ETHZ, psteger@phys.ethz.ch

import os, sys, time, glob
import gl_params as gp
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


def list_files():
    search_dir = gp.files.dir
    # print('system: ', search_dir)
    
    from stat import S_ISREG, ST_CTIME, ST_MODE
    dirs = list(filter(os.path.isdir, glob.glob(search_dir + "201*")))
    
    # move *all* files *.dat to dat in all dirs of this walkercase
    if gp.investigate == 'walker':
        for x in dirs:
            files = list(filter(os.path.isfile, glob.glob(x + "/*")))
            for y in files:
                ending = y.split('.')[-1]
                # renaming only if not prof in name!
                if y == ending or ending == 'pdf' or ending == 'conf':
                    # only once! if no . there, get garbage
                    continue
                os.rename(y, x+"/"+ending) 
 
    
    from datetime import datetime
    fdl = [(x,\
            datetime.strptime(x[x.find('201'):x.find('201')+14], '%Y%m%d%H%M'),\
            bufcount(x+'/ev.dat')) for x in dirs]
    fdl.sort(key=lambda x: x[1])

    for i in range(len(fdl)):
        print("%2d"%(i+1),': ',\
              datetime.strftime(fdl[i][1],'%Y-%m-%d %H:%M'),\
              ' with ', "%7d"%fdl[i][2], ' iterations')
    import numpy as np
    return np.transpose(np.array(fdl))[:][0]
## \fn list_files
# return all working or completed MCMC run directories


def get_run(default):
    # interactive input
    invalid=True
    while(invalid):
        try:
            user_input = input('Input: (default '+str(default)+"): ")
            if not user_input:
                user_input = str(default)
            selection = int(user_input)
        except ValueError:
            print("error in input")
        if 0<selection and selection<=default:
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
            user_input = input('population: 1 (default), 2 (if available)')
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
            user_input = input('profile: rho (default), nr, M, beta, betastar, nu, sig: ')
            if not user_input:
                user_input = default
            sel = user_input
        except ValueError:
            print("error in input")
        if sel=='rho' or sel=='M' or sel=='beta' or sel=='betastar' or sel=='nu' or sel=='sig' or sel=='nr':
            invalid = False
    return sel
## \fn get_prof
# interactive input for prof: rho, M, beta, betastar, nu, sig
# @return string


def run():
    action = 'k'
    while(action=='k'):
        fdl = list_files()
        sel = get_run(len(fdl))
        action = get_action()
        if action == 'k':
            import shutil
            shutil.rmtree(fdl[sel])
        
    prof = get_prof()
    if prof == 'rho' or prof=='nr':
        return fdl[sel]+'/', prof, 0
    pop  = get_pop()
    return fdl[sel]+'/', prof, pop
## \fn run
# display possible runs of the current investigation method, select one
# @return basename, prof (string)

if __name__ == '__main__':
    run()
