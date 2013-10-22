#!/usr/bin/env ipython-python3.2
import os, sys, time, glob
import gl_params as gp
import pdb

def bufcount(filename):
    '''determine no. lines optimally'''
    f = open(filename)
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read # loop optimization
    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)
    return lines


def list_files():
    search_dir = gp.files.dir
    # print 'system: ', search_dir
    
    from stat import S_ISREG, ST_CTIME, ST_MODE
    dirs = list(filter(os.path.isdir, glob.glob(search_dir + "20*")))
    
    # move *all* files *.dat to dat in all dirs of this walkercase
    for x in dirs:
        files = list(filter(os.path.isfile, glob.glob(x + "/*")))
        for y in files:
            ending = y.split('.')[-1]
            # renaming only if not prof in name!
            if y == ending or ending == 'pdf' or ending == 'conf': # only once! if no . there, get garbage
                continue
            os.rename(y, x+"/"+ending) 

    
    from datetime import datetime
    fdl = [(x,\
            datetime.strptime(x[x.find('201'):x.find('201')+14], '%Y%m%d%H%M%S'),\
            bufcount(x+'/profM')) for x in dirs]
    fdl.sort(key=lambda x: x[1])

    for i in range(len(fdl)):
        print("%2d"%(i+1),': ',\
              datetime.strftime(fdl[i][1],'%Y-%m-%d %H:%M:%S'),\
              ' with ', "%7d"%fdl[i][2], ' iterations')
    import numpy as np
    return np.transpose(np.array(fdl))[:][0]


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


def get_prof():
    # interactive input for prof: dens, delta1, delta2, M
    default = 'dens'
    invalid=True
    while(invalid):
        try:
            user_input = input('profile: dens (default), M, delta1, delta2, nu1, sig1: ')
            if not user_input:
                user_input = default
            selection = user_input
        except ValueError:
            print("error in input")
        if selection == 'dens' or selection == 'M' or selection == 'delta1' or selection == 'delta2' or selection == 'nu1' or selection == 'sig1':
            invalid = False
    return selection


def run():
    '''display possible runs of the current investigation method,
    select one, plot'''
    action = 'k'
    while(action=='k'):
        fdl = list_files()
        sel = get_run(len(fdl))
        action = get_action()
        if action == 'k':
            import shutil
            shutil.rmtree(fdl[sel])
        
    prof = get_prof()
    return fdl[sel]+'/', prof


if __name__ == '__main__':
    run()
