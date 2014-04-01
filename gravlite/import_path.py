#!/usr/bin/env python3

## @file
# import file with full path specification, only once, setting back sys.path

import os, sys, pdb
from imp import reload


def insert_sys_path(fullpath):
    sys.path.insert(0, fullpath)
    return sys.path


def import_path(fullpath):
    path, filename = os.path.split(fullpath)
    filename, ext = os.path.splitext(filename)
    sys.path.insert(0, path)
    module = __import__(filename)
    reload(module) # Might be out of date
    return module
## \fn import_path(fullpath)
# Import a file with full path specification.
# Allows one to import from anywhere, something __import__ does not do.
# @param fullpath path to filename
# @return module, use like module.some_fun()


def set_geometry(inv, machine):
    if machine == 'local':
        basepath = '/home/psteger/sci/gravlite/programs/'
    elif machine == 'darkside':
        basepath = '/home/ast/read/dark/gravlite/programs/'

    if inv == 'discmock' or inv == 'discsim':
        insert_sys_path(basepath + 'disc/')
    else:
        insert_sys_path(basepath + 'sphere/')
        from gl_priors import check_rho, check_beta
        from gl_project import rho_SUM_Mr
