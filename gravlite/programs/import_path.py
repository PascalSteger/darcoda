#!/usr/bin/env ipython3

## @file
# import file with full path specification, only once, setting back sys.path

import os, sys, ipdb
from imp import reload


def insert_sys_path(fullpath):
    sys.path.insert(0, fullpath)
    return sys.path
## \fn insert_sys_path(fullpath)
# insert full path at the beginning of sys.path
# and thus make it the first directory to be scanned
# @param fullpath string to be inserted


def remove_first():
    sys.path = sys.path[1:]
## \fn remove_first()
# remove first entry in sys.path
# and thus reverse insert_sys_path


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


def set_geometry(geom, machine):
    if machine == 'pstgnt332':
        basepath = '/home/psteger/sci/darcoda/gravlite/programs/'
    elif machine == 'darkside':
        basepath = '/home/ast/read/user/psteger/software/darcoda/gravlite/programs/'
    elif machine == 'lisa_HS':
        basepath = '/home/hsilverw/LoDaM/darcoda/gravlite/programs/'
    elif machine == 'lisa_SS':
        basebath = '/home/sofia/blah/darcoda/gravlite/programs/'

    insert_sys_path(basepath + 'datareduction/')
    insert_sys_path(basepath + geom)
## \fn set_geometry(geom, machine)
# get right directory for geometry-dependent imports
# @param geom string of investigation geometry: disc or sphere
# @param machine string for working machine. local or darkside
