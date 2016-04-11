#!/usr/bin/env ipython3

## @file
# import file with full path specification, only once, setting back sys.path
# you need to specify PYTHONPATH to include the gravimage/programs

import os
import sys
import pdb
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


def set_geometry(geom):
    #Find darcoda path
    relative_path=''
    while os.path.abspath(relative_path).split('/')[-1] != 'darcoda':
        relative_path += '../'
        if os.path.abspath(relative_path).split('/')[-1] == 'home':
            raise Exception('Cannot find darcoda folder in function set_geometry')
    darcoda_path = os.path.abspath(relative_path)
    programs_path = darcoda_path + '/gravimage/programs/'

    #print('programs_path = ', programs_path)
    insert_sys_path(programs_path + 'datareduction/')
    insert_sys_path(programs_path + geom)
## \fn set_geometry(geom, machine)
# get right directory for geometry-dependent imports
# @param geom string of investigation geometry: disc or sphere
# @param machine string for working machine. local or darkside
