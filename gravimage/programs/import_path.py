#!/usr/bin/env ipython3

## @file
# import file with full path specification, only once, setting back sys.path
# you need to specify PYTHONPATH to include the gravimage/programs

# (c) GPL v3 2015 Pascal Steger, pascal@steger.aero

import os
import sys
import pdb
from imp import reload

import gi_base as gb

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

def remove_third():
    sys.path[2] = ''
## \fn remove_third()
# remove third entry (where gi_params was stored for timestamped restart runs)

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
    print('Machine = ', machine)
    basepath = gb.get_basepath()+'programs/'
    insert_sys_path(basepath + 'reducedata/')
    insert_sys_path(basepath + geom)
## \fn set_geometry(geom, machine)
# get right directory for geometry-dependent imports
# @param geom string of investigation geometry: disc or sphere
# @param machine string for working machine. local or darkside
