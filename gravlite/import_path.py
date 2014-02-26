#!/usr/bin/env python3

## @file
# import file with full path specification, only once, setting back sys.path

import os, sys
from imp import reload

def import_path(fullpath):
    path, filename = os.path.split(fullpath)
    filename, ext = os.path.splitext(filename)
    sys.path.append(path)
    module = __import__(filename)
    reload(module) # Might be out of date
    del sys.path[-1]
    return module
## \fn import_path(fullpath)
# Import a file with full path specification.
# Allows one to import from anywhere, something __import__ does not do.
# @param fullpath path to filename
# @return module, to be used like "TODO"
