#!/usr/bin/python
'''Plot fortran images.'''

import numpy
import array
import struct
import pylab
import math
import sys
import os
import toolKit2 as TK

# file-dependent variables.
istart=26
istop =27

RM=4
mp=1.67e-24
kb=1.3801e-16

for i in xrange(istart,istop):
	number=repr(i).zfill(5)
	filename="MAPS/g_dirx"+number+".dat"         
	command = "./fortranimage.py "+filename+" 1"
	print command 
	os.system(command)
        

