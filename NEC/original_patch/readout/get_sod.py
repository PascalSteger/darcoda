#!/usr/bin/python
import sys
import os
import toolKit2 as TK

istart=270
istop =271


for i in xrange(istart,istop):
        number=repr(i).zfill(5)
        filename="../run/output_"+number
	command='./sod -inp '+filename+' -del 200. -sel 8'
	print command
	os.system(command)
