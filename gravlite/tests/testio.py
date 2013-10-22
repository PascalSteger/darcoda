#!/usr/bin/python
# (c) 2013 Pascal S.P. Steger

import numpy
import os

os.system("rm sample.txt")
d=open('sample.txt','w')
print>>d,"sample array with missing values"
print>>d,"second line to skip"
print>>d,"1 2 0 1 "
print>>d,"2 2 2 3"
print>>d,"1 8   6"
print>>d,"3 2 1 2"
d.close()
os.system('echo "";cat sample.txt')

#numpy.genfromtxt("sample.txt", dtype='float', comments='#', delimiter=None, skiprows=0,\
                    # skip_header=0, skip_footer=0, converters=None,\
                    # missing='', missing_values=None, filling_values=None,\
                    # usecols=None, names=None, excludelist=None, deletechars=None,\
                    # replace_space='_', autostrip=False, case_sensitive=True,\
                    # defaultfmt='f%i', unpack=None, usemask=False, loose=True)

R1,R2,R3,R4=numpy.genfromtxt("sample.txt",skiprows=2,usecols=(0,1,2,3),unpack=True,delimiter=(2,2,2,2),filling_values="-1")
print R1
print R2
print R3
print R4
