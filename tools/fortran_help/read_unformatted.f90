#! /usr/bin/python
" volume =getvol(filename,nx,ny,nz) --reads 3D-array from single precision Fortran binary file"
import sys
from numpy import *
from struct import unpack

def readint(r,knt):
    bs = unpack('i',r[knt+4:knt+8])
    return bs[0],knt+12

def readfloat(r,knt):
    bs = unpack('d',r[knt+4:knt+12])
    return bs[0],knt+12

def getvol(f,nx):
    uraw = open(f,'rb')
    r = uraw.read()
#    print r[1:10]
    knt=0
    vol=zeros(nx)

    ncpu,knt=readint(r,knt)
    nvar,knt=readint(r,knt)
    ndim,knt=readint(r,knt)
    nlevelmax,knt=readint(r,knt)
    nboundary,knt=readint(r,knt)
    gamma,knt=readfloat(r,knt)
    
    print "ncpu:     ",ncpu
    print "nvar:     ",nvar
    print "ndim:     ",ndim
    print "nlevelmax:",nlevelmax
    print "nboundary:",nboundary
    print "gamma:    ",gamma
    
    for k in range(nx):
        bs = unpack('i',r[knt:knt+4])
        vol[k]=bs[0]
        knt=knt+4
    return vol

file = sys.argv[1]
v=getvol(file,10)
print v
