#!/usr/bin/python
# zooms in on one snapshot, with factor 2 between the snapshots

import os
import sys
import numpy
import math
import pylab
from pylab import *
import matplotlib
matplotlib.use('Agg') # use png output
from matplotlib import rc
rc('text',usetex=True)# use latex in figures
import Image
import fortranfile    # for fortran data import
import initialize as my

calc = True
show = True
run  = False

# check syntax
i=len(sys.argv)
if i!=4:
    print "Usage: zoom_steps.py snap_start snap_stop nzoom"
    exit(1)

nstart = int(sys.argv[1])
nstop  = int(sys.argv[2])
nzoom  = int(sys.argv[3])

x,y,z,r= my.get_xyzr(nstop)
r = r*2**nzoom

for i in range(nzoom):
    nc = nstart+i*(nstop-nstart)/nzoom
    folder = "output_"+str(nc).zfill(5)
    xc,yc,zc,rc=my.get_xyzr(nc)

    r = r/2
    xmi= str(xc-r); xma = str(xc+r)
    ymi= str(yc-r); yma = str(yc+r)
    zmi= str(zc-r/(i+1)); zma = str(zc+r/(i+1))

    bndry="-xmi "+xmi+" -xma "+xma
    bndry=bndry+" -ymi "+ymi+" -yma "+yma
    bndry=bndry+" -zmi "+zmi+" -zma "+zma
    
    lma = 9+i
    stri = str(i)
    cmd1 = "amr2map -typ 1 -lma "+str(lma)+" -inp "+folder+" -out mov/gas/gas_zoomstep_"+stri
    cmd1 = cmd1 +".dat -dir z "+bndry
    cmd2 = "map2img.py -l --colormap=hot mov/gas/gas_zoomstep_"+stri+".dat "
    cmd2 = cmd2+"-o mov/gas/gas_zoomstep_"+stri+".png"
    my.run(cmd1,cmd2,calc,show,run)


