#/usr/bin/env python3

## \file
# check the distribution of LOS velocities

import numpy as np
from pylab import *
ion()

vlos=np.loadtxt('/home/psteger/sci/darcoda/gravimage/DTobs/0/table_all.bin', skiprows=39, usecols=(11,))

hist(vlos, 150)

xlabel(r'$v_{\text{LOS}}\quad\text{km/s}$')
ylabel(r'$f\quad[1]$')
pdb.set_trace()
