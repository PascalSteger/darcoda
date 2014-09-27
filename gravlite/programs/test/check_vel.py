#/usr/bin/env ipython3

## \file
# check the distribution of LOS velocities

import numpy as np
from pylab import *

vlos=np.loadtxt('/home/psteger/sci/darcoda/gravlite/DTobs/0/table_all.bin', skiprows=39, usecols=(11,))

hist(vlos, 150)

# TODO fit Maxwellian

xlabel(r'$v_{\text{LOS}}\quad\text{km/s}$')
ylabel(r'$f\quad[1]$')
ipdb.set_trace()
