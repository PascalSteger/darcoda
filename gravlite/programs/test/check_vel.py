#/usr/bin/env ipython3

## \file
# check the distribution of LOS velocities

import gl_plot as gpl
import numpy as np

vlos=np.loadtxt('/home/psteger/sci/darcoda/gravlite/DTobs/0/table_all.bin', skiprows=39, usecols=(11,))
fig, ax=gpl.start()

ax.hist(vlos, 150)

# TODO fit Maxwellian

ax.set_xlabel(r'$v_{\text{LOS}}\quad\text{km/s}$')
ax.set_ylabel(r'$f\quad[1]$')
plt.draw('gtk')
plt.show()
