#!/usr/bin/env ipython3

import numpy as np
import numpy.random as npr
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import pdb
import gi_params
gp = gi_params.Params()
import gravimage
gravimage.prepare_data(gp)
import gi_file as gf
gf.get_binned_data(gp)

import import_path as ip
ip.insert_sys_path('/home/psteger/sci/darcoda/gravimage/programs/sphere/')
import gi_analytic as ga
import gi_physics as phys
import gi_class_cube as gcc

def traf(x):
    return np.arctan(x)/np.pi+0.5
# \fn traf(x)
# transform [-10,10] interval into [0,1]

def invtraf(y):
    return 2000*(y-0.5)

def analytic_rho(x):
    nn = ga.rho_gaia(x, gp)[1] # 0 for DM, 1 for stars
    return nn

def modelrho(r0, logrho, nr0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16, v17, v18, nrinf):
    # direct n(r) measures: traf(...)
    vec = np.array([logrho,nr0,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18,nrinf])
    vec = traf(vec)
    vec = gcc.map_nr_data(vec, 1, gp)
    rh = phys.rho(r0, vec, 0, gp)
    return np.log(rh)

gp.debug = False
x = gp.xfine
y = np.log(analytic_rho(x))
popt3, pcov3 = curve_fit(modelrho, x, y, p0=npr.rand(gp.nrho), maxfev=10000)
nr01opt = traf(popt3[:])
nropt = gcc.map_nr_data(nr01opt, 1, gp)
print('nr01opt = ', nr01opt)
print('nropt = ',nropt)

fig = plt.figure(facecolor='white')
left, width = 0.25, 0.7
rect1 = [left, 0.4, width, 0.55]
rect2 = [left, 0.2, width, 0.2]
ax1 = fig.add_axes(rect1)  #left, bottom, width, height
ax2 = fig.add_axes(rect2, sharex=ax1)
ax1.plot(gp.xipol, np.log(gp.dat.nu[1]), 'k', lw=3, label='data')
ax1.plot(x, y, 'b', lw=2, label='analytic')
ax1.plot(x, modelrho(x, *popt3), 'r--', lw=2, label='fit on data')
ax1.plot(x, np.log(phys.rho(x, nropt, 1, gp)), 'g--', lw=1, label='phys.rho')
ax1.set_xscale('log')
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel('$\\rho$')
legend = ax1.legend(loc='lower left', shadow=False, borderpad=0.2, labelspacing=0.1, handletextpad=0.1, borderaxespad=0.3)
frame  = legend.get_frame()
frame.set_facecolor('1.0')
for label in legend.get_texts():
    label.set_fontsize(8)
for label in legend.get_lines():
    label.set_linewidth(2)  # the legend line width
ax2.set_xscale('log')
ax2.plot(x, y-modelrho(x, *popt3), 'r.-', alpha=0.8)
ax2.set_xlabel('$r\\,[{\\rm pc}]$')
plt.savefig('nu_fit_all.pdf')

pdb.set_trace()
nr01opt =  np.array([ 0.5167274,   0.73752469,  0.32992638,  0.34188365,  0.35801464,  0.37579106, 0.11849702,  0.33328699,  0.46511328,  0.53514588,  0.50219587,  0.58231034, 0.68821362,  0.82705712,  0.42085304,  0.1773261,   0.2354134,   0.46580655,  0.52848736,  0.22904507,  0.5694699 ])
nropt =  np.array([ 0.15649498,  6.65618254,  0.10293663,  0.1087109,   0.13849277,  0.24371261, 0.62633345,  1.05913181,  1.43774113,  1.82346043,  2.20091446,  2.60007997,  2.98745825,  3.423104,    3.80766658,  4.2089698,   4.62950843,  4.91166037,  4.97380638,  4.99718073,  5.2277589 ])

gp.dat.nrnu = [np.array([ 0.15476906,  0.85086798,  0.9342867 ,  0.88161169,  0.83254241, 0.85086798,  0.99930431,  1.22211638,  1.47184763,  1.78910057, 2.1987677 ,  2.51961046,  2.80345393,  3.10336133,  3.88504346,  4.52442727,  4.88817769,  5.07880404,  4.83455511,  6.32165657,  4.88817769]),
               np.array([ 0.15476906,  0.85086798,  0.9342867 ,  0.88161169,  0.83254241,  0.85086798,  0.99930431,  1.22211638,  1.47184763,  1.78910057, 2.1987677 ,  2.51961046,  2.80345393,  3.10336133,  3.88504346,  4.52442727,  4.88817769,  5.07880404,  4.83455511,  6.32165657,  4.88817769]),
               np.array([ 0.15476906,  0.85086798,  0.9342867 ,  0.88161169,  0.83254241, 0.85086798,  0.99930431,  1.22211638,  1.47184763,  1.78910057,  2.1987677 ,  2.51961046,  2.80345393,  3.10336133,  3.88504346,  4.52442727,  4.88817769,  5.07880404,  4.83455511,  6.32165657,    4.88817769]),
               np.array([ 0.15476906,  0.85086798,  0.9342867 ,  0.88161169,  0.83254241, 0.85086798,  0.99930431,  1.22211638,  1.47184763,  1.78910057,  2.1987677 ,  2.51961046,  2.80345393,  3.10336133,  3.88504346,  4.52442727,  4.88817769,  5.07880404,  4.83455511,  6.32165657,  4.88817769])]

gp.dat.nrnuerr = [np.array([  0.05158969,  12.22044422,   2.44408884,   2.44408884, 2.44408884,   2.44408884,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777, 0.48881777,   2.44408884,   2.44408884,   2.44408884,   2.44408884]),
                  np.array([  0.05158969,  12.22044422,   2.44408884,   2.44408884, 2.44408884,   2.44408884,   0.48881777,   0.48881777,    0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777, 0.48881777,   2.44408884,   2.44408884,   2.44408884,   2.44408884]),
                  np.array([  0.05158969,  12.22044422,   2.44408884,   2.44408884, 2.44408884,   2.44408884,   0.48881777,   0.48881777, 0.48881777,   0.48881777,   0.48881777,   0.48881777,  0.48881777,   0.48881777,   0.48881777,   0.48881777, 0.48881777,   2.44408884,   2.44408884,   2.44408884,   2.44408884]),
                  np.array([  0.05158969,  12.22044422,   2.44408884,   2.44408884, 2.44408884,   2.44408884,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777,   0.48881777, 0.48881777,   2.44408884,   2.44408884,   2.44408884,   2.44408884])]
