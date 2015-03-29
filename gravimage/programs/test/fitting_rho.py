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
    nn = ga.rho_gaia(x, gp)[0] # 0 for DM, 1 for stars
    return nn

def modelrho(r0, rhohalf, nr00, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16, v17, v18, nrinf):
    # direct n(r) measures: traf(...)
    vec = np.array([nr00,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18, nrinf])
    #print('vec=',vec)
    #pdb.set_trace()
    vec = traf(vec)
    vec = gcc.map_nr(np.hstack([rhohalf, vec]), 'rho', 0, gp)
    #rhoparam = np.hstack([rhohalf,vec[0],vec,vec[-1]])
    rh = phys.rho(r0, vec, 0, gp)
    return np.log(rh)

gp.debug = False
x = gp.xfine
y = np.log(analytic_rho(x))
popt3, pcov3 = curve_fit(modelrho, x, y, p0=npr.rand(gp.nrho), maxfev=10000)
nr01opt = traf(popt3[1:])
nropt = gcc.map_nr(np.hstack([popt3[0], nr01opt]), 'rho', 0, gp)
print('nr01opt = ', nr01opt)
print('nropt = ', nropt)

fig = plt.figure(facecolor='white')
left, width = 0.25, 0.7
rect1 = [left, 0.4, width, 0.55]
rect2 = [left, 0.2, width, 0.2]
ax1 = fig.add_axes(rect1)  #left, bottom, width, height
ax2 = fig.add_axes(rect2, sharex=ax1)
ax1.plot(x, y, 'b', lw=2, label='analytic')
ax1.plot(x, modelrho(x, *popt3), 'r--', lw=2, label='fit on data')
ax1.plot(x, np.log(phys.rho(x, nropt, 0, gp)), 'g--', lw=1, label='phys.rho')
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
plt.savefig('rho_fit_all.pdf')


nr01opt =  np.array([ 0.32162355,  0.79255368,  0.50936768,  0.52601379,  0.54526758,  0.57555981,  0.5790085,   0.60252249,  0.60668556,  0.62252676,  0.63173716,  0.64555501,0.65777126,  0.67083571,  0.68506604,  0.69139872,  0.66304763,  0.61462276,  0.73605097,  0.67611218])
nropt =  np.array([ 0.18235592,  0.31022168,  0.,          0.03977087,  0.0641871,   0.13203218, 0.24268292,  0.32983781,  0.39875718,  0.46647496,  0.53362559,  0.60953939,  0.68948985,  0.79792504,  0.91214721,  1.08483156,  1.36074729,  1.88041945,  2.31792881,  2.62089113,  3.001     ])
