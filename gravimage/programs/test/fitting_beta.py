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
import gi_file as gf
gf.get_binned_data(gp)

import import_path as ip
ip.insert_sys_path('/home/psteger/sci/darcoda/gravimage/programs/sphere/')

import gravimage
gravimage.prepare_data(gp)

import gi_analytic as ga
import gi_physics as phys
import gi_class_cube as gcc

gp.debug = True

def traf(x):
    return np.arctan(x)/np.pi+0.5
# \fn traf(x)
# transform [-10,10] interval into [0,1]

def analytic_betastar(x):
    ba = ga.beta_gaia(x, gp)[0]
    ba = phys.beta2betastar(ba)
    return ba

def modelbeta4(r0, vec0, vec1, vec2, vec3):
    vec = np.array([vec0, vec1, vec2, vec3])
    vec = traf(vec)
    vec = gcc.map_betastar_sigmoid(vec, gp)
    ba = phys.betastar(r0, vec, gp)
    return ba

x = gp.xepol
y = analytic_betastar(x)
popt3, pcov3 = curve_fit(modelbeta4, x, y, p0=npr.rand(4))
beta01opt = traf(popt3)
betaopt = gcc.map_betastar_sigmoid(beta01opt, gp)
print('beta01opt = ', beta01opt)
print('betaopt = ',betaopt)

fig = plt.figure(facecolor='white')
left, width = 0.25, 0.7
rect1 = [left, 0.4, width, 0.55]
rect2 = [left, 0.2, width, 0.2]
ax1 = fig.add_axes(rect1)  #left, bottom, width, height
ax2 = fig.add_axes(rect2, sharex=ax1)
ax1.plot(x, y, 'b', lw=2, label='analytic')
ax1.plot(x, modelbeta4(x, *popt3), 'r--', lw=2, label='fit')
ax1.plot(x, phys.betastar(x, betaopt, gp), 'g--', lw=1, label='phys.betastar')
ax1.set_ylim([-0.2, 1.2])
ax1.set_xscale('log')
ax1.set_yticks(np.linspace(0.0, 1.0, 6,endpoint=True))
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel('$\\beta*$')
legend = ax1.legend(loc='upper left', shadow=False, borderpad=0.2, labelspacing=0.1, handletextpad=0.1, borderaxespad=0.3)
frame  = legend.get_frame()
frame.set_facecolor('1.0')
for label in legend.get_texts():
    label.set_fontsize(8)
for label in legend.get_lines():
    label.set_linewidth(2)  # the legend line width

ax2.set_xscale('log')
ax2.plot(x, y-modelbeta4(x, *popt3), 'r.-', alpha=0.8)
ax2.set_xlabel('$r\\,[{\\rm pc}]$')
plt.savefig('beta_fit_all.pdf')

beta01opt =  np.array([ 0.50066057,1., 0.51280261, 0.78100506])
betaopt =  np.array([ 1.30792522e-03, 9.89999992e-01, 2.05121046e+00, 5.85522159e+00])
