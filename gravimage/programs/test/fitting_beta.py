#!/usr/bin/env ipython3

import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import pdb
import gl_params
gp = gl_params.Params()
import gl_file as gf
gf.get_binned_data(gp)

import import_path as ip
ip.insert_sys_path('/home/psteger/sci/darcoda/gravlite/programs/sphere/')
import gl_analytic as ga
import gl_physics as phys

def analytic_beta(x):
    # gaia 1 pop
    # get betastar!
    ba = ga.beta_gaia(x, gp)[0]
    ba = phys.beta2betastar(ba)
    return ba
    #     r_a1 = 1000
    # hern return x**2/(x**2+r_a1**2)


x = gp.xipol
#xturn = gp.Xscale[0] #max(gp.xipol)/2
y = analytic_beta(x)

def modelbetaj(r0, vec0, vec1, vec2, vec3):
    betastartmp = (vec0-vec1)*np.exp(-(r0/vec2)**vec3)+vec1
    return betastartmp


def modelbeta_sigmoid(r0, vec0, vec1, vec2, vec3, vec4):
    # sigmoid with offsets
    s=np.log(r0/vec3)
    betastartmp = vec0/(1+np.exp(vec1*s+vec2))+vec3*np.ones(len(r0))
    return betastartmp


def modelbeta4(r0, vec0, vec1, vec2, vec3):
    s0=np.log(r0/vec3)
    a0 = vec0
    a1 = vec1
    alpha = vec2
    betastartmp = (a0-a1)/(1+np.exp(alpha*s0))+a1
    return betastartmp

fig = plt.figure(facecolor='white')
left, width = 0.25, 0.7
rect1 = [left, 0.4, width, 0.55]
rect2 = [left, 0.2, width, 0.2]
ax1 = fig.add_axes(rect1)  #left, bottom, width, height
ax2 = fig.add_axes(rect2, sharex=ax1)

ax1.plot(x, y, 'k.-', color='black', lw=5, label='data')
popt3, pcov3 = curve_fit(modelbeta4, x, y)
#poptj, pcovj = curve_fit(modelbetaj, x, y, p0=[0,1,100])

x = gp.xfine
y = analytic_beta(x)
ax1.plot(x, y, 'b', lw=2, label='analytic')
#ax1.axvline(xturn)
ax1.set_ylim([-0.2, 1.2])
ax1.set_xscale('log')
#ax1.plot(x, modelbeta4(x, *popt3), 'g.-', alpha=0.8, label='new sigmoid')
ax1.plot(x, phys.betastar(x, popt3, gp), 'r--', lw=2, label='fit on data')
#ax1.plot(x, modelbetaj(x, *poptj), alpha=0.8, label='new 4 param model')

ax1.set_yticks(np.linspace(0.0, 1.0, 6,endpoint=True))
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel('$\\beta*$')
legend = ax1.legend(loc='upper left', shadow=False, borderpad=0.2, labelspacing=0.1, handletextpad=0.1, borderaxespad=0.3)
# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame  = legend.get_frame()
frame.set_facecolor('1.0')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize(8)
for label in legend.get_lines():
    label.set_linewidth(2)  # the legend line width

ax2.set_xscale('log')
ax2.plot(x, y-modelbeta4(x, *popt3), 'r.-', alpha=0.8)
#ax2.plot(x, y-modelbetaj(x, *poptj), 'g.-', alpha=0.8)

#ax1.set_xticks(np.linspace(0, 2000, 3, endpoint=True))
#ax2.set_yticks(np.linspace(-0.05, 0.05, 3,endpoint=True))
#ax2.set_ylabel('$\\Delta\\beta$')
ax2.set_ylim([-0.001, 0.001])
ax2.set_xlabel('$r\\,[{\\rm pc}]$')
#ax.draw()
plt.savefig('beta_fit_all.pdf')

pdb.set_trace()
