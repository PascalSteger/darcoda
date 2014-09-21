#!/usr/bin/env ipython3

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pdb
import gl_params
gp = gl_params.Params()
import gl_file as gf
gf.get_data(gp)

import import_path as ip
ip.insert_sys_path('/home/psteger/sci/gravlite/programs/sphere/')
import gl_analytic as ga
import gl_physics as phys

def analytic_beta(x):
    # gaia 1 pop
    # use betastar!
    ba = ga.beta_gaia(x, gp)[0]
    ba = phys.beta2betastar(ba)
    return ba
    #     r_a1 = 1000
    # hern return x**2/(x**2+r_a1**2)


def modelbeta0(r0, arr0):
    betatmp = arr0 * (r0/max(r0))**0
    return betatmp
def modelbeta1(r0, arr0, arr1):
    betatmp = arr0 * (r0/max(r0))**0
    betatmp += arr1 * (r0/max(r0))**1
    return betatmp
def modelbeta2(r0, arr0, arr1, arr2):
    betatmp = arr0 * (r0/max(r0))**0
    betatmp += arr1 * (r0/max(r0))**1
    betatmp += arr2 * (r0/max(r0))**2
    return betatmp
def modelbeta3(r0, arr0, arr1, arr2, arr3):
    betatmp = arr0 * (r0/max(r0))**0
    betatmp += arr1 * (r0/max(r0))**1
    betatmp += arr2 * (r0/max(r0))**2
    betatmp += arr3 * (r0/max(r0))**3
    return betatmp
def modelbeta4(r0, arr0, arr1, arr2, arr3, arr4):
    betatmp = arr0 * (r0/max(r0))**0
    betatmp += arr1 * (r0/max(r0))**1
    betatmp += arr2 * (r0/max(r0))**2
    betatmp += arr3 * (r0/max(r0))**3
    betatmp += arr4 * (r0/max(r0))**4
    return betatmp
def modelbeta5(r0, arr0, arr1, arr2, arr3, arr4, arr5):
    betatmp = arr0 * (r0/max(r0))**0
    betatmp += arr1 * (r0/max(r0))**1
    betatmp += arr2 * (r0/max(r0))**2
    betatmp += arr3 * (r0/max(r0))**3
    betatmp += arr4 * (r0/max(r0))**4
    betatmp += arr5 * (r0/max(r0))**5
    return betatmp


fig = plt.figure(facecolor='white')
left, width = 0.25, 0.7
rect1 = [left, 0.4, width, 0.55]
rect2 = [left, 0.2, width, 0.2]
ax1 = fig.add_axes(rect1)  #left, bottom, width, height
ax2 = fig.add_axes(rect2, sharex=ax1)

x = np.hstack([0, gp.xepol])
# np.linspace(0, 3e3, 20)
y = analytic_beta(x)

ax1.plot(x, y, 'b--', color='black', lw=3, label='analytic beta*')
popt0, pcov0 = curve_fit(modelbeta0, x, y)
popt1, pcov1 = curve_fit(modelbeta1, x, y)
popt2, pcov2 = curve_fit(modelbeta2, x, y)
popt3, pcov3 = curve_fit(modelbeta3, x, y)
popt4, pcov4 = curve_fit(modelbeta4, x, y)
popt5, pcov5 = curve_fit(modelbeta5, x, y)

ax1.plot(x, modelbeta0(x, *popt0), alpha=0.8, label='$0th order$')
ax1.plot(x, modelbeta1(x, *popt1), alpha=0.8, label='$1st order$')
ax1.plot(x, modelbeta2(x, *popt2), alpha=0.8, label='$2nd order$')
ax1.plot(x, modelbeta3(x, *popt3), alpha=0.8, label='$3rd order$')
#
ax1.plot(x, modelbeta4(x, *popt4), alpha=0.8, label='$4th order$')
ax1.plot(x, modelbeta5(x, *popt5), alpha=0.8, label='$5th order$')
ax1.plot(x, phys.betastar(x, max(x), popt5, gp), 'r--', lw=2, label='phys.beta')

ax1.set_yticks(np.linspace(0.0, 1.0, 6,endpoint=True))
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel('$\\beta*$')
legend = ax1.legend(loc='lower right', shadow=False, borderpad=0.2, labelspacing=0.1, handletextpad=0.1, borderaxespad=0.3)
# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame  = legend.get_frame()
frame.set_facecolor('1.0')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize(8)
for label in legend.get_lines():
    label.set_linewidth(2)  # the legend line width

#ax2.plot(x, 0.0*x, 'b--', lw=2)
ax2.plot(x, y-modelbeta0(x, *popt0), alpha=0.8)
ax2.plot(x, y-modelbeta1(x, *popt1), alpha=0.8)
ax2.plot(x, y-modelbeta2(x, *popt2), alpha=0.8)
ax2.plot(x, y-modelbeta3(x, *popt3), alpha=0.8)
#ax2.plot(x, y-phys.betastar(x, popt3, gp), 'r--', lw=2, label='phys.beta')
ax2.plot(x, y-modelbeta4(x, *popt4), alpha=0.8)
ax2.plot(x, y-modelbeta5(x, *popt5), alpha=0.8)
ax1.set_xticks(np.linspace(0, 2000, 3, endpoint=True))
ax2.set_yticks(np.linspace(-0.05, 0.05, 3,endpoint=True))
#ax2.set_ylabel('$\\Delta\\beta$')
ax2.set_xlabel('$r\\,[{\\rm pc}]$')
#ax.draw()
plt.savefig('beta_fit_all.png')

from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('beta_fit_all.pdf')
pp.savefig(fig)
pp.close()

pdb.set_trace()

