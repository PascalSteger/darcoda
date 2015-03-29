#!/usr/bin/env ipython3

import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import pdb
import gi_params
gp = gi_params.Params()

import import_path as ip
ip.insert_sys_path('/home/psteger/sci/darcoda/gravimage/programs/sphere/')
import gi_analytic as ga

def analytic_beta(x, gp):
    #r_a1 = 1000
    # gaia 1 pop
    return ga.beta_gaia(x, gp)[0]
    # Hern: return x**2/(x**2+r_a1**2)


def model_beta(r0, vec):
    betatmp = np.zeros(len(r0))
    for k in range(len(vec)):
        betatmp += vec[k] * (r0/max(r0))**k

    # clipping beta* to the range [-1,1]
    # thus not allowing any unphysical beta,
    # but still allowing parameters to go the the max. value
    # for i in range(len(r0)):
    #     if betatmp[i] > 1.:
    #         betatmp[i] = 1.
    #     if betatmp[i] <= -1.:
    #         betatmp[i] = -0.99999999999 # excluding -inf values in beta

    return 2.*betatmp/(1.+betatmp)
## \fn model_beta(r0, vec)
# map [0,1] to [-1,1] with a polynomial
# @param r0 radii [pc]
# @param vec normalized ai, s.t. abs(sum(ai)) = 1

if __name__=="__main__":
    fig = plt.figure()
    ax  = fig.add_subplot(1,1,1)

    x = np.linspace(0, 3e3, 20)
    # y = func(x, 2.5, 1.3, 0.5)
    y = analytic_beta(x, gp)

    ax.plot(x, y, 'b--')
    #yn = y*(1 + 0.1*np.random.normal(size=len(x)))
    #ax.plot(x, yn, 'b')
    #ax.plot(x, model_beta(x, 0.5), 'r')

    pdb.set_trace()
    popt, pcov = curve_fit(model_beta, x, y)

    ax.plot(x, model_beta(x, *popt), 'g')
    #ax.draw()
    plt.savefig('beta_fit.pdf')

    print(popt)
    print(pcov)
beta01opt =  [ 0.50066057  0.99999999  0.51280263  0.78100505]
betaopt =  [  1.30792775e-03   9.89999983e-01   2.05121052e+00   5.85522158e+00]
