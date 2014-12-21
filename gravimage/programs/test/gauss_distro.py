#!/usr/bin/env ipython3

import numpy as np
import matplotlib.pyplot as plt

if __name__=="__main__":
    fig = plt.figure()
    ax  = fig.add_subplot(1,1,1)
    sigma=1.
    x = np.linspace(-5., 5., 100)
    y = np.exp(-x**2/(2*sigma))/np.sqrt(2*sigma)

    ax.set_xscale('linear')
    ax.set_yscale('log')

    sigma = 1
    ax.plot(x, y, lw=2)
    ax.axvline(x = -2*sigma, alpha=0.5)
    ax.axvline(x = -sigma,   alpha=0.5)
    ax.axvline(x =  sigma,   alpha=0.5)
    ax.axvline(x =  2*sigma, alpha=0.5)
