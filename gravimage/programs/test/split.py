#!/usr/bin/env ipython3
import numpy.random as npr
import numpy as np
import ipdb
import matplotlib
#matplotlib.use('pdf')

import matplotlib.pyplot as plt




def myerrorbar(models, nbins, binmin, binmax):
    histo, yedges = np.histogram(models, nbins)
    locminy = yedges[0]
    locmaxy = yedges[-1]
    plt.imshow(np.array([histo.T,histo.T]).T,
               extent=[binmin, binmax, yedges.min(), yedges.max()],
               cmap=plt.cm.Blues)
    return locminy, locmaxy

if __name__=="__main__":
    fig = plt.figure()
    N = 3
    dummy, bins = np.histogram(npr.uniform(0,1,100),N)
    binmin = bins[:-1]
    binmax = bins[1:]
    ax0 = fig.add_subplot(int(str(N)+str(1)+str(1)))
    miny = 1.5
    maxy = 1.5
    for i in range(N):
        model = npr.normal(1.5,0.3, 100)
        fig.add_subplot(int(str(N)+str(i+1)+str(1)), sharey=ax0)
        numiny, numaxy = myerrorbar(model, 10, binmin[i], binmax[i])
        miny = min(numiny, miny)
        maxy = max(numaxy, maxy)
        ax0.set_ylim([miny, maxy])
        plt.draw()
        plt.show()
