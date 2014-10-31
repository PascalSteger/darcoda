#!/usr/bin/env ipython

import numpy as np
import numpy.random as npr

import matplotlib.pyplot as plt
fig = plt.figure()


xval = npr.normal(9, 3, 100)
yval = npr.normal(1, 0.5, 100)

nbins = 12
countx, xbins = np.histogram(xval, nbins)
for i in range(nbins):
    histo, yedges = np.histogram(models, nbins)


plt.imshow(np.array([histo.T,histo.T]).T,
           extent=[binmin, binmax, yedges.min(), yedges.max()],
           cmap=plt.cm.Blues)
