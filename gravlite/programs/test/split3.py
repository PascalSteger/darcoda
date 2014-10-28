#!/usr/bin/env ipython
import pdb
import matplotlib.pyplot as plt
plt.ion()
import numpy as np

# Generate some data
Npoints = 690000
y = np.random.normal(1, 0.5, Npoints)
x = np.random.uniform(0,1000,Npoints)

Nbin = 12
countinbinx, xbins = np.histogram(x, Nbin)
ybins = np.linspace(min(y), max(y), num=100)

H, xedges, yedges = np.histogram2d(x, y, bins=[xbins, ybins])

# Plot the results
fig, ax = plt.subplots(figsize=(8, 8), dpi=80, facecolor='w', edgecolor='k')

extent = [min(xbins), max(xbins), min(ybins), max(ybins)]

im = ax.imshow(H.T, extent=extent, aspect='auto', cmap=plt.cm.Blues, interpolation='nearest')
fig.colorbar(im)
plt.draw()
plt.show()
pdb.set_trace()
