#!/usr/bin/env ipython3

import numpy as np
import numpy.random as npr
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.ioff()

fig = plt.figure()
ax  = fig.add_subplot(111)

# first create a single histogram
mu, sigma = 200, 25
x = mu + sigma*npr.randn(10000)

bins, edges = np.histogram(x, range=[0.,220.], bins=30, density=True)

ax.step(edges[1:], bins, where='pre')

plt.show(block=True)
plt.savefig('myhist.png')
