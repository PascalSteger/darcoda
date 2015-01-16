#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import scoreatpercentile

mu, sigma = 100, 20
N = 100
x = mu + sigma*np.random.randn(N)

q1    = scoreatpercentile(x,25)
q3    = scoreatpercentile(x,75)
iqd   = q3 - q1

dx    = 2*iqd*N**(-1/3.)
nbins = (max(x)-min(x))/dx

hist, bins = np.histogram(x,bins=nbins)
width = 0.7 * (bins[1]-bins[0])
center = (bins[:-1]+bins[1:])/2
plt.bar(center, hist, align='center', width=width, facecolor='green')
#plt.hist(x,facecolor='red',bins=nbins*2,histtype='step',color='blue')
#plt.hist(x,facecolor='red',bins=nbins/2,histtype='step',color='red')
plt.show()
