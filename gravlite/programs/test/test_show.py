#!/bin/env ipython3

from pylab import *
ion()

import ipdb
import numpy as np

x = np.linspace(0, 2*np.pi, 100)
y = np.sin(x)


plot(x, y, label='sin')
ipdb.set_trace()

print('too late')
