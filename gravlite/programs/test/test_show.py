#!/bin/env ipython3

from pylab import *
ion()

import pdb
import numpy as np

x = np.linspace(0, 2*np.pi, 100)
y = np.sin(x)


plot(x, y, label='sin')
pdb.set_trace()

print('too late')
