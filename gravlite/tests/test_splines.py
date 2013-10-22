#!/usr/bin/env python3

import numpy as np
import numpy.random as npr
import pdb
from scipy.interpolate import splrep, splev
import gl_plot as gpl


xnew = np.array([   0.   ,   40.156,   58.552,   74.426,   91.87 ,  110.166,
        130.44 ,  155.967,  192.224,  248.297,  255.215,  262.119,
        269.01 ,  275.887])

ynew = np.array([ 2426.404,  2294.581,  1828.642,  1144.447,   918.291,   681.606,
         321.61 ,   200.627,   116.21 ,    19.241,    17.438,    14.488,
          12.037,    10.   ])

errorinpercent = npr.uniform(0.05, 0.15, len(xnew))

gpl.plot(xnew,ynew,'.',lw=2)
gpl.fill_between(xnew,ynew*(1.-errorinpercent),ynew*(1.+errorinpercent),alpha=0.5,color='blue')

# gives warning
pdb.set_trace()
tcknu  = splrep(xnew, ynew, k=2, s=0.1) # interpolation in real space
gpl.plot(xnew, splev(xnew,splrep(xnew,ynew, k=2, s=0.1)))

# try with third order
tcknu  = splrep(xnew, ynew, k=3, s=0.1) # interpolation in real space

# better, uses weights

w = 1/(errorinpercent*ynew)
tcknu  = splrep(xnew, ynew, w=w, k=2, s=0.1) # interpolation in real space


# 