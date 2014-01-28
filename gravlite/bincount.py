#!/usr/bin/env python3

## @file
# bincount takes an array, r, and counts the number of elements 
# in r bins of size bin. 

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import numpy as np

def bincou(r, low, high, nbin):
    bin = (high-low)/(1.*nbin) # TODO: check off by one error
    Pnts = int(np.round((1.*high-low)/(1.*bin)))
    rout = np.arange(Pnts)*bin + low
    arrayout = np.zeros(Pnts)
    count_bin = np.zeros(Pnts)
    error = np.zeros(Pnts)
    j = 0
    siz = len(r)

    for i in range(Pnts):
        while (rout[i] > r[j]):
            arrayout[i] = arrayout[i] + 1.
            if (j < siz-1):
                j = j + 1
            else:
                break

        count_bin[i] = arrayout[i]
    
    return rout, arrayout, count_bin
## \fn bincou(r, low, high, nbin)
# take an array, r, and count the number of elements in r bins of size bin
# WARNING!! THIS ROUTINE REQUIRES SORTED ACSENDING r ARRAYS.
# @param r array of floats
# @param low lower limit
# @param high upper limit
# @param nbin number of bins
# @return rout, arrrayout, count_bin
