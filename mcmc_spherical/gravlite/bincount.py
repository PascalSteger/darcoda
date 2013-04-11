#!/usr/bin/python2.7
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

# ============================ BINCOUNT ============================
# This routine takes an array, r, and counts the number of elements 
# in r bins of size bin. 
# WARNING!! THIS ROUTINE REQUIRES SORTED ACSENDING r ARRAYS.

import numpy as np

#PRO bincount,r,low,high,bin,rout,arrayout,count_bin
# return rout, arrrayout, count_bin
def bincou(r, low, high, nbin):

    bin = (high-low)/(1.*nbin) # TODO: check off by one error
    Pnts = int(np.round((1.*high-low)/(1.*bin)))
    rout = np.arange(Pnts)*bin + low
    arrayout = np.zeros(Pnts)
    count_bin = np.zeros(Pnts)
    error = np.zeros(Pnts)
    j=0
    siz = len(r)

    for i in range(Pnts):
        while (rout[i] > r[j]):
            arrayout[i]=arrayout[i]+1.
            if (j < siz-1):
                j=j+1
            else:
                break

        count_bin[i] = arrayout[i]
    
    return rout, arrayout, count_bin
