#!/usr/bin/python
import numpy as np

# ============================ BINSMOOTH ============================
# This routine takes an array(r) and bins it in r bins of size bin. If
# there is no data in a particular bin, it assigns a value of
# array(r)=nanflag.
# WARNING!! THIS ROUTINE REQUIRES SORTED ASCENDING r ARRAYS.

#PRO binsmooth,r,array,low,high,bin,rout,arrayout,nanflag,count_bin
# rout, arrayout, count_bin = 
def binsmoo(r, array, low, high, nbin, nanflag):

    bin = (high-low)/(1.*nbin) # binsize, TODO: check off-by-one errors
    Pnts = int(np.round((1.*high-low)/(1.*bin)))
    rout = np.arange(Pnts)*bin + low
    arrayout = np.zeros(Pnts)
    arrayout1 = np.zeros(Pnts)
    arrayout2 = np.zeros(Pnts)
    count_bin = np.zeros(Pnts)
    j=0
    siz = len(r)
    
    for i in range(Pnts):
        count=0
        while (rout[i] > r[j]):
            arrayout1[i]=arrayout1[i]+array[j]
            arrayout2[i]=arrayout2[i]+array[j]**2
            if (j < siz-1):
                j=j+1
                if (array[j] != 0):
                    count=count+1
            else:
                break

        if (count > 0):
            arrayout[i]=arrayout2[i]/count-(arrayout1[i]/count)**2
        else:
            arrayout[i]=nanflag
        count_bin[i] = count
    return rout, arrayout, count_bin
