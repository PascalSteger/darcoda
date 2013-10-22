#!/usr/bin/python
# from Numerical Recipes in C, chapter 14
# straight forward determination of mean, average deviation, standard deviation, variance, skewness and kurtosis

import numpy as np

def moment(data):

    n = 1.*len(data)
    if n <= 1:
        print "n must be at least 2 in moment"
        exit(1)


    s = 0. # First pass to get the mean
    for j in range(n):
        s += data[j]
    ave = s/n

    ep = 0.0; adev = 0.; var = 0.; skew = 0.; curt = 0.
    for j in range(n):
        s = data[j] - ave
        adev += fabs(s)                 # TODO: C++ fabs() in python?
        ep   += s
        p = s*s
        var  += p
        p *= s
        skew += p
        p *= s
        curt += p

    adev /= n;
    var=(var-ep*ep/n)/(n-1) # Corrected two-pass formula.
    sdev = np.sqrt(var)
    if var > 0.:
        skew /= (n*var*sdev);
        curt=(curt)/(n*var*var)-3.0;
    else:
        print "No skew/kurtosis when variance = 0 (in moment)"
        exit(1)

    return ave, adev, sdev, var, skew, curt
