#!/usr/bin/python
'''all helper functions from gl_funs'''
import numpy as np
import gl_helper as gh

def myfill(x):
    return str(int(x)).zfill(3)

def derivative(f):
    """
    Computes the numerical derivative of a function.
    """
    def df(x, h=0.1e-5):
        return ( f(x+h/2) - f(x-h/2) )/h
    return df

def deriv(y,x):
    '''gives numeric differentiation'''
    # assumes constant x spacings
    step = x[1]-x[0]
    xnew = []
    xnew.append(x[0]-step/2.)
    for i in range(len(x)):
        xnew.append(x[i]+step/2.)

    ynew = ipol(x,y,xnew)
    
    dydx = []
    for i in range(len(x)):
        dydx.append((ynew[i+1]-ynew[i])/step)

    return dydx

def derivcoarse(y,x):
    y0 = gh.ipol(x[:3],y[:3],x[0]-(x[1]-x[0]))
    y = np.hstack([y0, y])
    x = np.hstack([x[0]-(x[1]-x[0]), x])
    return (y[1:]-y[:-1])/(x[1:]-x[:-1])

def expDtofloat(s):
    return float(s.replace('D','E'))

def readcol(filena):
    a,b,c = np.loadtxt(filena,skiprows=1,unpack=True)
    return a,b,c

def readcoln(filena):
    return np.loadtxt(filena,skiprows=1,unpack=True)

def pretty(arr):
    if(type(arr) is 'float'):
        return "%.2f"%arr
    else:
        return "%.2f "*len(arr) % tuple(arr)

from scipy.interpolate import Rbf, InterpolatedUnivariateSpline

def smooth(xin,yin,smooth=1.e-9):
    return ipol(xin,yin,xin,smooth)

def smoothlog(xin,yin,smooth=1.e-9):
    return ipollog(xin,yin,xin,smooth)

def ipol(xin,yin,xout,smooth=1.e-9):
    rbf = Rbf(xin, yin, smooth=smooth)
    return rbf(xout)

def ipollog(xin,yin,xout,smooth=1.e-9):
    rbf = Rbf(xin, np.log10(yin), smooth=smooth)
    return 10**(rbf(xout))

def expol(xin,yin,xout):
    rbf = Rbf(xin, yin)
    return rbf(xout)

def wait():
    raw_input("Press Enter to continue...")
    return

