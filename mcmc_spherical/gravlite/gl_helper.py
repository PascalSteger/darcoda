#!/usr/bin/python2.7
# (c) 2013 Pascal S.P. Steger
'''all helper functions from gl_funs'''

import numpy as np
import gl_helper as gh
import gl_plot as gpl

def myfill(x):
    return str(int(x)).zfill(3)



def checknan(arr):
    if np.isnan(np.sum(arr)):
        print 'NaN found! Go check where it occured!'
        import pdb
        pdb.set_trace()
    return



def derivative(f):
    'Computes the numerical derivative of a function.'
    def df(x, h=0.1e-5):
        return ( f(x+h/2) - f(x-h/2) )/h
    return df



def deriv(y,x):
    'gives numeric differentiation, assumes constant x spacings'
    step = x[1]-x[0]
    xnew = []
    flick = 5.
    for i in range(len(x)):
        xnew.append(x[i]-step/flick)
        xnew.append(x[i]+step/flick)

    ynew = ipol(x,y,xnew,smooth=1.5e-9)
    
    dydx = []
    for i in range(len(x)):
        dydx.append((ynew[2*i+1]-ynew[2*i])/(step/flick))

    return dydx




def derivcoarse(y,x): #[y], [x]
    'numerical derivative, with h=const=stepsize of radial bins'
    y0 = gh.ipol(x[:3],y[:3],x[0]-(x[1]-x[0])) #[y]
    y = np.hstack([y0, y]) #[y]
    x = np.hstack([x[0]-(x[1]-x[0]), x]) #[x]
    return (y[1:]-y[:-1])/(x[1:]-x[:-1]) #[y/x]




def derivipol(y,x):
    'numerical derivated with function from interpolation'
    step = x[1]-x[0]
    rbf = Rbf(x, y, smooth=1.5e-9)

    dydx = np.zeros(len(x))
    import scipy.misc as sm
    for i in range(len(x)):
        dydx[i] = sm.derivative(rbf, x[i], dx=step, n=1, args=(), order=3)
    return dydx




def expDtofloat(s):
    'replace D floats with E floats'
    return float(s.replace('D','E'))




def readcol(filena):
    'read in 3 columns of data'
    a,b,c = np.loadtxt(filena,skiprows=1,unpack=True)
    return a,b,c




def readcoln(filena):
    'read in n columns of data'
    return np.loadtxt(filena,skiprows=1,unpack=True)



def pretty(arr):
    'clip floats after 2 digits'
    if(type(arr) is 'float'):
        return "%.2f"%arr
    else:
        return "%.2f "*len(arr) % tuple(arr)



from scipy.interpolate import Rbf, InterpolatedUnivariateSpline

def smooth(xin,yin,smooth=1.e-9):
    'interpolate function in lin space, smooth it'
    return ipol(xin,yin,xin,smooth)


def smoothlog(xin,yin,smooth=1.e-9):
    'interpolate function in log space, smooth it. use only if yin[i]>0'
    return ipollog(xin,yin,xin,smooth)



def ipol(xin,yin,xout,smooth=1.e-9):
    'interpolate function in lin space'
    if np.isnan(np.sum(yin)):
        print 'NaN found! Go check where it occured!'
        import pdb
        pdb.set_trace()
    rbf = Rbf(xin, yin, smooth=smooth)
    return rbf(xout)



def ipollog(xin, yin, xout, smooth=1.e-9):
    'interpolate function in log space'
    if min(yin)<0.:
        print 'negative value encountered in ipollog, working with ipol now'
        return ipol(xin,yin,xout,smooth)
    
    rbf = Rbf(xin, np.log10(yin), smooth=smooth)
    return 10**(rbf(xout))



def expol(xin,yin,xout):
    'extrapolate data to xout'
    rbf = Rbf(xin, yin)
    return rbf(xout)



def wait():
    'wait for key press'
    raw_input("Press Enter to continue...")
    return

