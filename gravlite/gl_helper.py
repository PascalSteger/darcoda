#!/usr/bin/env ipython-python3.2

##
# @file
# all helper functions from gl_funs
# (c) 2013 Pascal S.P. Steger


import numpy as np
import pdb
import gl_file as gf
# import gl_plot as gpl
from scipy.interpolate import splrep, splev
from scipy.integrate import quad

## print integer with leading zeros
# @param x integer
# @param N number of paddings
# @return 00x
def myfill(x,N=3):
    return str(int(x)).zfill(N)

## integrate y over x, using splines
# @param x free variable
# @param integrand
# @param A left boundary
# @param B right boundary
def quadinf(x, y, A, B):
    # work in logarithmic space to enable smoother interpolating functions
    # scale to avoid any values <= 1
    minlog = min(np.log(y))
    shiftlog = np.exp(1.5-minlog) # 1.-minlog not possible, other wise we get divergent integrals for high radii (low length of x and y)
    tcknulog = splrep(x,np.log(y*shiftlog),k=1,s=0.1)
    invexp = lambda x: min(y[0],np.exp(splev(x,tcknulog,der=0))/shiftlog)
    # pdb.set_trace() # TADA, no new errors :)
    dropoffint = quad(invexp, A, B)
    return dropoffint[0]


## integrate y over x, using splines in log(log(y)) space
# @param x free variable
# @param integrand
# @param A left boundary
# @param B right boundary
def quadinfloglog(x, y, A, B):
    # work in logarithmic space to enable smoother interpolating functions
    # scale to avoid any values <= 1 (such that a second log step can be taken)
    minlog = min(np.log(y))
    shiftlog = np.exp(1.5-minlog) # 1.-minlog not possible, other wise we get divergent integrals for high radii (low length of x and y)
    tcknulog = splrep(x,np.log(np.log(y*shiftlog)),k=1,s=0.1)
    invexp = lambda x: min(y[0],np.exp(np.exp(splev(x,tcknulog,der=0)))/shiftlog)
    # gpl.plot(x,y)
    # xnew = max(x)+np.arange(0.,max(x),max(x)/30.)
    # gpl.plot(xnew,invexp(xnew))
    # pdb.set_trace()
    dropoffint = quad(invexp, A, B)
    return dropoffint[0]

## abort if NaN encountered
# @param arr array of float values
# @return Exception if NaN found
def checknan(arr):
    if np.isnan(np.sum(arr)):
        print('NaN found! Go check where it occured!')
        gf.get_working_pars()
        raise Exception('NaN', 'found')
        return True
    else:
        return False

## abort if negative value found
# @param arr array of float values
# @return Exception if (non-physical) negative value found
def checkpositive(arr):
    if checknan(arr):
        print('found non-physical NaN value!')
        return True
    if min(arr) <= 0.:
        print('found non-physical negative value!')
        gf.get_working_pars()
        raise Exception('negative','found')
        return True
    else:
        return False

##  Computes the numerical derivative of a function.
# @param f function (ex. sin  or  cos)
def derivative(f):
    def df(x, h=0.1e-5):
        return ( f(x+h/2) - f(x-h/2) )/h
    return df

## gives numeric differentiation, assumes constant x spacings
# @param y dependent variable, array of size N
# @param x free variable, array of same size N
# @return dy/dx
def deriv(y,x):
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

    return np.array(dydx)

##  numerical derivative, with h=const=stepsize of radial bins
# @param y dependent variable, array of size N
# @param x free variable, array of same size N
# @return quotient of differences
def derivcoarse(y,x):
    y0 = ipol(x[:3],y[:3],x[0]-(x[1]-x[0])) #[y]
    y = np.hstack([y0, y]) #[y]
    x = np.hstack([x[0]-(x[1]-x[0]), x]) #[x]
    return (y[1:]-y[:-1])/(x[1:]-x[:-1])

## numerical derivated with function from interpolation
# @param y dependent variable, array of size N
# @param x free variable, array of same size N
def derivipol(y,x):
    step = x[1]-x[0]
    rbf = Rbf(x, y, smooth=1.5e-9)

    dydx = np.zeros(len(x))
    import scipy.misc as sm
    for i in range(len(x)):
        dydx[i] = sm.derivative(rbf, x[i], dx=step, n=1, args=(), order=3)
    return dydx

## replace D floats with E floats
# @param s string
# @return s in 'xxxExx' format
def expDtofloat(s):
    return float(s.replace('D','E'))

## read in 3 columns of data
# @param filena filename
# @return a,b,c: columns of a 3-column file
def readcol3(filena):
    a,b,c = np.loadtxt(filena,skiprows=1,unpack=True)
    return a,b,c

## read in 5 columns of data
# @param filena filename
# @return 5 columns
def readcol5(filena):
    a,b,c,d,e = np.loadtxt(filena,skiprows=1,unpack=True)
    return a,b,c,d,e

## read in n columns of data
# @param filena filename
# @return array
def readcoln(filena):
    return np.loadtxt(filena,skiprows=1,unpack=True)

    
## clip floats after dig=3 digits
# @param arr array of floats
# @param dig number of digits, default is 3
def pretty(arr,dig=3):
    return ("%."+str(dig)+"f")%arr

from scipy.interpolate import Rbf, InterpolatedUnivariateSpline
## interpolate function in lin space, smooth it
# @param xin free variable array
# @param yin dependent array
# @smooth smoothing number. very small as a default. can be enhanced for smaller derivatives
def smooth(xin,yin,smooth=1.e-9):
    return ipol(xin,yin,xin,smooth)

## interpolate function in log space, smooth it. use only if yin[i]>0
# @param xin free variable array
# @param yin dependent array
# @smooth smoothing number. very small as a default. can be enhanced for smaller derivatives
def smoothlog(xin,yin,smooth=1.e-9):
    return ipollog(xin,yin,xin,smooth)


## interpolate function in lin space, based on radial basis functions
# @param xin free variable array
# @param yin dependent array
# @param xout interpolation points
# @smooth smoothing number. very small as a default. can be enhanced for smaller derivatives
def ipol(xin,yin,xout,smooth=1.e-9):
    if np.isnan(np.sum(yin)):
        print('NaN found! Go check where it occured!')
        pdb.set_trace()
    rbf = Rbf(xin, yin, smooth=smooth)
    return rbf(xout)

## interpolate function in log space, based on radial basis functions
# @param xin free variable array
# @param yin dependent array
# @param xout interpolation points
# @smooth smoothing number. very small as a default. can be enhanced for smaller derivatives
def ipollog(xin, yin, xout, smooth=1.e-9):
    if min(yin)<0.:
        print(file='negative value encountered in ipollog, working with ipol now')
        return ipol(xin,yin,xout,smooth)
    
    rbf = Rbf(xin, np.log10(yin), smooth=smooth)
    return 10**(rbf(xout))


## extrapolate data to xout, based on radial basis functions
# @param xin free variable array
# @param yin dependent array
# @param xout interpolation points
def expol(xin,yin,xout):
    rbf = Rbf(xin, yin)
    return rbf(xout)

## wait for key press
def wait():
    input("Press Enter to continue...")
    return

## split interval linearly into bins
# @param rmin starting point of interval
# @param rmax endpoint of interval
# @param nbin number of bins
# @return arrays of (beginning of bins, end of bins, position of bins)
def bin_r_linear(rmin, rmax, nbin):
    binlength = (rmax - rmin)/(1.*nbin) #[rcore]
    binmin = np.zeros(nbin);  binmax = np.zeros(nbin)
    rbin = np.zeros(nbin)
    for k in range(nbin):
        binmin[k] = rmin+k*binlength #[rcore]
        binmax[k] = rmin+(k+1)*binlength #[rcore]
        rbin[k]   = binmin[k]+0.5*binlength #[rcore]
    return binmin, binmax, rbin


## split interval into bins, in a logarithmic fashion
# @param rmin starting point of interval
# @param rmax endpoint of interval
# @param nbin number of bins
# @return arrays of (beginning of bins, end of bins, position of bins)
def bin_r_log(rmin, rmax, nbin):
    bins   = np.logspace(np.log10(rmin), np.log10(rmax), nbin+1, endpoint=True)
    binmin = bins[:-1]; binmax = bins[1:]
    rbin = np.zeros(nbin)
    for k in range(nbin):
        rbin[k]   = np.logspace(np.log10(binmin[k]),np.log10(binmax[k]),3, endpoint=True)[1]
    return binmin, binmax, rbin

## split interval into bins of constant particle number
# @param r0 radii from all particles in an array
# @param no integer, number of bins
# @return arrays of (beginning of bins, end of bins, position of bins)
def bin_r_const_tracers(r0,no):
    # procedure: get all particles in bin i
    #            get minimum, maximum radius. get radius of min/max before/after bin i
    #            get mean of (half of max bin/min next bin) for bin radius
    order = np.argsort(r0)
    r0 = np.array(r0)[order]

    no = int(no)
    
    mini = list(range(1,len(r0),no))
    maxi = list(range(no,len(r0),no))
    maxi.append(len(r0))
    mini = np.array(mini)-1; maxi = np.array(maxi)-1
    minri = [];   maxri = []
    nbin = int(1.*len(r0)/no)
    for i in range(nbin):
        minri.append(r0[mini[i]])
        maxri.append(r0[maxi[i]])

    midri = []
    # TODO: check fine to take *not* half the first particle radius for midri
    for i in range(nbin):
        midri.append((minri[i]+maxri[i])/2.)
    return minri, maxri, midri

