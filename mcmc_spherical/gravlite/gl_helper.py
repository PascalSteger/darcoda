#!/usr/bin/python
# (c) 2013 Pascal S.P. Steger
'''all helper functions from gl_funs'''

import numpy as np
import pdb
import gl_file as gf
import gl_plot as gpl


def myfill(x):
    return str(int(x)).zfill(3)


def checknan(arr):
    if np.isnan(np.sum(arr)):
        print('NaN found! Go check where it occured!')
        gf.get_working_pars()
        # TODO: jump to next main iteration in gravlite.py
    return True



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
    y0 = ipol(x[:3],y[:3],x[0]-(x[1]-x[0])) #[y]
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



def pretty(arr,dig=3):
    'clip floats after dig=3 digits'
    return ("%."+str(dig)+"f")%arr



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
        print >> 'negative value encountered in ipollog, working with ipol now'
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




def bin_r_linear(rmin, rmax, nbin):

    binlength = (rmax - rmin)/(1.*nbin) #[rcore]
    # print >> 'binlength [rcore] = ', binlength
    binmin = np.zeros(nbin);  binmax = np.zeros(nbin)
    rbin = np.zeros(nbin)
    for k in range(nbin):
        binmin[k] = rmin+k*binlength #[rcore]
        binmax[k] = rmin+(k+1)*binlength #[rcore]
        rbin[k]   = binmin[k]+0.5*binlength #[rcore]
    return binmin, binmax, rbin




def bin_r_log(rmin, rmax, nbin):

    bins   = np.logspace(np.log10(rmin), np.log10(rmax), nbin+1, endpoint=True)
    binmin = bins[:-1]; binmax = bins[1:]
    rbin = np.zeros(nbin)
    for k in range(nbin):
        rbin[k]   = np.logspace(np.log10(binmin[k]),np.log10(binmax[k]),3, endpoint=True)[1]
    return binmin, binmax, rbin


def bin_r_const_tracers(r0,no): # radii, number of tracers
    # TODO HALO: specify number of particles, build bins accordingly
    # procedure: get all particles in bin i
    #            get minimum, maximum radius. get radius of min/max before/after bin i
    #            get mean of (half of max bin/min next bin) for bin radius
    order = np.argsort(r0)
    r0 = np.array(r0)[order]

    mini = range(1,len(r0),no)
    maxi = range(no,len(r0),no)
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

