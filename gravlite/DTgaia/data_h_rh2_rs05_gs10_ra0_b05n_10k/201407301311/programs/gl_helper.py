#!/usr/bin/env ipython3

##
# @file
# all helper functions

# (c) 2013 Pascal S.P. Steger, psteger@phys.ethz.ch
import sys, traceback, pdb#, ipdb
import numpy as np
from scipy.interpolate import splrep, splev, interp1d
from scipy.integrate import quad
import gl_plot as gpl
from pylab import *

def myfill(x, N=3):
    if x==np.inf or x==-np.inf:
        return "inf"
    return str(int(x)).zfill(N)
## \fn myfill(x, N=3)
#  print integer with leading zeros
# @param x integer
# @param N number of paddings
# @return 00x

def ipol_rhalf_log(X, Y, rhalf):
    for i in range(len(X)):
        if X[i] < rhalf:
            lowr = X[i]
            lowi = i
        if X[len(X)-i-1] > rhalf:
            highr = X[len(X)-i-1]
            highi = len(X)-i-1
    dr = np.log(highr) - np.log(lowr)
    dy = np.log(Y[highi]) - np.log(Y[lowi])

    Yhalf = np.exp(np.log(Y[lowi]) + (np.log(rhalf)-np.log(lowr))*dy/dr)
    return Yhalf
# \fn ipol_rhalf_log(X, Y, rhalf)
# interpolate function in loglog space to value in between data points,
# used for determination of Sig(rhalf)
# @param X radii, usually gp.xipol in our case, in [pc]
# @param Y values, usually Sig(gp.xipol) here, in [arb. units]
# @param rhalf half-light radius in [pc]


def print_summary(Rscale, Rc):
    print('Rscale = ', Rscale, ' pc')
    print('max(R) = ', max(Rc), ' pc')
    print('last element of R : ', Rc[-1], ' pc')
    print('total number of stars: ', len(Rc))
    return
## \fn print_summary(Rscale, Rc)
# print Radius information
# @param Rscale
# @param Rc


def draw_random_subset(x, ntracer):
    ind = np.arange(len(x))
    np.random.shuffle(ind)     # random.shuffle already changes ind
    return ind[:ntracer]
## \fn draw_random_subset(x, ntracer)
# take an array, and return shuffled array of ntracer entries
# no entry will be repeated


def add_errors(R, error):
    return R + error*np.random.randn(len(R))
## \fn add_errors(R, error)
# add procentual error to any quantity
# @param R array
# @param error in percent


def quadinflog(x, y, A, B, stop = False):
    if max(y) <= 0:
        print('integration over 0')
        return 0
    # work in log(y) space to enable smoother interpolating functions
    # scale to avoid any values <= 1
    minlog = min(np.log(y[y>0])) # exclude y==0 to circumvent error
    shiftlog = np.exp(1.5-minlog) # 1.-minlog not possible
    # otherwise we get divergent integrals for high radii (low length of x and y)
    
    # replace 0 values with 1e-XX
    yshift = y * shiftlog
    for i in range(len(yshift)):
        if yshift[i] == 0:
            yshift[i] = 10**(-10*i)

    splpar_nul = splrep(x, np.log(yshift), k=1, s=0.1)
    maxy = max(y)
    # invexp = lambda x: min(maxy, np.exp(splev(x,splpar_nul))/shiftlog)
    invexp = lambda x: np.exp(splev(x, splpar_nul))/shiftlog
    #invexparr= lambda x: np.minimum(maxy*np.ones(len(x)), np.exp(splev(x, splpar_nul))/shiftlog)
    dropoffint = quad(invexp, A, B, epsrel=1e-4, limit=100, full_output=1)
    # for debugging:
    # dropoffint = quad(invexp, A, B, epsrel=1e-3, limit=50, full_output=1)
    if stop and len(dropoffint)>3:
        print('warning by quad in quadinflog')
    return dropoffint[0]
## \fn quadinflog(x, y, A, B)
# integrate y over x, using splines
# @param x free variable
# @param y integrand
# @param A left boundary
# @param B right boundary
# @param stop = True


def quadinfloglog(x, y, A, B, stop = True):
    # work in log(x), log(y) space to enable smoother interpolating functions
    # scale to avoid any values <= 1
    # otherwise we get divergent integrals for high radii (low length of x and y)
    splpar_nul = splrep(np.log(x), np.log(y), k=1, s=0.1) # tunable k=2; s=0, ..
    invexp = lambda x: np.exp(splev(np.log(x), splpar_nul))
    # invexp = lambda x: min(y[0], np.exp(splev(np.log(x), splpar_nul)))
    dropoffint = quad(invexp, A, B, epsrel=1e-3, limi=100, full_output=1)
    if stop and len(dropoffint)>3:
        print('warning in quad in quadinfloglog')
    return dropoffint[0]
## \fn quadinflog(x, y, A, B)
# integrate y over x, using splines
# @param x free variable
# @param y integrand
# @param A left boundary
# @param B right boundary


def quadinflog2(x, y, A, B):
    # work in logarithmic space to enable smoother interpolating functions
    # scale to avoid any values <= 1 (such that a second log step can be taken)
    minlog = min(np.log(y))
    shiftlog = np.exp(1.5-minlog) # 1.-minlog not possible,
    # otherwise we get divergent integrals for high radii (low length of x and y)
    
    splpar_nul = splrep(x,np.log(np.log(y*shiftlog)), k=1, s=0.1)
    invexp = lambda x: min(y[0], np.exp(np.exp(splev(x, splpar_nul)))/shiftlog)
    dropoffint = quad(invexp, A, B)
    return dropoffint[0]
## \fn quadinflog2(x, y, A, B)
# integrate y over x, using splines in log(log(y)) space not used
# @param x free variable
# @param y integrand
# @param A left boundary
# @param B right boundary
# @return integrated y


def checknan(arr, place=''):
    if np.isnan(np.sum(arr)):
        print('NaN found! '+place)
        raise Exception('NaN', 'found')
        # not executed anymore :)
        traceback.print_tb(sys.exc_info()[2])
        return True
    else:
        return False
## \fn checknan(arr, place):
# NaN encountered at any position in arr?
# @param arr array of float values
# @param place = '' show user where to search
# @return True if NaN found


def checkpositive(arr, place=''):
    if checknan(arr, place):
        return True
    if min(arr) < 0.:
        print(place+" < 0 !")
        # traceback.print_tb(sys.exc_info()[2])
        raise Exception('negative','found')
        return True
    if max(arr) >= np.inf:
        print('infinity in checkpositive '+place)
        raise Exception('infinity', 'found')
        return True
    else:
        return False
## \fn checkpositive(arr, place)
# abort if negative value found
# @param arr array of float values
# @param place = '' tell user what to do better
# @return Exception if (non-physical) negative value found



def expolpowerlaw(R0, Sigdat, Rnuright, minp = -2.001):
    alpha = (np.log(Sigdat[-3])-log(Sigdat[-1]))/(np.log(R0[-3])-np.log(R0[-1]))
    alpha = min(alpha, minp) # assert finite mass
    Sig0  = Sigdat[-1]
    logSigdatright = np.log(Sig0)+alpha*(np.log(Rnuright)-np.log(R0[-1]))
    return np.exp(logSigdatright)
## \fn expolpowerlaw(R0, Sigdat, Rnuright, minp)
# extrapolate to high radii, assuming a minimum powerlaw
# @param R0 starting radii
# @param Sigdat starting profile
# @param Rnuright extrapolation radii
# @param minp minimum powerlaw .. extrapolation has to fall at least as much


def complete_nu(R0, Sigdat, Sigerr, Rnu):
    Rnuleft = Rnu[Rnu<R0[0]]   # extension of radii to the left
    Rnuright = Rnu[Rnu>R0[-1]] # extension of radii to the right
    R0nu = Rnu[(Rnu>=R0[0]) * (Rnu<=R0[-1])] # [pc]

    # TODO: use powerlaw to smaller radii
    Sigdatleft = np.exp(expol(R0, np.log(Sigdat), Rnuleft, 'linear'))
    Sigerrleft = (Sigerr[0]/Sigdat[0])*Sigdatleft

    Sigdatnu = linipollog(R0, Sigdat, R0nu)
    Sigerrnu = linipollog(R0, Sigerr, R0nu)

    #Sigdatright = np.exp(expol(R0, np.log(Sigdat), Rnuright, 'linear'))
    #Sigerrright = (Sigerr[-1]/Sigdat[-1])*Sigdatright
    # powerlaw extension to highest radii from last 3 binned points
    Sigdatright = expolpowerlaw(R0, Sigdat, Rnuright)
    Sigerrright = (Sigerr[-1]/Sigdat[-1])*Sigdatright

    Sigdatnu = np.hstack([Sigdatleft, Sigdatnu, Sigdatright])
    Sigerrnu = np.hstack([Sigerrleft, Sigerrnu, Sigerrright])
    #ipdb.set_trace()

    return Sigdatnu, Sigerrnu
## \fn complete_nu(R0, Sigdat, Sigerr, Rnu)
# inter- and extrapolate Sigdat and Signu to Rnu (gp.xfine, generally)
# @param R0 radii in [pc]
# @param Sigdat surface density
# @param Sigerr error
# @param Rnu new radii [pc]


def derivative(f):
    def df(x, h=0.1e-5):
        return ( f(x+h/2) - f(x-h/2) )/h
    return df
## \fn derivative(f)
# Computes the numerical derivative of a function.
# @param f function (ex. sin  or  cos)

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
## \fn deriv(y,x)
# gives numeric differentiation, assumes constant x spacings
# @param y dependent variable, array of size N
# @param x free variable, array of same size N
# @return dy/dx

def derivcoarse(y,x):
    y0 = ipol(x[:3],y[:3],x[0]-(x[1]-x[0])) #[y]
    y = np.hstack([y0, y]) #[y]
    x = np.hstack([x[0]-(x[1]-x[0]), x]) #[x]
    return (y[1:]-y[:-1])/(x[1:]-x[:-1])
## \fn derivcoarse(y,x)
# numerical derivative, with h=const=stepsize of radial bins
# @param y dependent variable, array of size N
# @param x free variable, array of same size N
# @return quotient of differences


def derivipol(y,x):
    step = x[1]-x[0]
    rbf = Rbf(x, y, smooth=1.5e-9)

    dydx = np.zeros(len(x))
    import scipy.misc as sm
    for i in range(len(x)):
        dydx[i] = sm.derivative(rbf, x[i], dx=step, n=1, args=(), order=3)
    return dydx
## \fn derivipol(y,x)
# numerical derivated with function from interpolation
# @param y dependent variable, array of size N
# @param x free variable, array of same size N


def err(a, gp):
    import numpy.random as npr
    return gp.err/a*(1.+npr.rand()/100.)
## \fn err(a)
# penalize prior violations with almost constant low likelihood
# a small change makes sure that if all nlive points in MultiNest
# get the same error, it does not abort there

def expDtofloat(strun):
    return float(strun.decode('utf-8').replace('D','E'))
## \fn expDtofloat(s)
# replace D floats with E floats
# @param s string
# @return s in 'xxxExx' format


def extend_radii(r0):
    r0ext = np.array([r0[0]*0.25, r0[0]*0.50, r0[0]*0.75])
    dR = r0[1:] - r0[:-1]
    r0nu = np.hstack([r0ext, r0, dR/2.+r0[:-1]])
    r0nu.sort()
    return r0nu
## \fn extend_radii(r0)
# insert equally spaced points at 1/4, 1/2, and 3/4 of a bin span
# @param r0 original bin centers


def readcol3(filena):
    a,b,c = np.loadtxt(filena,skiprows=1,unpack=True)
    return a,b,c
## \fn readcol3(filena)
# read in 3 columns of data
# @param filena filename
# @return a,b,c: columns of a 3-column file


def readcol5(filena):
    a,b,c,d,e = np.loadtxt(filena,skiprows=1,unpack=True)
    return a,b,c,d,e
## \fn readcol5(filena)
# read in 5 columns of data
# @param filena filename
# @return 5 columns


def readcoln(filena):
    return np.loadtxt(filena,skiprows=1,unpack=True)
## \fn readcoln(filena)
# read in n columns of data
# @param filena filename
# @return array

    
def pretty(arr,dig=3):
    return ("%."+str(dig)+"f")%arr
## \fn pretty(arr,dig=3)
# clip floats after dig=3 digits
# @param arr array of floats
# @param dig number of digits, default is 3


from scipy.interpolate import Rbf, InterpolatedUnivariateSpline
def ipol(xin,yin,xout,smooth=1.e-9):
    if np.isnan(np.sum(yin)):
        print('NaN found in interpolation! Go check where it occured!')
    rbf = Rbf(xin, yin, smooth=smooth)
    return rbf(xout)
## \fn ipol(xin, yin, xout, smooth=1.e-9)
# interpolate function in lin space, based on radial basis functions
# @param xin free variable array
# @param yin dependent array
# @param xout interpolation points
# @param smooth smoothing number. very small as a default. can be enhanced for smaller derivatives


def ipollog(xin, yin, xout, smooth=1.e-9):
    if min(yin)<0.:
        print('negative value encountered in ipollog, working with ipol now')
        return ipol(xin,yin,xout,smooth)
    
    rbf = Rbf(xin, np.log10(yin), smooth=smooth)
    return 10**(rbf(xout))
## \fn ipollog(xin, yin, xout, smooth=1.e-9)
# interpolate function in log10 space, based on radial basis functions
# @param xin free variable array
# @param yin dependent array
# @param xout interpolation points
# @param smooth smoothing number. very small as a default. can be enhanced for smaller derivatives


def linipollog(xin, yin, xout):
    f = interp1d(np.log(xin), np.log(yin), kind='linear')
    # sometimes, we have very small numerics giving slightly too
    # big xout (not in range of xin anymore)
    # here we correct that by setting the boundaries equal
    if np.abs(xin[0]-xout[0])/xin[0] < xin[0]/1e6:
        xout[0] = xin[0]
    if np.abs(xin[-1]-xout[-1])/xin[-1] < xin[-1]/1e6:
        xout[-1] = xin[-1]
    return np.exp(f(np.log(xout)))
## \fn linipollog(xin, yin, xout)
# interpolate linearly in ln-ln-space
# @param xin input positions in linear space
# @param yin input function in linear space
# @param xout output positions in linear space


def expol(xin, yin, xout, fn='multiquadric'):
    rbf = Rbf(xin, yin, function=fn)
    return rbf(xout)
## \fn expol(xin, yin, xout, fn)
# extrapolate data to xout, based on radial basis functions
# @param xin free variable array
# @param yin dependent array
# @param xout interpolation points
# @param function of interpolation, see http://docs.scipy.org/doc/scipy-0.7.x/reference/generated/scipy.interpolate.Rbf.html


def wait():
    input("Press Enter to continue...")
    return
## \fn wait()
# wait for key press


def bin_r_linear(rmin, rmax, nbin):
    binlength = (rmax - rmin)/(1.*nbin) #[rscale]
    binmin = np.zeros(nbin);  binmax = np.zeros(nbin)
    rbin = np.zeros(nbin)
    for k in range(nbin):
        binmin[k] = rmin+k*binlength #[rscale]
        binmax[k] = rmin+(k+1)*binlength #[rscale]
        rbin[k]   = binmin[k]+0.5*binlength #[rscale]
    return binmin, binmax, rbin
## \fn bin_r_linear(rmin, rmax, nbin)
# split interval linearly into bins
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
## \fn bin_r_log(rmin, rmax, nbin)
# split interval into bins, in a logarithmic fashion
# @param rmin starting point of interval
# @param rmax endpoint of interval
# @param nbin number of bins
# @return arrays of (beginning of bins, end of bins, position of bins)


def bin_r_const_tracers(r0, no):
    # procedure: get all particles in bin i
    #            get minimum, maximum radius. get radius of min/max before/after bin i
    #            get mean of (half of max bin/min next bin) for bin radius
    order = np.argsort(r0)
    r0 = np.array(r0)[order]

    no = int(no)
    
    mini = list(range(1,  len(r0), no))
    maxi = list(range(no, len(r0), no))
    maxi.append(len(r0))
    mini = np.array(mini)-1
    maxi = np.array(maxi)-1
    minri = [];   maxri = []
    nbin = int(1.*len(r0)/no)
    for i in range(nbin):
        minri.append(r0[mini[i]])
        maxri.append(r0[maxi[i]])

    midri = []
    for i in range(nbin):
        midri.append((minri[i]+maxri[i])/2.)
    return minri, maxri, midri
## \fn bin_r_const_tracers(r0, no)
# split interval into bins of constant particle number
# @param r0 radii from all particles in an array
# @param no integer, number of bins
# @return arrays of (beginning of bins, end of bins, position of bins)



def sort_profiles_binwise(profs):
    for i in range(len(profs)):
        # sort all mass models bin by bin
        profs[i] = np.sort(profs[i])
    return profs
## \fn sort_profiles_binwise(profs)
# take array of profiles, sort along each bin
# @param profs array [prof1, prof2, ..., profN]
# @return [bin1 bin2 .. binN]_min .. [bin1 bin2 .. binN]_max


def get_median_1_2_sig(profs):
    bins = len(profs)
    Mmedi = np.zeros(bins); Mmax  = np.zeros(bins); Mmin  = np.zeros(bins)
    M95hi = np.zeros(bins); M95lo = np.zeros(bins)
    M68hi = np.zeros(bins); M68lo = np.zeros(bins)
    mlen = len(profs[0])
    for i in range(bins):
        Mmax[i]  = profs[i, mlen-1]
        M95hi[i] = profs[i, 0.95 * mlen]
        M68hi[i] = profs[i, 0.68 * mlen]
        Mmedi[i] = profs[i, 0.50 * mlen]
        M68lo[i] = profs[i, 0.32 * mlen]
        M95lo[i] = profs[i, 0.05 * mlen]
        Mmin[i]  = profs[i, 0]
    return Mmin, M95lo, M68lo, Mmedi, M68hi, M95hi, Mmax
## \fn get_median_1_2_sig(profs)
# return envelopes of profiles at median, min, max, and 1, 2 sigma interval for a range of profiles
# @param profs bin-wise sorted profiles
# @return  Mmin, M95lo, M68lo, Mmedi, M68hi, M95hi, Mmax


def binsmooth(r, array, low, high, nbin, nanreplace):
    # sort r and array in ascending r order
    index = np.argsort(r)
    r = r[index]
    array = array[index]

    # determine linearly spaced bins
    binsize = (high-low)/(1.*nbin)
    rout = np.arange(nbin)*binsize + low
    binmin = rout - binsize/2.
    binmax = rout + binsize/2.

    # prepare run
    arrayout = np.zeros(nbin)
    arrayout1 = np.zeros(nbin)
    arrayout2 = np.zeros(nbin)
    count_bin = np.zeros(nbin)
    j=0
    siz = len(r)
    
    for i in range(nbin):
        count = 0
        while (binmax[i] > r[j]):
            arrayout1[i] = arrayout1[i]+array[j]
            arrayout2[i] = arrayout2[i]+array[j]**2
            if (array[j] != 0):
                count = count + 1
            if (j < siz-1):
                j = j + 1
            else:
                break

        if (count > 0):
            arrayout[i] = np.sqrt(arrayout2[i]/count-(arrayout1[i]/count)**2) # def of sig
        else:
            arrayout[i] = nanreplace
        count_bin[i] = count
    return binmin, binmax, rout, arrayout, count_bin
## \fn binsmooth(r, array, low, high, nbin, nanreplace)
# This routine takes an array(r) and bins it in r bins of size bin, def. sigma from array in bins
# @param r [pc]
# @param array corresponding values for e.g. v_z
# @param low [pc] lower bound on z range
# @param high [pc] high bound on z range
# @param nbin int, number of bins
# @param nanreplace if there is no data in a particular bin, it assigns a value of array(r)=nanreplace
# @return sqrt(avg(x**2)-avg(x)**2) from each bin, with dim [array]


def bincount(r, rmax):
    # sort radial bins
    index = np.argsort(r)
    r = r[index]

    # prepare for run
    nbin = len(rmax)
    arrayout  = np.zeros(nbin)
    count_bin = np.zeros(nbin)
    std       = np.zeros(nbin)
    j = 0
    siz = len(r)

    for i in range(nbin):
        while (rmax[i] > r[j]):
            arrayout[i] = arrayout[i] + 1.
            if (j < siz-1):
                j = j + 1
            else:
                break
        count_bin[i] = arrayout[i]
        std[i] = np.sqrt(count_bin[i])
    
    return arrayout, std
## \fn bincount(r, rmax)
# take an array, r, and count the number of elements in r bins of size bin
# @param r array of floats
# @param rmax upper bound of bins
# @return arrrayout, std
