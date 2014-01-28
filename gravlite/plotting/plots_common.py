#!/usr/bin/env python3

##
# @file
# common functions needed in plot_* programs

# (c) 2013 ETHZ, Pascal Steger, psteger@phys.ethz.ch

from pylab import subplots, savefig, show, ion, ioff

def M_anf(r):
    return r**2/(r+1.)**2
## \fn M_anf(r)
# analytic mass profile
# @param r radius in [pc]


def ipol(xin,yin,xout,smooth=1.e-9):
    rbf = Rbf(xin, yin, smooth=smooth)
    return rbf(xout)
## \fn ipol(xin, yin, xout, smooth=1.e-9)
# interpolation with radial basis function
# @param xin free variable
# @param yin dependant variable
# @param xout new free variables
# @param smooth smoothing parameter, 1e-9 is small


def ipollog(xin,yin,xout,smooth=1.e-9):
    rbf = Rbf(xin, np.log10(yin), smooth=smooth)
    return 10**(rbf(xout))
## \fn ipollog(xin, yin, xout, smooth=1.e-9)
# interpolation in logarithmic space
# @param xin free variable
# @param yin dependant variable
# @param xout new free variable
# @param smooth smoothin parameter, 1e-9 is small


def prepare_plots():
    global f,ax1
    ion()
    f, ax1 = subplots(1,1)
    draw()
    return f,ax1
## \fn prepare_plots()
# set up plotting window
# @return figure and axis object



def setlabels(ax, xtext, ytext):
    ax.set_xlabel(r'$'+xtext+'$')
    ax.set_ylabel(r'$'+ytext+'$')
## \fn setlabels(ax, xtext, ytext)
# set labels of axes
# @param ax axis object
# @param xtext label of x axis
# @param ytext label of y axis


def setlims(ax, xlim, ylim):
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
## \fn setlims(ax, xlim, ylim)
# set limits on axes objects
# @param ax axis object
# @param xlim 2D array: limits on x axis
# @param ylim 2D array: limits on y axis


def save_profile(basename, prof, M95lo, M68lo, Mmedi, M68hi, M95hi):
    fout = open(basename+".prof"+prof+".conf",'w')
    print(M95lo,file=fout)
    print(M68lo,file=fout)
    print(Mmedi,file=fout)
    print(M68hi,file=fout)
    print(M95hi,file=fout)
    fout.close()
    return
## \fn save_profile(basename, prof, M95lo, M68lo, Mmedi, M68hi, M95hi)
# save profile sigma data to file
# @param basename string
# @param prof string: profile to be saved
# @param M95lo
# @param M68lo
# @param Mmedi
# @param M68hi
# @param M95hi
    

def save_plot(nam):
    savefig(nam)
    return
## \fn save_plot(nam)
# save plot to file
# @param nam filename


def show_plots():
    ion()
    show()
    ioff()
    return
## \fn show_plots()
# stop program, show plots


def readcol(filena):
    a,b,c = np.loadtxt(filena,skiprows=1,unpack=True)
    return a,b,c
## \fn readcol(filena)
# read 3 columns of data
# @param filena filename
