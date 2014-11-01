#!/usr/bin/python

from scipy import *
from scipy.optimize import leastsq
from scipy.optimize import curve_fit

import pylab
import sys
import numpy as np

def residuals(p, y, x): 
	err = y-peval(x,p) 
	return err

def peval(x, p0,p1,p2): 
	#return p[0]*(1-exp(-(p[2]*x)**p[4])) + p[1]*(1-exp(-(p[3]*(x))**p[5] ))
        return np.log10(p0*(1+x**2/p1**2)**(-5/2)/p1**3+p2)

if(len(sys.argv)!=2):
	print 'use: fit_plummer.py prof.dat'
	exit(1)
	
filename=sys.argv[1]
data = loadtxt(filename)

y     = np.log10(data[:,2])
x     = (data[:,0]+data[:,1])/2*1000
sigma = np.log10(data[:,3])

M=10
b=0.1
off=1
pname = (['M','b','off'])
p0 = array([M, b,off])
#plsq = leastsq(residuals, p0, args=(y, x), maxfev=2000)
popt,pcov = curve_fit(peval,x,y,p0,sigma,maxfev=2000)
print pcov
pylab.plot(x,y,'x')#,'title "Meas" with points')
#pylab.errorbar(x,y,sigma)
pylab.plot(x,peval(x,popt[0],popt[1],popt[2]))#,'title "Fit" with lines lt -1')
pylab.plot(x,peval(x,popt[0],popt[1],popt[2]))#,'title "Fit" with lines lt -1')
pylab.savefig('kinfit.png')

print "Final parameters"
print popt
