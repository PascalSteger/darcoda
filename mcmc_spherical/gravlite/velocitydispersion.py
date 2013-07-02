#!/usr/bin/python
# calculate velocity dispersion of wedges
####################################################

import numpy
from BiWeight import meanbiweight

import sys
if(len(sys.argv)<2):
    print "use: velocitydispersion.py [car,scl,sex,for]"
    exit(1)
dwarf = sys.argv[1]
dir = "/home/ast/read/dark/dwarf_data/data_obs/"
dir = '/home/psteger/sci/dwarf_data/data_obs/'
print dir+dwarf+"/table_merged.bin"

#distance error, set for sim. TODO: change
rerror  = 0.1
vrerror = 0.1

# get radius, used for all binning
x,y,vlos = numpy.loadtxt(dir+dwarf+'/centerpos.txt',skiprows=1,unpack=True)
totmass = len(x)
r = numpy.sqrt(x*x+y*y)

#set binning
bins = 12
#nbins = (max - min)*N^(1/3)/(2*(Q3-Q1)) #(method of Wand)
rmin = 0.;   rmax = max(r)
binlength = (rmax-rmin)/bins
binmin = numpy.zeros(bins)
binmax = numpy.zeros(bins)
rbin=numpy.zeros(bins)
for i in range(bins):
    binmin[i] = rmin +     i*binlength
    binmax[i] = rmin + (i+1)*binlength
    rbin[i]   = binmin[i] + 0.5*binlength

rs = rerror*numpy.random.randn(len(r))+r
vlos = vrerror*numpy.random.randn(len(vlos))+vlos
d=open(dir+dwarf+'/velocitydispersionlos.txt','w')
print>>d,'r','sigma_r(r)','error'
d.close()
#iterations for drawing a given radius in bin
n = 1000
dispvelocity=numpy.zeros((bins,n))
a=numpy.zeros((bins,n))
p_dvlos = numpy.zeros(bins)
p_edvlos= numpy.zeros(bins)

for k in range(n):
    rsi = rerror*numpy.random.randn(len(rs))+rs
    vlosi = vrerror*numpy.random.randn(len(vlos))+vlos
    for i in range(bins):
        ind1 = numpy.argwhere(numpy.logical_and(rsi>binmin[i],rsi<binmax[i])).flatten()
        a[i][k] = 1.*len(ind1)
        vlos1 = vlosi[ind1]
        if(len(ind1)<=1):
            dispvelocity[i][k] = 0.
        else:
            dispvelocity[i][k] = meanbiweight(vlos1,ci_perc=68.4,ci_mean=True,ci_std=True)[1]

for i in range(bins):
    dispvel = 1.*numpy.sum(dispvelocity[i])/n
    ab = 1.*numpy.sum(a[i])/n
    if(ab == 0):
        dispvelerror = 0.0
    else:
        dispvelerror = dispvel/numpy.sqrt(ab)
    p_dvlos[i]  = dispvel
    p_edvlos[i] = dispvelerror

maxvlos = max(p_dvlos)
d=open(dir+dwarf+'/velocitydispersionlos.txt','a')
for i in range(bins):
    print>>d,rbin[i],p_dvlos[i],p_edvlos[i]
d.close()

from pylab import *
ion(); subplot(111)
print 'rbin = ',rbin
print 'p_dvlos = ',p_dvlos
print 'p_edvlos = ',p_edvlos
plot(rbin,p_dvlos,'b',linewidth=3)
fill_between(rbin,p_dvlos-p_edvlos,p_dvlos+p_edvlos,alpha=0.5,color='r')

# vertical line
#axvline(x=4.0)
######################################################################################
# overall line-of-sight velocity of the whole dwarf galaxy, in [km/s]
kms = 1
vLOS= {
    'for': lambda x: x * (53),#+/- 3
    'car': lambda x: x * (224),#+/- 3
    'sex': lambda x: x * (227), #+/- 3
    'scl': lambda x: x * (108)  #+/- 3
    }[dwarf](kms)
print 'vLOS = ',vLOS
#axhline(y=vLOS)

xlabel(r'$r [pc]$')
ylabel(r'$\langle\sigma_{LOS}\rangle [km/s]$')
ylim([-5,30])
#plt.legend(['\rho','\rho'],'lower left')
title(dwarf)
savefig(dir+dwarf+"/vlosdispersion.png")
ioff();show()


