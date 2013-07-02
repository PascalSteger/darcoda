#!/usr/bin/python
import numpy as np
import pdb
# calculate approximative center of mass, assuming constant stellar mass
# read in [some].bin.[MV,V-I], output ID,Mstar

#choose simulation

import sys
if(len(sys.argv)<2):
    print "use: centreofmass.py [car,scl,sex,for]"
    exit(1)
    
dwarf = sys.argv[1] #TODO: read from first command line parameter
dir = "/home/ast/read/dark/dwarf_data/data_obs/"
dir = '/home/psteger/sci/dwarf_data/data_obs/'
print dir+dwarf+"/table_merged.bin"

delim = [0,22,3,3,6,4,3,5,6,6,7,5,6,5,6,5,6]
ID = np.genfromtxt(dir+dwarf+"/table_merged.bin",skiprows=29,unpack=True,usecols=(0,1),delimiter=delim,dtype="string")
RAh,RAm,RAs,DEd,DEm,DEs,Vmag,VI,VHel,e_VHel,SigFe,e_SigFe,SigMg,e_SigMg,PM=np.genfromtxt(dir+dwarf+"/table_merged.bin",skiprows=29,unpack=True,usecols=tuple(range(2,17)),delimiter=delim,filling_values=-1)

# only use stars which are members of the dwarf
pm = (PM>=0.95)
print "fraction of members = ",1.*sum(pm)/len(pm)
ID=ID[1][pm]; RAh=RAh[pm]; RAm=RAm[pm]; RAs=RAs[pm]; DEd=DEd[pm]; DEm=DEm[pm]; DEs=DEs[pm]; Vmag = Vmag[pm]; VI=VI[pm]; VHel=VHel[pm]; e_VHel=e_VHel[pm]; SigFe=SigFe[pm]; e_SigFe=e_SigFe[pm]; SigMg=SigMg[pm]; e_SigMg=e_SigMg[pm];PM=PM[pm]



sig = abs(RAh[0])/RAh[0]
print 'RAh: signum = ',sig
RAh = RAh/sig
xs = 15*(RAh*3600+RAm*60+RAs)*sig       # [arcsec/15]

sig = abs(DEd[0])/DEd[0]
print 'DEd: signum = ',sig
DEd = DEd/sig
ys = (DEd*3600+DEm*60+DEs)*sig          # [arcsec]

# determine COM weighted by probability of membership
split = 0.99
ms = (PM-split)/(1-split)               # [1], thought of as propto [Msun]
for i in range(len(ms)):
    if ms[i]<0.: ms[i]=1.-split


centreofmassx = np.sum(xs*ms)/sum(ms)   # [arcsec]
centreofmassy = np.sum(ys*ms)/sum(ms)   # [arcsec]

print '(xc,yc)=',centreofmassx,centreofmassy

xsnew = xs-centreofmassx                # [arcsec]
ysnew = ys-centreofmassy                # [arcsec]


####################################################################################
# calculate v_{LOS} after subtracting bulk motion of dwarf
# overall line-of-sight velocity of the whole dwarf galaxy, in [km/s]
kms = 1
vLOS= {
    'for': lambda x: x * (53),#+/- 3   [km/s]
    'car': lambda x: x * (224),#+/- 3  [km/s]
    'sex': lambda x: x * (227), #+/- 3 [km/s]
    'scl': lambda x: x * (108)  #+/- 3 [km/s]
    }[dwarf](kms)
print 'vLOS = ',vLOS

vlos = VHel - vLOS                      # [km/s]
# TODO: error in VHel


arcsec = 60
rcore= {
    'for': lambda x: x * (13.8),#+/- 0.8   [arcsec]
    'car': lambda x: x * (8.8),#+/- 1.2    [arcsec]
    'sex': lambda x: x * (16.6), #+/- 1.2  [arcsec]
    'scl': lambda x: x * (5.8)  #+/- 1.6   [arcsec]
    }[dwarf](arcsec)
print 'rcore = ',rcore,' arcsec'

xsnew = xsnew/rcore                     # [rcore]
ysnew = ysnew/rcore                     # [rcore]

rhalf = {
    'for': lambda x: x * (339),#+/- 36   [pc]
    'car': lambda x: x * (137),#+/- 22   [pc]
    'sex': lambda x: x * (294),#+/- 38   [pc]
    'scl': lambda x: x * ( 94) #+/- 26   [pc]
    }[dwarf](1.)
print 'rhalf = ',rhalf,' pc'

# now determine rhalf from our dataset, assuming brightness propto Vmag

Vtot = np.sum(Vmag)
r0 = np.sqrt(xsnew**2+ysnew**2)         # [rcore]
order = np.argsort(r0)
r0 = np.array(r0)[order];   Vmag = np.array(Vmag)[order]
for i in range(len(r0)):
    if np.sum(Vmag[r0<r0[i]])>Vtot/2.:
        rhalfdata = (r0[i]+r0[i-1])/2.  # [rcore]
        break

xsnew = xsnew / rhalfdata * rhalf       # [pc]
ysnew = ysnew / rhalfdata * rhalf       # [pc]

c=open(dir+dwarf+'/centerpos.txt','w')
# print "x y z" on first line, to interprete data later on
print>>c,'x','y','z','vLOS'
c.close()
# print x,y coordinate wrt center of all stars
c=open(dir+dwarf+'/centerpos.txt', 'a')
for k in range(len(xsnew)):
    print>>c,xsnew[k],ysnew[k],vlos[k] # [pc, km/s]
c.close()

from pylab import *
ion(); subplot(111)
en = len(xsnew)
scatter(xsnew[:en], ysnew[:en], c=PM[:en], s=35, vmin=0.95, vmax=1.0, lw=0.0, alpha=0.5)
colorbar()
circ=Circle((0,0), radius=rhalf, fc='None', ec='b', lw=3)
ax=gca()
ax.add_patch(circ)

# visible region
maxr = max(np.abs(xsnew))
mayr = max(np.abs(ysnew))
width2 = max([maxr,mayr])
xlim([-width2,width2])
ylim([-width2,width2])
axes().set_aspect('equal')

xlabel(r'$RA [pc]$'); ylabel(r'$DE [pc]$')
#legend(['\rho','\rho'],dwarf)
title(dwarf)
savefig(dir+dwarf+"/centerpos.png")
ioff(); show()

