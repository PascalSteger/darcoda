#!/usr/bin/python2.7
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import numpy
# calculate approximative stellar masses from V-band magnitude, V-I color, distance to dwarf galaxy
# read in [some].bin.[MV,V-I], output ID,Mstar

# choose simulation

import sys
if(len(sys.argv)<2):
    print "use: stellarmass.py [car,scl,sex,for]"
    exit(1)
    
dwarf=sys.argv[1]
dir="/home/ast/read/dark/dwarf_data/"
print dir+dwarf+"/table_merged.bin"

delim=[0,22,3,3,6,4,3,5,6,6,7,5,6,5,6,5,6]
ID=numpy.genfromtxt(dir+dwarf+"/table_merged.bin",skiprows=29,unpack=True,usecols=(0,1),delimiter=delim,dtype="string")
RAh,RAm,RAs,DEd,DEm,DEs,Vmag,VI,VHel,e_VHel,SigFe,e_SigFe,SigMg,e_SigMg,PM=numpy.genfromtxt(dir+dwarf+"/table_merged.bin",skiprows=29,unpack=True,usecols=tuple(range(2,17)),delimiter=delim,filling_values=-1)

print 'Vmag = ',Vmag[0:10]

MCMD,MvCMD,VICMD=numpy.genfromtxt(dir+dwarf+"/../SCMD/SCMD.dat",skiprows=1,unpack=True,usecols=(0,4,8),filling_values=-1)

# luminosity distance, measured in [parsecs]
kpc = 1000
DL= {
    'for': lambda x: x * (138),#+/- 8
    'car': lambda x: x * (101),#+/- 5
    'sex': lambda x: x * (86), #+/- 4
    'scl': lambda x: x * (79)  #+/- 4
    }[dwarf](kpc)

# print "DL = ",DL

from pylab import *
ion();subplot(111)
n,bins,rectangles = hist(PM, bins=20, normed=True)
axvline(x=0.95,color='r')
xlabel(r'PM')
ylabel(r'number')
xlim([0,1])
savefig(dir+dwarf+"/PM.png")
ioff();clf()

# only use stars which are members of the dwarf
pm = (PM>0.95)*(VI<70)
print pm
print "fraction of members = ",1.0*sum(pm)/len(pm)
ID=ID[1][pm]; RAh=RAh[pm]; RAm=RAm[pm]; DEd=DEd[pm]; DEm=DEm[pm]; DEs=DEs[pm]; Vmag = Vmag[pm]; VI=VI[pm]; VHel=VHel[pm]; e_VHel=e_VHel[pm]; SigFe=SigFe[pm]; e_SigFe=e_SigFe[pm]; SigMg=SigMg[pm]; e_SigMg=e_SigMg[pm];PM=PM[pm]


VMag = Vmag-5.0*(numpy.log10(DL)-1.0)
minVMag,maxVMag=numpy.min(VMag),numpy.max(VMag)
minVI,maxVI    =numpy.min(VI),numpy.max(VI)
print "min, max of VMag= ",minVMag,maxVMag

windowCMD = (minVMag<=MvCMD)*(MvCMD<=maxVMag)*(minVI<=VICMD)*(VICMD<=maxVI)
MCMD = MCMD[windowCMD]
MvCMD = MvCMD[windowCMD]
VICMD = VICMD[windowCMD]



ion(); subplot(111)
# set_xaxis('log')
plot(VI,VMag,'b.',linewidth=1)
# errorbar(rmean,rho,xerr=rspan,yerr=err,linewidth=3)
plot(VICMD,MvCMD,'r+',linewidth=3)
xyset = [[str(MCMD[i]),VICMD[i],MvCMD[i]] for i in range(len(VICMD))]
print xyset[0:3]
for label, x, y in xyset:
    annotate(label,xy = (x, y))


# visible region
# plt.xlim([10**0,3*10**1])
# plt.ylim([10**-2,10**2])
xlabel(r'$V-I$')
ylabel(r'$M_V$')
# legend(['\rho','\rho'],'lower left')
# title('z=11.7')
ioff();savefig(dir+dwarf+"/HRD.png")
show();clf()
