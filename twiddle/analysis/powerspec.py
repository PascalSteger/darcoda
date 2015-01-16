#!/usr/bin/env python3

## \file
# plot scatter plots from halos file

import sys
i=len(sys.argv)
if i<4:
    print("Usage: scatter.py infile outfile columnx(start from 1) columny [log]")
    print("Example: scatter.py halos scatter.png 12 4 1")
    exit(1)
    
infile  = sys.argv[1]; outfile = sys.argv[2]
colx     = [int(sys.argv[3])-1]
coly     = [int(sys.argv[4])-1]
cole     = [6-1]
boolog  = False
if(i>5):
    boolog  = int(sys.argv[5])>0

from matplotlib import pyplot as PLT
fig = PLT.figure()
ax1 = fig.add_subplot(111)

import numpy as NP
with open(infile) as f:
    vx = NP.loadtxt(f, dtype='float', comments="#", skiprows=0, usecols=colx)
    vx = NP.ravel(vx)   # 'flatten' v, cumulative sum
with open(infile) as g:
    vy = NP.loadtxt(g, dtype='float', comments="#", skiprows=0, usecols=coly)
    vy = NP.ravel(vy)   # 'flatten' v, cumulative sum
with open(infile) as h:
    ey = NP.loadtxt(h, dtype='float', comments="#", skiprows=0, usecols=cole)
    ey = NP.ravel(ey)   # 'flatten' v, cumulative sum

vxhi= vx;  vyhi= vy/(10**6/0.702)**3
eyhi = ey*vyhi

if(boolog):
    vxhi = NP.log10(vxhi);    vyhi = NP.log10(vyhi)

import math; from math import *
k = 1./((1./0.702)/(vxhi+1)); pk = vyhi
kfine=NP.logspace(0,NP.log10(max(vxhi)),1000)
Ob = 0.0455; Oc = 0.728; Om = 0.272; h=0.702
keq= 7.46/100.*Om

zeq= 3233. #http://en.wikipedia.org/wiki/Wilkinson_Microwave_Anisotropy_Probe
b1 = 0.313*(Om*h**2)**(-0.419)*(1.+0.607*(Om*h**2)**0.674)
b2 = 0.238*(Om*h**2)**0.223
zd = 1291.*((Om*h**2)**0.251/(1.+0.659*(Om*h**2)**0.828))*(1.+b1*(Ob*h**2)**b2) 
Req= 31.5*Ob/(zeq/1000.)*h**2
Rd = 31.5*Ob/(zd/1000.)*h**2
s  = 2./(3*keq)*math.sqrt(6./Req)*NP.log((sqrt(1+Rd)+sqrt(Rd+Req))/(1.+sqrt(Req)))
#f  = 1./(1.+(kfine*s/5.4)**4)

a1 = (46.9*Om*h**2)**0.670*(1.+(32.1*Om*h**2)**-0.532)
a2 = (12.0*Om*h**2)**0.424*(1.+(45.0*Om*h**2)**-0.582)
ac = a1**(-Ob/Om)*a2**(-(Ob/Om)**3)
b1 = 0.944*(1.+(458*Om*h**2)**-0.708)**-1.0
b2 = (0.395*Om*h**2)**-0.0266
bc = 1/(1.+b1*((Oc/Om)**b2-1.0))
q  = kfine/(13.41*keq)
Tck= ac * NP.log(1.8*bc*q)/(14.2*q**2)

y  = (1.+zeq)/(1.+zd)
Gy = y*(-6.*sqrt(1.+y)+(2.+3.*y)*NP.log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1)))
ab = ab=2.07*keq*s*(1.+Rd)**-0.75*Gy
Dk = 1.0*NP.exp(-kfine/20)
Tbk= ab*NP.sin(kfine*s)/(kfine*s)*Dk
Tk = Oc/Om*Tck+Ob/Om*Tbk


#Peacock, BBKS
O  = Oc
q  = kfine/(O*h*NP.exp(-Ob*(1.+sqrt(2.*h)/O)))/h
Tkadi= NP.log(1.+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)**2+(5.46*q)**3+(6.71*q)**4)**(-0.25)
Tkiso= (5.6*q)**2*(1.+(15.*q+(0.9*q)**(3./2.)+(5.6*q)**2)**1.24)**(-1/1.24)

Or = 8.24E-5
z  = 1000.
Oz = Om*(1.+z)**3/(Oc+Or*(1+z)**2+Om*(1+z)**3)
OLz= Oc/(Oc+Or*(1+z)**2+Om*(1+z)**3)
D1z= 1./(1.+z)*(5.*Oz)/2*(Oz**(4./7.)-OLz+(1.+Oz/2.)*(1.+OLz/70.))**(-1.)

delh=1.94*10**(-5)*Om**(-0.785-0.05*NP.log(Om))
cH0 =4270.55
fac =delh**2*(cH0)**4*2*math.pi**2
#fac = (10.**6/(h*2.*math.pi))**3

pkbbks= Tkadi**2*fac*D1z
pkehu = Tk**2*fac*D1z

PLT.errorbar(k,pk,eyhi,fmt='o',alpha=0.5,label='IC by grafic')
PLT.plot(kfine,pkbbks,alpha=0.5,label='BBKS adiabatic CDM')
PLT.plot(kfine,pkehu,alpha=0.5,label='Eisenstein,Hu 1997')
PLT.xscale('log'); PLT.yscale('log')
PLT.xlim(0.5,512)
PLT.ylim(10**-8,1000)
PLT.xlabel(r'k [h/Mpc]')
PLT.xticks(NP.logspace(0,2,3),['1','10','100'])
PLT.ylabel(r'P(k) [(Mpc/h)$^3$]')
PLT.yticks(NP.logspace(-7,3,6))
PLT.legend(loc=1)
PLT.savefig(outfile)

