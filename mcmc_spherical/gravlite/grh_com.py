#!/usr/bin/python2.7
import numpy as np
import params as ps

print ps.simpos
x,y,z=np.loadtxt(ps.simpos,skiprows=1,unpack=True)
ndm=len(x)

centreofmassx=np.sum(x)/ndm
centreofmassy=np.sum(y)/ndm
centreofmassz=np.sum(z)/ndm

print centreofmassx,centreofmassy,centreofmassz

xnew=x-centreofmassx
ynew=y-centreofmassy
znew=z-centreofmassz

print ps.fileposcenter
filepos=open(ps.fileposcenter,'w')
print>>filepos,'x','y','z'
for k in range(ndm):
    print>>filepos,xnew[k],ynew[k],ynew[k]
filepos.close()
