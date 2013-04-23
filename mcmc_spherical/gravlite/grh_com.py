#!/usr/bin/python2.7
import numpy as np
import gl_params as gp
import gr_params as gpr

print 'input: ', gpr.simpos
x,y,z=np.loadtxt(gpr.simpos,skiprows=1,unpack=True)
ndm=len(x)

centreofmassx=np.sum(x)/ndm
centreofmassy=np.sum(y)/ndm
centreofmassz=np.sum(z)/ndm

print centreofmassx,centreofmassy,centreofmassz

xnew = (x-centreofmassx)*gp.ascale
ynew = (y-centreofmassy)*gp.ascale
znew = (z-centreofmassz)*gp.ascale

print 'output: ', gpr.fileposcenter
filepos=open(gpr.fileposcenter,'w')
print>>filepos,'x','y','z'
for k in range(ndm):
    print>>filepos,xnew[k],ynew[k],ynew[k]
filepos.close()
