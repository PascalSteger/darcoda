#!/usr/bin/python
# move centered positions to spherical coordinates
import numpy as np
import gl_params as gp
import gr_params as gpr
import pdb
import random

print 'input:'
print gpr.fileposcenter
xall,yall,zall    = np.loadtxt(gpr.fileposcenter, skiprows=1,unpack=True)
print gpr.simvel
vxall,vyall,vzall = np.loadtxt(gpr.simvel,        skiprows=1,unpack=True)

rall = np.sqrt(xall**2+yall**2)
rest = ( rall < gpr.rmax )
xall = xall[rest];  yall = yall[rest];  zall = zall[rest]
vxall= vxall[rest]; vyall= vyall[rest]; vzall= vzall[rest]

#  get ntracer random points
if gpr.ntracers1>0:
  n = min(gpr.ntracers1,len(xall)-1)
  trace = random.sample(xrange(len(xall)), n)
else:
  n = min(gpr.ntracers2,len(xall)-1)
  trace = random.sample(xrange(len(xall)), n)

x  = [ xall[i]    for i in trace ]
y  = [ yall[i]    for i in trace ]

print 'output:'
print gpr.fileposcartesian
fileposcartesian=open(gpr.fileposcartesian,'w')
print>>fileposcartesian,'x','y'
for k in range(n):
  print>>fileposcartesian,x[k],y[k]
fileposcartesian.close()

vz = [ vzall[i]   for i in trace ]

print gpr.filevelcartesian
filevelcartesian=open(gpr.filevelcartesian,'w')
print>>filevelcartesian,'vz'
for k in range(n):
  print>>filevelcartesian,vz[k]
filevelcartesian.close()

# r
x = np.array(x);  y = np.array(y)
r = np.sqrt(x**2+y**2)

frac = np.sum((r<gpr.rmax))*1./len(r)
print 'fraction of particles nearer than rmax: ',frac*100,'%'

# phi (azimuthal angle [-pi,pi])
phi   = np.zeros(n)
for k in range(n):
  if x[k]>0:
    phi[k] = np.arctan(y[k]/x[k])
  elif x[k]==0:
    phi[k] = np.sign(y[k])*np.pi/2
  elif y[k]>=0:
    phi[k] = np.arctan(y[k]/x[k])+np.pi
  elif y[k]<0:
    phi[k] = np.arctan(y[k]/x[k])-np.pi

print gpr.fileposspherical
fileposspherical=open(gpr.fileposspherical,'w')
print>>fileposspherical,'r','phi'
for k in range(n):
  print>>fileposspherical,r[k],phi[k]
fileposspherical.close()

# no more filevelspherical, only LOS needed
