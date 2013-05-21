#!/usr/bin/python
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
# move centered positions to spherical coordinates
import numpy as np
import params as ps
import pdb
import random

print 'input:'
print ps.fileposcenter
xall,yall,zall    = np.loadtxt(ps.fileposcenter, skiprows=1,unpack=True)
print ps.simvel
vxall,vyall,vzall = np.loadtxt(ps.simvel,        skiprows=1,unpack=True)

rall = np.sqrt(xall**2+yall**2+zall**2)
rest = ( rall < 2.*ps.rmax )
xall = xall[rest];  yall = yall[rest];  zall = zall[rest]
vxall= vxall[rest]; vyall= vyall[rest]; vzall= vzall[rest]

#  get ntracer random points
if ps.ntracers1>0:
  n = min(ps.ntracers1,len(xall)-1)
  trace = random.sample(xrange(len(xall)), n)
else:
  n = min(ps.ntracers2,len(xall)-1)
  trace = random.sample(xrange(len(xall)), n)

x  = [ xall[i]    for i in trace ]
y  = [ yall[i]    for i in trace ]
z  = [ zall[i]    for i in trace ]

print 'output:'
print ps.fileposcartesian
fileposcartesian=open(ps.fileposcartesian,'w')
print>>fileposcartesian,'x','y','z'
for k in range(n):
  print>>fileposcartesian,x[k],y[k],z[k]
fileposcartesian.close()

vx = [ vxall[i]   for i in trace ]
vy = [ vyall[i]   for i in trace ]
vz = [ vzall[i]   for i in trace ]

print ps.filevelcartesian
filevelcartesian=open(ps.filevelcartesian,'w')
print>>filevelcartesian,'x','y','z'
for k in range(n):
  print>>filevelcartesian,vx[k],vy[k],vz[k]
filevelcartesian.close()

# r
x = np.array(x);  y = np.array(y);  z=np.array(z)
r = np.sqrt(x**2+y**2+z**2)

frac = np.sum((r<ps.rmax))*1./len(r)
print 'fraction of particles nearer than rmax: ',frac*100,'%'

# theta (polar angle [0,pi]); phi (azimuthal angle [-pi,pi])
theta = np.arccos(z/r)
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

print ps.fileposspherical
fileposspherical=open(ps.fileposspherical,'w')
print>>fileposspherical,'r','phi (azimuthal)','theta (polar)'
for k in range(n):
  print>>fileposspherical,r[k],phi[k],theta[k]
fileposspherical.close()


# Velocities #
vx = np.array(vx);  vy = np.array(vy);  vz = np.array(vz)
vr = (x*vx+y*vy+z*vz)/np.sqrt(x**2+y**2+z**2)

thetadot = -(r*vz-z*vr)/np.sqrt(r**4-z**2*r**2)
vtheta   = r*thetadot

phidot = (x*vy-y*vx)/(x**2+y**2)
vphi   = r*phidot*np.sin(theta)

print ps.filevelspherical
filevelspherical=open(ps.filevelspherical,'w')
print>>filevelspherical,'vr','vphi (azimuthal)','vtheta (polar)'
for k in range(n):
  print>>filevelspherical,vr[k],vphi[k],vtheta[k]
filevelspherical.close()
