#!/usr/bin/python2.7
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import params as ps
import numpy as np
import pdb
import multiprocessing as mp

# calculate density falloff
binmin,binmax,rbin = ps.binparams()

# volume of a bin with height binlength
vol=np.zeros(ps.bins)
for i in range(ps.bins):
  vol[i]=4./3.*np.pi*(binmax[i]**3-binmin[i]**3)

print 'input:'
print ps.fileposspherical
r,phi,theta = np.loadtxt(ps.fileposspherical, unpack=True, skiprows=1)
# rs=rerror*np.random.randn(len(r))+r
# rs not changed later on, ever: misplacement, for all realizations. wanted?
# yes, we scatter radii in foo_pool
rs = r

def foo_pool(k):
  rsi = ps.rerror*np.random.randn(len(rs))+rs
  locdens = []; loca = []
  for i in range(ps.bins):
    ind1 = np.argwhere( np.logical_and(rsi > binmin[i], rsi < binmax[i])).flatten()
    locdens.append(1.*len(ind1)/vol[i]*ps.munit)
    loca.append(1.*len(ind1))
  return locdens,loca

density=[]; alog=[]
def log_result(result):
  # This is called whenever foo_pool(i) returns a result.
  # density/alog are modified only by the main process, not the pool workers.
  dtmp,atmp=result
  density.append(dtmp)
  alog.append(atmp)

pool = mp.Pool(processes=ps.procs)
for k in range(ps.nit):
  pool.apply_async(foo_pool, args = (k, ), callback = log_result)
pool.close()
pool.join()

denarr = np.array(density)
dens0  = np.sum(denarr[:,0])#/ps.nit

aarr = np.array(alog)
ab0  = np.sum(aarr[:,0])/ps.nit
denserr0 = dens0/np.sqrt(ab0)

print 'output:'
print ps.filedenfalloff
file=open(ps.filedenfalloff,'w')
print>>file,'r','nu(r)','error'
for i in range(ps.bins):
  dens=np.sum(denarr[:,i])#/ps.nit
  ab=np.sum(aarr[:,i])/ps.nit
  if i==0:
    denserror=0
    denserr=dens/np.sqrt(ab)
  else:
    denserr=dens/np.sqrt(ab)
    denserror=np.sqrt((denserr/dens0)**2+(dens/(dens0**2)*denserr0)**2) 
  # print>>file,rbin[i],dens/dens0,denserror
  print>>file,rbin[i],dens,denserr
  print rbin[i],dens,denserr
file.close()


