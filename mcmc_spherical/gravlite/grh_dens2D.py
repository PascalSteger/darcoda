#!/usr/bin/python
import gl_params as gp
import gr_params as gpr
import numpy as np
import pdb
import multiprocessing as mp

# calculate density falloff
binmin,binmax,rbin = gpr.binparams()

# volume of a bin with height binlength, 2D
vol=np.zeros(gpr.bins)
for i in range(gpr.bins):
  vol[i] = np.pi*(binmax[i]**2-binmin[i]**2)

print 'input:'
print gpr.fileposspherical
r,phi = np.loadtxt(gpr.fileposspherical, unpack=True, skiprows=1)
# rs=rerror*np.random.randn(len(r))+r
# rs not changed later on, ever: misplacement, for all realizations. wanted?
# yes, we scatter radii in foo_pool
rs = r

def foo_pool(k):
  rsi = gpr.rerror*np.random.randn(len(rs))+rs
  locdens = []; loca = []
  for i in range(gpr.bins):
    ind1 = np.argwhere( np.logical_and(rsi > binmin[i], rsi < binmax[i])).flatten()
    locdens.append(1.*len(ind1)/vol[i]*gpr.munit)
    loca.append(1.*len(ind1))
  return locdens,loca

density=[]; alog=[]
def log_result(result):
  # This is called whenever foo_pool(i) returns a result.
  # density/alog are modified only by the main process, not the pool workers.
  dtmp, atmp = result
  density.append(dtmp)
  alog.append(atmp)

pool = mp.Pool(processes=gpr.procs)
for k in range(gpr.nit):
  pool.apply_async(foo_pool, args = (k, ), callback = log_result)
pool.close()
pool.join()

denarr = np.array(density)
dens0  = np.sum(denarr[:,0])#/gpr.nit

aarr = np.array(alog)
ab0  = np.sum(aarr[:,0])/gpr.nit
denserr0 = dens0/np.sqrt(ab0)

print 'output:'
print gpr.filedenfalloff
file = open(gpr.filedenfalloff,'w')
print >> file,'r','nu(r)','error'
for i in range(gpr.bins):
  dens=np.sum(denarr[:,i])#/gpr.nit
  ab=np.sum(aarr[:,i])/gpr.nit
  if i==0:
    denserror = 0
    denserr   = dens/np.sqrt(ab)
  else:
    denserr   = dens/np.sqrt(ab)
    denserror = np.sqrt((denserr/dens0)**2+(dens/(dens0**2)*denserr0)**2) 
  # print >> file,rbin[i],dens/dens0,denserror
  print >> file,rbin[i],dens,denserr
  print rbin[i],dens,denserr
file.close()
