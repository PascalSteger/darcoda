#!/usr/bin/python2.7
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
'''calculate line of sight velocity dispersion (line of sight parallel to x-axis)'''

import params as ps
import numpy as np
import multiprocessing as mp
import pdb
from BiWeight import meanbiweight

print 'input:'
print ps.fileposcartesian
x,y,z    = np.loadtxt(ps.fileposcartesian, unpack=True, skiprows=1)
print ps.filevelcartesian
vx,vy,vz = np.loadtxt(ps.filevelcartesian, unpack=True, skiprows=1)

cut = np.sqrt(x**2+y**2+z**2)
ind = np.argwhere(cut <= ps.rcut).flatten() #cut everything outside of rcut

r   = np.sqrt(y**2+z**2) #perpendicular distance to x-axis

rs  = r #ps.rerror*np.random.randn(len(r))+r
vxs = vx # ps.vrerror*np.random.randn(len(vx))+vx

rs   = rs[ind]
vxs  = vxs[ind]

binmin,binmax,rbin = ps.binparams()

def foo_pool(k):
  rsi  = ps.rerror*np.random.randn(len(rs))+rs
  vxsi = ps.vrerror*np.random.randn(len(vxs))+vxs
  locdisp=[]; loca = []
  for i in range(ps.bins):
    ind1 = np.argwhere( np.logical_and(rsi > binmin[i], rsi < binmax[i])).flatten()
    vx1  = vxsi[ind1]
    locdisp.append(meanbiweight(vx1,ci_perc=68.4,ci_mean=True,ci_std=True)[1])
    loca.append(len(ind1))
  return locdisp,loca

dispvel=[]
alog=[]
def log_result(result):
  # This is called whenever foo_pool(i) returns a result.
  # result_list is modified only by the main process, not the pool workers.
  d,a=result
  dispvel.append(d)
  alog.append(a)

pool = mp.Pool(processes=ps.procs)
for k in range(ps.nit):
  pool.apply_async(foo_pool, args = (k, ), callback = log_result)

pool.close()
pool.join()

sigarr = np.array(dispvel)
abarr  = np.array(alog)

print 'output:'
print ps.filesig
filesig=open(ps.filesig,'w')
print>>filesig,'r','sigma_los(r)','error'
for i in range(ps.bins):
  sigma = np.sum(sigarr[:,i])/ps.nit
  ab    = np.sum(abarr[:,i])/ps.nit
  dispvelerror = sigma/np.sqrt(ab)
  print>>filesig,rbin[i],sigma,dispvelerror
  print rbin[i],sigma,dispvelerror
filesig.close()


