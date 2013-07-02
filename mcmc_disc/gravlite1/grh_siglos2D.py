#!/usr/bin/python
# calculate line of sight velocity dispersion (line of sight parallel to x-axis)
import gl_params as gp
import gr_params as gpr
import gl_helper as gh
import numpy as np
import multiprocessing as mp
import pdb
from BiWeight import meanbiweight

print 'input:'
print gpr.fileposcartesian
x,y    = np.loadtxt(gpr.fileposcartesian, unpack=True, skiprows=1)
print gpr.filevelcartesian
vz = np.loadtxt(gpr.filevelcartesian, unpack=True, skiprows=1)
vz /= 0.482126                          # [km/s], for G = L = M = 1

cut = np.sqrt(x**2+y**2)
ind = np.argwhere(cut <= gpr.rcut).flatten() #cut everything outside 2D radius rcut

r   = np.sqrt(x**2+y**2) #perpendicular distance to x-axis

rs  = r #gpr.rerror*np.random.randn(len(r))+r
vzs = vz # gpr.vrerror*np.random.randn(len(vx))+vx

rs   = rs[ind]
vzs  = vzs[ind]

if gp.lograd:
  binmin,binmax,rbin = gh.bin_r_log(gpr.rmax/gpr.nbins,gpr.rmax,gp.nbins)
else:
  binmin,binmax,rbin = gh.bin_r_linear(gpr.rmin,gpr.rmax,gp.nbins)



def foo_pool(k):
  rsi  = gpr.rerror*np.random.randn(len(rs))+rs
  vzsi = gpr.vrerror*np.random.randn(len(vzs))+vzs
  locdisp=[]; loca = []
  for i in range(gpr.bins):
    ind1 = np.argwhere( np.logical_and(rsi > binmin[i], rsi < binmax[i])).flatten()
    vz1  = vzsi[ind1]
    locdisp.append(meanbiweight(vz1,ci_perc=68.4,ci_mean=True,ci_std=True)[1])
    loca.append(len(ind1))
  return locdisp,loca

dispvel=[]
alog=[]
def log_result(result):
  # This is called whenever foo_pool(i) returns a result.
  # result_list is modified only by the main process, not the pool workers.
  d, a = result
  dispvel.append(d)
  alog.append(a)

pool = mp.Pool(processes=gpr.procs)
for k in range(gpr.nit):
  pool.apply_async(foo_pool, args = (k, ), callback = log_result)

pool.close()
pool.join()

sigarr = np.array(dispvel)
abarr  = np.array(alog)

print 'output:'
print gpr.filesig
filesig=open(gpr.filesig,'w')
print>>filesig,'r','sigma_los(r)','error'
for i in range(gpr.bins):
  sigma = np.sum(sigarr[:,i])/gpr.nit
  ab    = np.sum(abarr[:,i])/gpr.nit
  dispvelerror = sigma/np.sqrt(ab)
  print>>filesig,rbin[i],sigma,dispvelerror
  print rbin[i],sigma,dispvelerror
filesig.close()


