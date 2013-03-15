#!/usr/bin/python2.7
# calculate density falloff

import numpy as np
import multiprocessing as mp
import params as ps

# set binning
bins = 2*ps.bins
rmin =   ps.rmin
rmax = 2*ps.rmax
binlength=(rmax-rmin)/bins
rbin=np.zeros(bins+1)
for i in range(1,bins+1):
  rbin[i]=rmin+(i-1.)*binlength+0.5*binlength

print 'input:'
print ps.fileposspherical
r,phi,theta = np.loadtxt(ps.fileposspherical,unpack=True,skiprows=1)
ndm = ps.ntot;#len(r)
rs  = r       #ps.rerror*np.random.randn(ndm)+r

def foo_pool(k):
  rsi = ps.rerror*np.random.randn(len(rs))+rs
  locmass = []; loca = []
  for i in range(bins):
    ind1=np.argwhere(np.logical_and( rsi > rbin[i], rsi < rbin[i+1])).flatten()
    locmass.append(1.*len(ind1)/ndm)
    loca.append(1.*len(ind1))
  return locmass,loca

mass=[]
alog=[]
def log_result(result):
  m,a=result
  mass.append(m)
  alog.append(a)

pool=mp.Pool(processes=24)
for k in range(ps.nit):
  pool.apply_async(foo_pool, args=(k, ), callback=log_result)
pool.close()
pool.join()

massarr = np.array(mass)
aarr = np.array(alog)

print 'output:'
print ps.filemass
filemass = open(ps.filemass,'w')
print>>filemass,'r','M(r)','error'
totmass=0
ac=0
for i in range(bins):
  totmass+=np.sum(massarr[:,i])/ps.nit
  ac+=np.sum(aarr[:,i])/ps.nit
  masserror=totmass/np.sqrt(ac)
  print>>filemass,rbin[i+1],totmass,masserror
  print rbin[i+1],totmass,masserror
filemass.close()


