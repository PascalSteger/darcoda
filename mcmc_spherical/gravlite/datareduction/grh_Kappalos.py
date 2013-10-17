#!/usr/bin/env ipython-python3.2
# calculate line of sight velocity kurtosis (line of sight parallel to z-axis)
import numpy as np
import multiprocessing as mp
import pdb

import gl_params as gp
import gr_params as gpr
from gl_helper import expDtofloat, bin_r_linear, bin_r_log, bin_r_const_tracers
from scipy.stats import kurtosis

def run():
    for comp in range(gpr.ncomp):
        print('input:', gpr.fileposspherical[comp])
        R,Phi,vlos = np.loadtxt(gpr.fileposspherical[comp],\
                                comments='#', unpack=True)

        Rs  = R[:] # gpr.Rerror*np.random.randn(len(r))+r # for initial offset
        vloss = vlos[:]  # gpr.vrerror*np.random.randn(len(vx))+vx # same same

        Rmin = min(Rs); Rmax = max(Rs)
        if gp.lograd:
            print(gpr.nbins,' bins in log spacings')
            Binmin, Binmax, Rbin = bin_r_log(Rmax/gpr.nbins, Rmax, gpr.nbins)
        elif gp.consttr:
            print(len(R)/gpr.nbins,' particles per bin')
            Binmin, Binmax, Rbin = bin_r_const_tracers(Rs, len(Rs)/gpr.nbins)
        else:
            print(gpr.nbins, ' bins in linear spacings')
            Binmin, Binmax, Rbin = bin_r_linear(Rmin, Rmax, gpr.nbins)

        # if Dispvel is [] still after pool call,
        # some error occured inside following function:
        def foo_pool(k):
            Rsi  = gpr.Rerror  * np.random.randn(len(Rs))+Rs
            vlossi = gpr.vrerror * np.random.randn(len(vloss))+vloss
            lockappa=[]; loca = []
            for b in range(gpr.bins):
                ind1 = np.argwhere( np.logical_and(Rsi > Binmin[b],\
                                                   Rsi <= Binmax[b])).flatten()
                vlos1  = vlossi[ind1]
                lockappa.append(kurtosis(vlossi, axis=0,\
                                         fisher=False, bias=False))
                loca.append(len(ind1))
            return lockappa,loca

        Kappavel=[]; alog=[]
        def log_result(result):
            # This is called whenever foo_pool(i) returns a result.
            # result_list is modified only by main process, not pool workers.
            kappa, alo = result
            Kappavel.append(kappa)
            alog.append(alo)

        pool = mp.Pool(processes=gpr.procs)
        for k in range(gpr.nit):
            pool.apply_async(foo_pool, args = (k, ), callback = log_result)
        pool.close()
        pool.join()

        Kappaarr = np.array(Kappavel)
        abarr  = np.array(alog)

        print('output:', gpr.filekappa[comp])
        filekappa = open(gpr.filekappa[comp],'w')
        print('# R [pc]','Binmin [pc]','Binmax [pc]',\
              'Kappa_los(R) [1]','error [1]', file=filekappa)

        for b in range(gpr.bins):
            Kappa = np.sum(Kappaarr[:,b])/gpr.nit
            ab    = np.sum(abarr[:,b])/gpr.nit
            Kappavelerror = Kappa/np.sqrt(ab)
            print(Rbin[b], Binmin[b], Binmax[b], \
                  Kappa, Kappavelerror, file=filekappa)
            print(Rbin[b], Kappa, Kappavelerror)
        filekappa.close()
    return

if __name__=="__main__":
    run()
