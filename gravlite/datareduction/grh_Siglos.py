#!/usr/bin/env python3

##
# @file
# calculate line of sight velocity dispersion (line of sight parallel to x-axis)

# (c) 2013 ETHZ, Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import multiprocessing as mp
import pdb
from pylab import *

import gl_params
gp = gl_params.Params()
import gr_params as gpr
from gl_helper import expDtofloat, bin_r_linear, bin_r_log, bin_r_const_tracers
from BiWeight import meanbiweight

def run():
    for comp in range(gpr.ncomp):
        print('input:')
        print(gpr.fileposspherical[comp])
        R,Phi,vz = np.loadtxt(gpr.fileposspherical[comp], unpack=True, skiprows=1)
        vz /= 0.482126                          # [km/s], for G = L = M = 1

        Rs  = R[:] # gpr.Rerror*np.random.randn(len(r))+r # for initial offset
        vzs = vz[:]      # gpr.vrerror*np.random.randn(len(vx))+vx # same same

        Rmin = min(Rs); Rmax = max(Rs)
        if gp.lograd:
            Binmin, Binmax, Rbin = bin_r_log(Rmax/gp.nipol, Rmax, gp.nipol)
            print(gp.nipol,' bins in log spacings')
        elif gp.consttr:
            Binmin, Binmax, Rbin = bin_r_const_tracers(Rs, len(Rs)/gp.nipol)
            print(len(R)/gp.nipol,' particles per bin')
        else:
            Binmin, Binmax, Rbin = bin_r_linear(Rmin, Rmax, gp.nipol)
            print(gp.nipol, ' bins in linear spacings')

        # if Dispvel is [] still after pool call,
        # some error occured inside following function:
        def foo_pool(k):
            Rsi  = gpr.Rerror  * np.random.randn(len(Rs))+Rs
            vzsi = gpr.vrerror * np.random.randn(len(vzs))+vzs
            locdisp=[]; loca = []
            for b in range(gpr.bins):
                ind1 = np.argwhere( np.logical_and(Rsi > Binmin[b],\
                                                   Rsi <= Binmax[b])).flatten()
                vz1  = vzsi[ind1]
                locdisp.append(meanbiweight(vz1,\
                                            ci_perc=68.4,ci_mean=True,\
                                            ci_std=True)[1])
                loca.append(len(ind1))
            return locdisp,loca

        Dispvel=[]; alog=[]
        def log_result(result):
            # This is called whenever foo_pool(i) returns a result.
            # result_list is modified only by the main process, not the pool workers
            dis, alo = result
            Dispvel.append(dis)
            alog.append(alo)

        pool = mp.Pool(processes=gpr.procs)
        for k in range(gpr.nit):
            pool.apply_async(foo_pool, args = (k, ), callback = log_result)
        pool.close()
        pool.join()

        Sigarr = np.array(Dispvel)
        abarr  = np.array(alog)

        filesig = open(gp.files.sigfiles[comp],'w')
        print('R [pc]','Binmin [pc]','Binmax [pc]',\
              'Sigma_los(R) [km/s]','error [km/s]', file=filesig)

        P_Sigma = np.zeros(gpr.bins); P_ESigma = np.zeros(gpr.bins)
        for b in range(gpr.bins):
            P_Sigma[b] = np.sum(Sigarr[:,b])/gpr.nit
            ab    = np.sum(abarr[:,b])/gpr.nit
            P_ESigma[b] = P_Sigma[b]/np.sqrt(ab)
            print(Rbin[b],Binmin[b],Binmax[b],P_Sigma[b],P_ESigma[b], file=filesig)
            print(Rbin[b],P_Sigma[b],P_ESigma[b])
        filesig.close()

        if not gpr.showplots: continue
        ion(); subplot(111)
        plot(Rbin,P_Sigma,'b',lw=1)
        fill_between(Rbin,P_Sigma-P_ESigma,P_Sigma+P_ESigma,alpha=0.5,color='r')
        # [rscale],2*[km/s]

        xlabel(r'$R [\mathrm{Rscale}]$')
        ylabel(r'$\langle\sigma_{\mathrm{LOS}}\rangle [\mathrm{km/s}]$')
        # ylim([-5,30])
        # xlim([0,3])
        savefig(gpr.get_siglos_png(comp))
        if gpr.showplots:
            ioff();show();clf()

    return

if __name__=="__main__":
    run()
