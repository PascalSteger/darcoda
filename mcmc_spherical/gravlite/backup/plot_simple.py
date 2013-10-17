#!/usr/bin/python2.7
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
from pylab import *
import pdb
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline
from gl_analytic import Mwalkertot

global f, ax1

dir = '/home/psteger/sci/dwarf_data/data_disc_simple/'
#dir = '/home/ast/read/dark/dwarf_data/data_disc_simple/'
nampart = '20130527113603_cprior_nulog_mslope_bprior_rprior' # first, failed after 56k its
nampart = '20130527192624_cprior_nulog_mslope_bprior_rprior' # 27k, stuck with Surf
nampart = '20130527194538_cprior_nulog_mslope_bprior_rprior' # 16k, sent Justin and Dave
nampart = '20130528081522_cprior_nulog_mslope_bprior_rprior' # 100k, overestimation?
nampart = '20130528090813_cprior_nulog_mslope_bprior_rprior_quad' # 50k, quadratic, overestimation?
nampart = '20130528101433_cprior_nulog_mslope_bprior_rprior'      # some k, linear, scaling sqrt(2) for nu
nampart = '20130529141312_cprior_nulog_mslope_bprior_rprior'      # 45k, kappa prior, 10k stars
nampart = '20130531121342_cprior_nulog_mslope_bprior_rprior'      # 15k, kappa prior, 2k stars, overestimate
nampart = '20130531131908_cprior_nulog_mslope_bprior_rprior'      # 500, kappa, 2k stars
basename = dir + nampart + '/' + nampart




def M_anf(r):
    return r**2/(r+1.)**2

def ipol(xin,yin,xout,smooth=1.e-9):
    rbf = Rbf(xin, yin, smooth=smooth)
    return rbf(xout)

def ipollog(xin,yin,xout,smooth=1.e-9):
    rbf = Rbf(xin, np.log10(yin), smooth=smooth)
    return 10**(rbf(xout))

def prepare_plots():
    global f,ax1
    ion()
    f, ax1 = plt.subplots(1,1)
    draw()
    return

def setlabels(ax,xtext,ytext):
    ax.set_xlabel(r'$'+xtext+'$')
    ax.set_ylabel(r'$'+ytext+'$')

def setlims(ax,xlim,ylim):
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

def plot_data(radii,prof,col='blue',lw=2,ls='-'):
    ax1.plot(radii, prof, c=col, lw=lw, ls=ls)
    setlabels(ax1,'r [pc]','\Sigma_z [M_\odot/pc^2]')#[10^5 M_{\odot}]')
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_xlim([min(radii),max(radii)])
    #ax1.errorbar(dp.dat.Mx, dp.dat.Mdat, yerr=dp.dat.Merr, c='green')
    #setlabels(ax1,'r [R_s]','beta')
    # ax1.set_xscale('log'); ax1.set_yscale('log')
    
    # setlims(ax1,[0,2*rplmax],[0.1,10.**6])
    draw()

def save_plot(nam):
    plt.savefig(nam)
    return

def show_plots():
    ioff();show()
    return


print 'input'
print basename
M = np.loadtxt(basename+'.profM',skiprows=0,unpack=False)
#M = np.loadtxt(basename+'.profbeta1',skiprows=0,unpack=False)

radii = M[0]
print radii
radii = radii
profs = M[1:]

prepare_plots()
#for i in range(10):# len(profs)):
#    plot_data(radii,profs[i])
#plot_data(radii,profs[len(profs)/2],'yellow')

Mprofbins = np.transpose(profs)

for i in range(len(Mprofbins)):
    # sort all mass models bin by bin
    Mprofbins[i] = np.sort(Mprofbins[i])

bins=len(radii)
Mmedi = np.zeros(bins); Mmax  = np.zeros(bins); Mmin  = np.zeros(bins)
M95hi = np.zeros(bins); M95lo = np.zeros(bins)
M68hi = np.zeros(bins); M68lo = np.zeros(bins)
mlen = len(Mprofbins[0])
for i in range(len(radii)):
    Mmax[i]  = Mprofbins[i,mlen-1]
    M95hi[i] = Mprofbins[i,0.95 * mlen]
    M68hi[i] = Mprofbins[i,0.68 * mlen]
    Mmedi[i] = Mprofbins[i,0.50 * mlen]
    M68lo[i] = Mprofbins[i,0.32 * mlen]
    M95lo[i] = Mprofbins[i,0.05 * mlen]
    Mmin[i]  = Mprofbins[i,0]

ax1.fill_between(radii, M95lo, M95hi, color='k', alpha = 0.5)
ax1.fill_between(radii, M68lo, M68hi, color='k', alpha = 0.3)
setlabels(ax1,'r [pc]','\Sigma_z [M_\odot/pc^2]')#[10^5 M_{\odot}]')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlim([min(radii),max(radii)])

    



def readcol(filena):
    a,b,c = np.loadtxt(filena,skiprows=1,unpack=True)
    return a,b,c

G1  = 6.67398e-11                       # [m^3 kg^-1 s^-2]
pc  = 3.08567758e16                     # [m]
msun= 1.981e30                          # [kg]
km  = 1000.                             # [m]
G1  = G1*msun/km**2/pc                  # [pc msun^-1 km^2 s^-2]

K = 1650.
D = 250.
F = 0.1*1.650

Mdat = (K*radii/np.sqrt(radii**2.+D**2.)+2.*F*radii) / (2.0*np.pi*G1) / 1000.

ax1.plot(radii,Mdat,'k',lw=2)

ax1.plot(radii,Mmedi,'r',ls='--',lw=2)


fout = open(basename+".profM.conf",'w')
print >> fout,M95lo
print >> fout,M68lo
print >> fout,Mmedi
print >> fout,M68hi
print >> fout,M95hi
fout.close()

analyt = M_anf(radii)

print '# radii  lower 95%    lower 68%   median      upper 68%   upper 95%   analytic'
for i in range(len(radii)):
    print radii[i],M95lo[i],M68lo[i],Mmedi[i],M68hi[i],M95hi[i],analyt[i]

save_plot(basename+".profM.png")
show_plots()

