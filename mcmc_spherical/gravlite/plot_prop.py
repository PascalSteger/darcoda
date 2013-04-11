#!/usr/bin/python2.7
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
from pylab import *
import pdb
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline
from gl_analytic import Mwalkertot

global f, ax1

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

def plot_data(radii,prof,col='blue',lw=2):
    ax1.plot(radii, prof, c=col, lw=lw)
    #ax1.errorbar(dp.dat.Mx, dp.dat.Mdat, yerr=dp.dat.Merr, c='green')
    setlabels(ax1,'r [pc]','M [M_\odot]')#[10^5 M_{\odot}]')
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

# sample data, not used anymore
# dir = "/home/ast/read/dark/dwarf_data/programs/"
# basename = dir+"/my.profM"

#dir = '/home/ast/read/dark/dwarf_data/hernquist_justin/'
#'10000_no_testplot_20000_models_run_1/'
#basename = dir + "/_pop0_cprior_beta_mslope_10000_0_24.profM"
#basename = dir + '/_pop0_cprior_beta_mslope_1000_0_12.profM' # working 1000 with beta=0
#basename = dir +'_pop0_cprior_mslope_1000_0_12.profM' # working 1000 with beta allowed to change slightly
#basename = dir + '_pop0_cprior_beta_mslope_10000_0_16.profM'
#basename = dir + '_pop0_cprior_mslope_10000_0_16.profM'
#basename = dir + '20130221092951_popboth_cprior_mslope_5000_5000_12/'
#basename = dir + '20130221111051_popboth_cprior_mslope_5000_5000_12/'

# ca = 0
dir = '/home/ast/read/dark/dwarf_data/data_walker/c1_100_050_100_100_core_c2_010_050_100_100_core_003_6d/'
# ca = 1
dir = '/home/ast/read/dark/dwarf_data/data_walker/c1_100_050_050_100_cusp_c2_100_050_100_100_cusp_003_6d/'

nampart = '20130327162354_cprior_mslope'
nampart = '20130327202714_cprior_mslope'
nampart = '20130403103916_cprior_mslope'

nampart = '20130408092450_cprior_mslope'
basename = dir + nampart + '/' + nampart

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

rsc = 1.#0.5
Msc = 1.#10.
sel = (radii<1500.)
radsc = radii[sel]*rsc
plot_data(radsc,M95hi[sel]*Msc,'g',lw=2)
plot_data(radsc,M68hi[sel]*Msc,'b',lw=2)
plot_data(radsc,Mmedi[sel]*Msc,'r',lw=2)
plot_data(radsc,M68lo[sel]*Msc,'b',lw=2)
plot_data(radsc,M95lo[sel]*Msc,'g',lw=2)

# plot_data(radsc,Mmax,'r',lw=2)
# plot_data(radsc,Mmin,'r',lw=2)

def readcol(filena):
    a,b,c = np.loadtxt(filena,skiprows=1,unpack=True)
    return a,b,c

#datMr,datMdat,datMerr = readcol('/home/ast/read/dark/dwarf_data/hernquist_justin/enclosedmass/unit_hern_1_enclosedmass_1000_0.txt')

#sel = (datMr<max(radii))
# plot_data(rsc*datMr[sel],Msc*datMdat[sel],'black',lw=3)
# print 'data ',ipollog(datMr,datMdat,radii)
# plot_data(rsc*radii,Msc*ipollog(datMr,datMdat,radii), lw=3)

#plot_data(rsc*radii,Msc*M_anf(radii),'black',lw=4)
#plot_data(rsc*radii,0.*ones(len(radii)),'black',lw=4)
plot_data(rsc*radii[sel],Msc*Mwalkertot(radii)[sel],'black',lw=4)
# setlims(ax1,[0.,3.],[-1.,1.])

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

