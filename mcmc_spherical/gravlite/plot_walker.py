#!/usr/bin/python
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
from pylab import *
import pdb
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline
from gl_analytic import Mwalkertot, rhowalkertot_3D
from matplotlib.backends.backend_pdf import PdfPages

# Walker data sets
import select_run
basename, prof = select_run.run()

def M_anf(r):
    return r**2/(r+1.)**2

def ipol(xin,yin,xout,smooth=1.e-9):
    rbf = Rbf(xin, yin, smooth=smooth)
    return rbf(xout)

def ipollog(xin,yin,xout,smooth=1.e-9):
    rbf = Rbf(xin, np.log10(yin), smooth=smooth)
    return 10**(rbf(xout))

def setlabels(ax,xtext,ytext):
    ax.set_xlabel(r'$'+xtext+'$')
    ax.set_ylabel(r'$'+ytext+'$')

def setlims(ax,xlim,ylim):
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

def save_plot(nam):
    plt.savefig(nam)
    return

def show_plots():
    show()
    return



print 'input'
print basename
M = np.loadtxt(basename+'prof'+prof,skiprows=0,unpack=False)

radii = M[0]
radii = radii
print 'len (M) = ',len(M)
profs = M[1:]                           # all saved profiles
#profs = M[-10:]
#profs = M[-10000::10] # only the last 1e5 profiles, thinned by 10

Mprofbins = np.transpose(profs)
radii = radii[:-1]
Mprofbins = Mprofbins[:-1]

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
Msc = 1.
sel = (radii<15000.)             # TODO: selection right?
radsc = radii[sel]*rsc


def plotGraph():
    fig = plt.figure(figsize=(3, 3), dpi=300)
    ### Plotting arrangements ###
    xlabel(r'$r\quad[\mathrm{pc}]$')
    if prof == 'M':
        ylabel(r'$M\quad[\mathrm{M}_{\odot}]$') #[10^5 M_{\odot}]')
    elif prof == 'dens':
        ylabel(r'$\rho\quad[\mathrm{M}_{\odot}/\mathrm{pc}^3]$') #[10^5 M_{\odot}]')
    elif prof=='delta1' or prof == 'delta2':
        ylabel(r'$\delta$') #[10^5 M_{\odot}]')
    elif prof == 'nu1':
        ylabel(r'$\nu_1$')
    elif prof == 'sig1':
        ylabel(r'$\sigma_1$')
        
    fill_between(radsc, M95lo[sel]*Msc, M95hi[sel]*Msc, color='black',alpha=0.2,lw=1)
    fill_between(radsc, M68lo[sel]*Msc, M68hi[sel]*Msc, color='black',alpha=0.4,lw=1)
    plot(radsc,Mmedi[sel]*Msc,'r',lw=1)
    # theoretical model
    if prof == 'M':
        plot(rsc*radii[sel],Msc*Mwalkertot(radii)[sel],'--',color='black',lw=1)
    elif prof == 'dens':
        plot(rsc*radii[sel],Msc*rhowalkertot_3D(radii)[sel],'--',color='black',lw=1)
    elif prof == 'sig1':
        # TODO: plot data as background
        rad, sig1, sigerr1 = gh.readcol(gp.files.sigfiles[1])
        fill_between(rad, sig1-sigerr1, sig1+sigerr1, color='blue', alpha=0.8,lw=1)
        plot(rad, sig1, color='blue', lw=1)
    if prof != 'delta1' and prof != 'delta2' and prof != 'sig1':
        xscale('log')
        yscale('log')
    xlim([min(radsc),max(radsc)])
    ylim([min(M95lo[sel]*Msc),max(M95hi[sel]*Msc)])
    return fig
    

def readcol(filena):
    a,b,c = np.loadtxt(filena,skiprows=1,unpack=True)
    return a,b,c


ion()
plot1 = plotGraph()
pp = PdfPages(basename + 'prof'+prof+'.pdf')
pp.savefig(plot1)

# We can also set the file's metadata via the PdfPages object:
d = pp.infodict()
d['Title'] = 'Multipage PDF'
d['Author'] = u'Pascal Steger'
d['Subject'] = 'dwarf spheroidal dark matter density profile'
d['Keywords'] = 'PdfPages multipage keywords author title subject'
d['CreationDate'] = datetime.datetime(2013,07,19)
d['ModDate'] = datetime.datetime.today()
pp.close()
ioff()

fout = open(basename+'prof'+prof+'.conf','w')
print >> fout,M95lo
print >> fout,M68lo
print >> fout,Mmedi
print >> fout,M68hi
print >> fout,M95hi
fout.close()

analyt = M_anf(radii)

import gl_helper as gh
print '# radii  lower 95%    lower 68%   median      upper 68%   upper 95%   analytic'
for i in range(len(radii)):
    print gh.pretty(radii[i]),\
          gh.pretty(M95lo[i]),\
          gh.pretty(M68lo[i]),\
          gh.pretty(Mmedi[i]),\
          gh.pretty(M68hi[i]),\
          gh.pretty(M95hi[i]),\
          gh.pretty(analyt[i])


show_plots()

