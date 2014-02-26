#!/usr/bin/env python3

##
# @file
# plotting Hernquist profiles of all accepted models

# (c) 2013 ETHZ, Pascal Steger, psteger@phys.ethz.ch

import numpy as np
from pylab import *
import pdb
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline
from gl_analytic import Mwalkertot
from matplotlib.backends.backend_pdf import PdfPages
from plots_common import *
import gl_params
gp = gl_params.Params()

dir = gp.files.machine+'/backup/'

nampart = 'working_10000_16_nobetaprior/'
nampart = 'working_2x5000_12_nobetaprior/'
nampart = 'working_1000_12_nobetaprior/'

basename = dir + nampart

print('input: ', basename)
M = np.loadtxt(basename+'profM',skiprows=1,unpack=False)
#M = np.loadtxt(basename+'.profbeta1',skiprows=0,unpack=False)

radii = M[0]
radii = radii
profs = M[1:]

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
sel = (radii<15000.)             # TODO: selection right?
radsc = radii[sel]*rsc


def plotGraph():
    fig = plt.figure()
    ### Plotting arrangements ###
    xlabel(r'$r\quad[\mathrm{pc}]$')
    ylabel(r'$M\quad[\mathrm{M}_{\odot}]$') #[10^5 M_{\odot}]')
    fill_between(radsc, M95lo[sel]*Msc, M95hi[sel]*Msc, color='black',alpha=0.2,lw=1)
    fill_between(radsc, M68lo[sel]*Msc, M68hi[sel]*Msc, color='black',alpha=0.4,lw=1)
    plot(radsc,Mmedi[sel]*Msc,'r',lw=2)
    # theoretical model
    datMr,datMdat,datMerr = readcol(dir+nampart+'unit_hern_1_enclosedmass.txt')

    plot(datMr,datMdat,'--',color='black',lw=2)
    xscale('log');    yscale('log')
    xlim([min(radsc),max(radsc)])
    ylim([min(M95lo[sel]*Msc),max(M95hi[sel]*Msc)])
    return fig
    
def readcol(filena):
    a,b,c = np.loadtxt(filena,skiprows=1,unpack=True)
    return a,b,c
ion()
plot1 = plotGraph()
pp = PdfPages(basename + '.profM.pdf')
pp.savefig(plot1)

# We can also set the file's metadata via the PdfPages object:
d = pp.infodict()
d['Title'] = 'Multipage PDF'
d['Author'] = u'Pascal Steger'
d['Subject'] = 'dwarf spheroidal dark matter density profile'
d['Keywords'] = 'PdfPages multipage keywords author title subject'
d['CreationDate'] = datetime.datetime(2013,05,06)
d['ModDate'] = datetime.datetime.today()
pp.close()
ioff()

save_profile(basename, 'M', M95lo, M68lo, Mmedi, M68hi, M95hi)

analyt = M_anf(radii)

print('# radii  lower 95%    lower 68%   median      upper 68%   upper 95%   analytic')
for i in range(len(radii)):
    print(radii[i],M95lo[i],M68lo[i],Mmedi[i],M68hi[i],M95hi[i],analyt[i])

show_plots()

