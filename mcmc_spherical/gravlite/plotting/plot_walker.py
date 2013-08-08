#!/usr/bin/python
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
from pylab import *
import pdb

from scipy.interpolate import Rbf, InterpolatedUnivariateSpline

import gl_params as gp
from gl_analytic import Mwalkertot, rhowalkertot_3D, betawalker
import gl_helper as gh
from gl_int import int_surfden
from matplotlib.backends.backend_pdf import PdfPages

# Walker data sets
import select_run
basename, prof = select_run.run()

def M_anf(r):
    return r**2/(r+1.)**2

def save_plot(nam):
    plt.savefig(nam)
    return

def show_plots():
    show()
    return


def halflightradii(case):
    if gp.walkercase == 0:
        gamma_star1 =   0.1;    gamma_star2 =   1.0 # 1. or 0.1
        beta_star1  =   5.0;    beta_star2  =   5.0 # fixed to 5
        r_star1     = 1000.;    r_star2     = 1000. # 500 or 1000
        r_a1        =   1.0;    r_a2        =   1.0
        gamma_DM    = 0 # 0 or 1
        
    elif gp.walkercase == 1:
        gamma_star1 =   1.0;    gamma_star2 =   1.0 # 1. or 0.1
        beta_star1  =   5.0;    beta_star2  =   5.0 # fixed to 5
        r_star1     =  500.;    r_star2     = 1000. # 500 or 1000
        r_a1        =   1.0;    r_a2        =   1.0
        gamma_DM    = 0 # core
        
    elif gp.walkercase == 2:
        gamma_star1 =   1.0;    gamma_star2 =   1.0 # 1. or 0.1
        beta_star1  =   5.0;    beta_star2  =   5.0 # fixed to 5
        r_star1     =  500.;    r_star2     = 1000. # 500 or 1000
        r_a1        =   1.0;    r_a2        =   1.0
        gamma_DM    = 1 # cusp

    return r_star1, r_star2



print 'input'
print basename
M = np.loadtxt(basename+'prof'+prof,skiprows=0,unpack=False)
print 'len (M) = ',len(M)

radii = M[0]
profs = M[1:]                           # all saved profiles
#profs = M[-10:]
#profs = M[-10000::10] # only the last 1e5 profiles, thinned by 10


Mprofbins = np.transpose(profs)
# radii = radii[:-1]
# Mprofbins = Mprofbins[:-1]


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

if prof == 'nu1':
    Mmax=int_surfden(radii, Mmax)
    M95hi=int_surfden(radii, M95hi)
    M68hi=int_surfden(radii, M68hi)
    Mmedi=int_surfden(radii, Mmedi)
    M68lo=int_surfden(radii, M68lo)
    M95lo=int_surfden(radii, M95lo)
    Mmin=int_surfden(radii, Mmin)

rsc = 1.#0.5
Msc = 1.
sel = (radii<15000.)             # TODO: selection right?
radsc = radii[sel]*rsc


def read_scale():
    for i in range(gp.pops+1):
        A = np.loadtxt(gp.files.get_scale_file(i), unpack=False, skiprows=1)
        gp.rcore_2D.append(A[0])
        gp.dens0rcore_2D.append(A[1])
        gp.dens0pc_2D.append(A[2])
        gp.totmass.append(A[3])
        gp.maxvlos.append(A[4])


def plotGraph():
    fig = plt.figure(figsize=(5, 3), dpi=300)
    ### Plotting arrangements ###
    xlabel(r'$r\quad[\mathrm{pc}]$')
    if prof == 'M':
        ylabel(r'$M\quad[\mathrm{M}_{\odot}]$') #[10^5 M_{\odot}]')
    elif prof == 'dens':
        ylabel(r'$\rho\quad[\mathrm{M}_{\odot}/\mathrm{pc}^3]$') #[10^5 M_{\odot}]')
    elif prof=='delta1' or prof == 'delta2':
        ylabel(r'$\beta$') #[10^5 M_{\odot}]')
    elif prof == 'nu1':
        ylabel(r'$\nu_1$')
    elif prof == 'sig1':
        ylabel(r'$\sigma_1$')
        
    fill_between(radsc, M95lo*Msc, M95hi*Msc, color='black',alpha=0.2,lw=1)
    fill_between(radsc, M68lo*Msc, M68hi*Msc, color='black',alpha=0.4,lw=1)
    plot(radsc,Mmedi*Msc,'r',lw=1)
    
    r1, r2 = halflightradii(gp.walkercase)
    if prof == 'dens' or prof == 'M' or prof=='delta1' or prof=='nu1' or prof=='sig1':
        axvline(x=r1, visible=True)
    if prof == 'dens' or prof == 'M' or prof=='delta2':
        axvline(x=r2, visible=True)
    
    # theoretical model
    if prof == 'M':
        plot(radsc,Msc*Mwalkertot(radsc),'--',color='black',lw=1)
    elif prof == 'dens':
        plot(radsc,Msc*rhowalkertot_3D(radsc),'--',color='black',lw=1)
        plot(radsc,Msc*profs[0],'.',color='orange',lw=1)
    elif prof == 'delta1':
        plot(radsc,betawalker(radsc)[0],color='black')
    elif prof == 'delta2':
        plot(radsc,betawalker(radsc)[1],color='black')

    elif prof == 'sig1':
        rad, sig1, sigerr1 = gh.readcol(gp.files.sigfiles[1])
        rad *= gp.rcore_2D[1]
        sig1*=gp.maxvlos[1]
        sigerr1*=gp.maxvlos[1]
        # plot(rad, sig1, '--', color='b', lw=2)
        # plot(rad, sig1-sigerr1, '--', color='b', lw=1)
        # plot(rad, sig1+sigerr1, '--', color='b', lw=1)
        fill_between(rad, sig1-sigerr1, sig1+sigerr1, color='blue', alpha=0.3,lw=0)
    elif prof == 'nu1':
        rad, nu1, nuerr1 = gh.readcol(gp.files.nufiles[1])
        rad *= gp.rcore_2D[1]
        nu1*=gp.dens0pc_2D[1]
        nuerr1*=gp.dens0pc_2D[1]
        fill_between(rad, nu1-nuerr1, nu1+nuerr1, color='blue', alpha=0.3,lw=0)


    if prof != 'delta1' and prof != 'delta2' and prof != 'sig1' and prof != 'nu1':
        # xscale('log')
        yscale('log')
    if prof == 'nu1':
        yscale('log')
    if prof == "M":
        xscale('log')
    xlim([min(radsc),max(radsc)])
    ylim([min(M95lo[sel]*Msc),max(M95hi[sel]*Msc)])
    return fig
    

def readcol(filena):
    a,b,c = np.loadtxt(filena,skiprows=1,unpack=True)
    return a,b,c


ion()
read_scale()
plot1 = plotGraph()
pp = PdfPages(basename + 'prof'+prof+'.pdf')
pp.savefig(plot1)


# We can also set the file's metadata via the PdfPages object:
d = pp.infodict()
d['Title'] = 'Profile for Dwarf Galaxy Mock Datasets'
d['Author'] = u'Pascal Steger'
d['Subject'] = 'dwarf spheroidal dark matter density and other profiles'
d['Keywords'] = 'PdfPages multipage keywords Pascal Steger 2013'
d['CreationDate'] = datetime.datetime.today()
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

