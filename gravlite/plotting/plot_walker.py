#!/usr/bin/env python3

##
# @file
# plot output profiles for Walker datasets

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
from pylab import *
import pdb
from scipy.interpolate import splrep, splev

import gl_params as gp
from gl_analytic import Mwalkertot, rhowalktot_3D, betawalker
import gl_helper as gh
from gl_project import rho_INT_Rho
from matplotlib.backends.backend_pdf import PdfPages
from plots_common import *

# Walker data sets
import select_run
basename, prof, pop = select_run.run()
# TODO: adapt prof = 'nu' and pop = 1


## get scale radii for each of the three defined cases
# @param case which set to work on. so far, three are defined
def halflightradii(case):
    if gp.case == 0:
        gamma_star1 =   0.1;    gamma_star2 =   1.0 # 1. or 0.1
        beta_star1  =   5.0;    beta_star2  =   5.0 # fixed to 5
        r_star1     = 1000.;    r_star2     = 1000. # 500 or 1000
        r_a1        =   1.0;    r_a2        =   1.0
        gamma_DM    = 0 # 0 or 1
        
    elif gp.case == 1:
        gamma_star1 =   1.0;    gamma_star2 =   1.0 # 1. or 0.1
        beta_star1  =   5.0;    beta_star2  =   5.0 # fixed to 5
        r_star1     =  500.;    r_star2     = 1000. # 500 or 1000
        r_a1        =   1.0;    r_a2        =   1.0
        gamma_DM    = 0 # core
        
    elif gp.case == 2:
        gamma_star1 =   1.0;    gamma_star2 =   1.0 # 1. or 0.1
        beta_star1  =   5.0;    beta_star2  =   5.0 # fixed to 5
        r_star1     =  500.;    r_star2     = 1000. # 500 or 1000
        r_a1        =   1.0;    r_a2        =   1.0
        gamma_DM    = 1 # cusp

    return r_star1, r_star2


print('input: ', basename)
M = np.loadtxt(basename+'prof'+prof+str(pop),skiprows=0,unpack=False)
print('len (M) = ',len(M))

radii = M[0]
profs = M[1:]                           # all saved profiles
#profs = M[-10:]
#profs = M[-10000::10] # only the last 1e5 profiles, thinned by 10

Mprofbins = np.transpose(profs)
newMprofbins = gh.sort_profiles_binwise(Mprofbins)

Mmin, M95lo, M68lo, Mmedi, M68hi, M95hi, Mmax = gh.get_median_1_2_sig(newMprofbins)

if prof == 'nu':
    Mmax=rho_INT_Rho(radii, Mmax)
    M95hi=rho_INT_Rho(radii, M95hi)
    M68hi=rho_INT_Rho(radii, M68hi)
    Mmedi=rho_INT_Rho(radii, Mmedi)
    M68lo=rho_INT_Rho(radii, M68lo)
    M95lo=rho_INT_Rho(radii, M95lo)
    Mmin=rho_INT_Rho(radii, Mmin)

# def extralast(rad, prof):
#     tck = splrep(rad[:-1],np.log(prof[:-1]),k=2,s=0.1)
#     prof = np.hstack([prof[:-1],np.exp(splev(rad[-1],tck))])
#     return prof

# M95hi = extralast(radii,M95hi)
# M68hi = extralast(radii,M68hi)
# Mmedi = extralast(radii,Mmedi)
# M68lo = extralast(radii,M68lo)
# M95lo = extralast(radii,M95lo)

    
rsc = 1.#0.5
Msc = 1.
sel = (radii<15000.)             # TODO: selection right?
radsc = radii[sel]*rsc


def read_scale():
    for i in range(gp.pops+1):
        A = np.loadtxt(gp.files.get_scale_file(i), unpack=False, skiprows=1)
        gp.Rscale.append(A[0])
        gp.Nu0rscale.append(A[1])
        gp.Nu0pc.append(A[2])
        gp.totmass.append(A[3])
        gp.maxvlos.append(A[4])


def plotGraph():
    fig = plt.figure()
    xlabel(r'$r\quad[\mathrm{pc}]$')
    if prof == 'M':
        ylabel(r'$M\quad[\mathrm{M}_{\odot}]$') #[10^5 M_{\odot}]')
    elif prof == 'rho':
        ylabel(r'$\rho\quad[\mathrm{M}_{\odot}/\mathrm{pc}^3]$') #[10^5 M_{\odot}]')
    elif prof=='beta':
        ylabel(r'$\beta_'+str(pop)+'$') #[10^5 M_{\odot}]')
    elif prof == 'nu':
        ylabel(r'$\nu_'+str(pop)+'$')
    elif prof == 'sig':
        ylabel(r'$\sigma_'+str(pop)+'$')
        
    fill_between(radsc, M95lo*Msc, M95hi*Msc, color='black',alpha=0.2,lw=1)
    fill_between(radsc, M68lo*Msc, M68hi*Msc, color='black',alpha=0.4,lw=1)
    plot(radsc,Mmedi*Msc,'r',lw=1)
    
    r1, r2 = halflightradii(gp.case)
    if prof == 'dens' or prof == 'M' or prof=='delta1' or prof=='nu1' or prof=='sig1':
        axvline(x=r1, visible=True)
    if prof == 'dens' or prof == 'M' or prof=='delta2':
        axvline(x=r2, visible=True)
    
    # theoretical model
    if prof == 'M':
        plot(radsc,Msc*Mwalkertot(radsc),'--',color='black',lw=1)
    elif prof == 'rho':
        plot(radsc,Msc*rhowalktot_3D(radsc),'--',color='black',lw=1) # 7 co
        # plot lowest profile as well
        # plot(radsc,Msc*profs[0],'.',color='orange',lw=1)
    elif prof == 'beta':
        plot(radsc,betawalker(radsc)[pop],color='black')
    elif prof == 'sig':
        rad, dummy, dummy, sig1, sigerr1 = gh.readcol5(gp.files.sigfiles[pop])
        rad *= gp.Rscale[pop]
        sigdat *= gp.maxvlos[pop]
        sigerr *= gp.maxvlos[pop]
        # plot(rad, sigdat, '--', color='b', lw=2)
        # plot(rad, sigdat-sigerr, '--', color='b', lw=1)
        # plot(rad, sigdat+sigerr, '--', color='b', lw=1)
        fill_between(rad, sigdat-sigerr, sigdat+sigerr, color='blue', alpha=0.3,lw=0)
    elif prof == 'nu':
        rad, dummy, dummy, nudat, nuerr = gh.readcol5(gp.files.nufiles[pop])
        rad *= gp.Rscale[pop]
        nudat *=gp.Nu0pc[pop]
        nuerr *=gp.Nu0pc[pop]
        fill_between(rad, nudat-nuerr, nudat+nuerr, color='blue', alpha=0.3,lw=0)


    if prof!='beta' and prof != 'sig':
        yscale('log')
    if prof == 'nu':
        yscale('log')
    if prof == "M" or prof == 'rho' or prof == 'nr':
        xscale('log')
        xlim([100.,1200.]);  ylim([0.005,1.5])
    if prof == 'beta':
        xlim([100.,1200.]);  ylim([-0.15,0.9])
    return fig



ion()
read_scale()
plot1 = plotGraph()
pp = PdfPages(basename + 'prof'+prof+str(pop)+'.pdf')
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

save_profile(basename, prof, M95lo, M68lo, Mmedi, M68hi, M95hi)

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

