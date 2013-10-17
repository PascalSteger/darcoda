#!/usr/bin/env python3
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
from pylab import *
import pdb
from scipy.interpolate import splrep, splev

import gl_params as gp
from gl_analytic import Mwalkertot, rhogaiatot_3D, betagaia
import gl_helper as gh
from gl_project import rho_INT_Rho
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

def readcol(filena):
    a,b,c = np.loadtxt(filena,skiprows=1,unpack=True)
    return a,b,c

def halflightradius(case):
    beta_star1 = 5; r_DM = 1000
    if gp.case == 1:
        gamma_star1=0.1;r_star1=100;r_a1=100;gamma_DM=1;rho0=0.064
    elif gp.case == 2:
        gamma_star1=0.1;r_star1=250;r_a1=250;gamma_DM=0;rho0=0.400
    elif gp.case == 3:
        gamma_star1=0.1;r_star1=250;r_a1=np.inf;gamma_DM=1;rho0=0.064
    elif gp.case == 4:
        gamma_star1=0.1;r_star1=1000;r_a1=np.inf;gamma_DM=0;rho0=0.400
    elif gp.case == 5:
        gamma_star1=1.0;r_star1=100;r_a1=100;gamma_DM=1;rho0=0.064
    elif gp.case == 6:
        gamma_star1=1.0;r_star1=250;r_a1=250;gamma_DM=0;rho0=0.400
    elif gp.case == 7:
        gamma_star1=1.0;r_star1=250;r_a1=np.inf;gamma_DM=1;rho0=0.064
    elif gp.case == 8:
        gamma_star1=1.0;r_star1=1000;r_a1=np.inf;gamma_DM=0;rho0=0.400
            
    return r_star1

print('input: ',basename)
M = np.loadtxt(basename+'prof'+prof,skiprows=0,unpack=False)
print('len (M) = ',len(M))

radii = M[0]
profs = M[1:]                           # all saved profiles
#profs = M[-10:]
#profs = M[-10000::10] # only the last 1e5 profiles, thinned by 10


Mprofbins = np.transpose(profs)
# radii = radii[1:] # TODO: debug: allowin for first bin suddenly makes
# overall resemblance change!
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
        gp.rcore_2D.append(A[0])
        gp.dens0rcore_2D.append(A[1])
        gp.dens0pc_2D.append(A[2])
        gp.totmass.append(A[3])
        gp.maxvlos.append(A[4])


def plotGraph():
    fig = plt.figure()
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
    
    r1 = halflightradius(gp.case)
    if prof == 'dens' or prof == 'M' or prof=='delta1' or prof=='nu1' or prof=='sig1':
        axvline(x=r1, visible=True)
    
    # theoretical model
    if prof == 'M':
        plot(radsc,Msc*Mwalkertot(radsc),'--',color='black',lw=1)
    elif prof == 'dens':
        plot(radsc,Msc*rhogaiatot_3D(radsc),'--',color='black',lw=1) # 7 co
        # plot lowest profile as well
        # plot(radsc,Msc*profs[0],'.',color='orange',lw=1)
    elif prof == 'delta1':
        plot(radsc, betagaia(radsc), color='black')
    elif prof == 'sig1':
        rad, dummy, dummy, sig1, sigerr1 = gh.readcol5(gp.files.sigfiles[1])
        rad *= gp.rcore_2D[1]
        sig1*=gp.maxvlos[1]
        sigerr1*=gp.maxvlos[1]
        # plot(rad, sig1, '--', color='b', lw=2)
        # plot(rad, sig1-sigerr1, '--', color='b', lw=1)
        # plot(rad, sig1+sigerr1, '--', color='b', lw=1)
        fill_between(rad, sig1-sigerr1, sig1+sigerr1, color='blue', alpha=0.3,lw=0)
    elif prof == 'nu1':
        rad, dummy, dummy, nu1, nuerr1 = gh.readcol5(gp.files.nufiles[1])
        rad *= gp.rcore_2D[1]
        # TODO: missing factor 10 somewhere..
        pdb.set_trace()
        nu1 *= gp.dens0pc_2D[1]
        nuerr1 *= gp.dens0pc_2D[1]
        fill_between(rad, nu1-nuerr1, nu1+nuerr1, color='blue', alpha=0.3,lw=0)

    if prof!='delta1' and prof!='delta2' and prof != 'sig1' and prof != 'nu1':
        # xscale('log')
        yscale('log')
    if prof == 'nu1':
        yscale('log')
    if prof == "M" or prof == 'dens':
        xscale('log')
        #xlim([100.,1200.]);  ylim([0.005,1.5])
    # if prof == 'delta1':
    #     xlim([100.,1200.]);  ylim([-0.15,0.9])

    return fig
    



ion()
read_scale()
plot1 = plotGraph()
pp = PdfPages(basename + 'prof'+prof+'.pdf')
pp.savefig(plot1)


# We can also set the file's metadata via the PdfPages object:
d = pp.infodict()
d['Title'] = 'Profile for Dwarf Galaxy Mock Datasets'
d['Author'] = 'Pascal Steger'
d['Subject'] = 'dwarf spheroidal dark matter density and other profiles'
d['Keywords'] = 'PdfPages multipage keywords Pascal Steger 2013'
d['CreationDate'] = datetime.datetime.today()
d['ModDate'] = datetime.datetime.today()
pp.close()
ioff()

fout = open(basename+'prof'+prof+'.conf','w')
print(M95lo, file=fout)
print(M68lo, file=fout)
print(Mmedi, file=fout)
print(M68hi, file=fout)
print(M95hi, file=fout)
fout.close()

analyt = M_anf(radii)

import gl_helper as gh
print('# radii  lower 95%    lower 68%   median      upper 68%   upper 95%   analytic')
for i in range(len(radii)):
    print(gh.pretty(radii[i]),\
          gh.pretty(M95lo[i]),\
          gh.pretty(M68lo[i]),\
          gh.pretty(Mmedi[i]),\
          gh.pretty(M68hi[i]),\
          gh.pretty(M95hi[i]),\
          gh.pretty(analyt[i]))


show_plots()

