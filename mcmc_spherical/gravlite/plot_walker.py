#!/usr/bin/python
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
from pylab import *
import pdb
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline
from gl_analytic import Mwalkertot, rhowalkertot_3D
from matplotlib.backends.backend_pdf import PdfPages

# Walker data sets
base = '/home/psteger/sci/dwarf_data/'
base = '/home/ast/read/dark/dwarf_data/'
# ca = 0:
dir = base + 'data_walker/c1_100_050_100_100_core_c2_010_050_100_100_core_003_6d/'
dir = base + 'data_walker/c1_010_050_100_100_core_c2_100_050_100_100_core_003_6d/'
# ca = 1:
dir = base + 'data_walker/c1_100_050_050_100_core_c2_100_050_100_100_core_003_6d/'
# ca = 2:
dir = base + 'data_walker/c1_100_050_050_100_cusp_c2_100_050_100_100_cusp_003_6d/'


# nampart = '20130425120348_cprior_mslope_rprior' # ca2 working fine for 1000 iterations
# nampart = '20130426090433_cprior_mslope_rprior' # ca2 too high mass at high radii, 50k
# nampart = '20130426120258_cprior_mslope_rprior' # ca2 better mass? no, too low overall
# nampart = '20130426161637_cprior_nulog_denslog_mslope_rprior' # ca2 new denslog:  too high mass, especially around middle radii
# nampart = '20130426165536_cprior_nulog_denslog_mslope_rprior' # ca2 and up to 100000 its: works
# nampart = '20130429110855_cprior_nulog_denslog_mslope_rprior' # 50k steps: works
# nampart = '20130502080536_cprior_nulog_denslog_mslope_rprior' # ca0 10: London failed 2.5k, too high mass
#nampart = '20130510090417_case_1_0_0_cprior_nulog_denslog_mslope_rprior' # ca2, core, running
# nampart = '20130621123935_case_2_10000_0_cprior_nulog_denslog_mslope_rprior' #ca2, core
# nampart = '20130426165536_cprior_nulog_denslog_mslope_rprior'
# nampart = '20130426133539_cprior_mslope_rprior'
# nampart = '20130702113536_case_2_10000_0_cprior_nulog_denslog_mslope_rprior' # new cusp
#nampart = '20130621123935_case_2_10000_0_cprior_nulog_denslog_mslope_rprior' # Surrey largest
#nampart = '20130621121716_case_2_10000_0_cprior_nulog_denslog_mslope_rprior' # Surrey 2nd
#nampart = '20130717081725_case_2_10000_0_cprior_nulog_denslog_mslope_rprior' # first mio cusp
#nampart = '20130717085139_case_1_10000_0_cprior_nulog_denslog_mslope_rprior'

#nampart = '20130717164820_case_2_10000_0_cprior_nulog_denslog_mslope_rprior' # cusp 2nd case

nampart = '20130718090029_case_2_10000_0_cprior_nulog_denslog_mslope_rprior' # new stepsize, cusp
nampart = '20130718103257_case_2_10000_0_cprior_nulog_denslog_mslope_rprior' # twiddle 1
nampart = '20130718120843_case_2_10000_0_cprior_nulog_denslog_delta_mslope_rprior' # fast conv
nampart = '20130718132442_case_2_10000_0_cprior_nulog_denslog_mslope_rprior' # finally fine

# nampart = '20130717164748_case_1_10000_0_cprior_nulog_denslog_mslope_rprior' # core 2nd case
# nampart = '20130718133920_case_1_10000_0_cprior_nulog_denslog_mslope_rprior' # finally core
# nampart = '20130718161434_case_1_10000_0_cprior_nulog_denslog_mslope_rprior' # debugged core
# nampart = '20130718183721_case_1_10000_0_cprior_nulog_denslog_mslope_rprior' # 2nd try debug

basename = dir + nampart + '/' + nampart

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
prof = 'delta1' # dens, M, delta1, delta2
M = np.loadtxt(basename+'.prof'+prof,skiprows=0,unpack=False)

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
    elif prof=='delta':
        ylabel(r'$\delta$') #[10^5 M_{\odot}]')
    fill_between(radsc, M95lo[sel]*Msc, M95hi[sel]*Msc, color='black',alpha=0.2,lw=1)
    fill_between(radsc, M68lo[sel]*Msc, M68hi[sel]*Msc, color='black',alpha=0.4,lw=1)
    plot(radsc,Mmedi[sel]*Msc,'r',lw=1)
    # theoretical model
    if prof == 'M':
        plot(rsc*radii[sel],Msc*Mwalkertot(radii)[sel],'--',color='black',lw=1)
    elif prof == 'dens':
        plot(rsc*radii[sel],Msc*rhowalkertot_3D(radii)[sel],'--',color='black',lw=1)
    # xscale('log')
    if prof != 'delta1' and prof != 'delta2':
        yscale('log')
    xlim([min(radsc),max(radsc)])
    ylim([min(M95lo[sel]*Msc),max(M95hi[sel]*Msc)])
    return fig
    
def readcol(filena):
    a,b,c = np.loadtxt(filena,skiprows=1,unpack=True)
    return a,b,c
ion()
plot1 = plotGraph()
pp = PdfPages(basename + '.prof'+prof+'.pdf')
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

fout = open(basename+'.prof'+prof+'.conf','w')
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


show_plots()

