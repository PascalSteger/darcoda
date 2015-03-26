#!/usr/bin/env ipython3

## @file
# collect profiles and perform actions on them

# (c) 2015 ETHZ Pascal Steger, pascal@steger.aero

import numpy as np
import numpy.random as npr
import pdb

import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev
plt.ioff()

from gi_class_profiles import Profiles
import gi_output as go
import gi_helper as gh
import gi_analytic as ga
import gi_project as gip

USE_ALL=False

def unit(prof):
    if prof == 'rho' or prof == 'nu':
        unit = '[Msun/pc^3]'
    elif prof == 'M':
        unit = '[Msun]'
    elif prof == 'nr':
        unit = '[ln(Msun/pc^3)/ln(pc)]'
    elif prof == 'beta':
        unit = '[1]'
    elif prof == 'betastar':
        unit = '[1]'
    elif prof == 'Sig':
        unit = '[(Msun/pc^3)]'
    elif prof == 'sig':
        unit = '[km/s]'
    else:
        unit = 'a.u.'
    return unit
## \fn unit(prof)
# return string with units for profile
# @param prof string for profile

class ProfileCollection():
    def __init__(self, pops, nepol):
        self.chis = []
        self.goodchi = []  # will contain all chi^2 for models in self.subset
        self.goodprof = [] # will contain the corresponding profiles
        self.profs = []
        self.subset = [0., np.inf]
        self.x0 = np.array([])
        self.binmin = np.array([])
        self.binmax = np.array([])
        self.Mmin  = Profiles(pops, nepol)
        self.M99lo = Profiles(pops, nepol)
        self.M95lo = Profiles(pops, nepol)
        self.M68lo = Profiles(pops, nepol)
        self.Mmedi = Profiles(pops, nepol)
        self.M68hi = Profiles(pops, nepol)
        self.M95hi = Profiles(pops, nepol)
        self.M99hi = Profiles(pops, nepol)
        self.Mmax  = Profiles(pops, nepol)
        self.analytic = Profiles(pops, 100)
        initrange = [1e30, -1e30]
        self.ranges = {'rho0': initrange,
                       'nr0': initrange,
                       'M0': initrange,
                       'nu0': initrange,
                       'nu1': initrange,
                       'nu2': initrange,
                       'Sig0': initrange,
                       'Sig1': initrange,
                       'Sig2': initrange,
                       'sig1': initrange,
                       'sig2': initrange,
                       'nrnu0': initrange,
                       'nrnu1': initrange,
                       'nrnu2': initrange,
                       'beta1': initrange,
                       'beta2': initrange,
                       'betastar1': initrange,
                       'betastar2': initrange}
    ## \fn __init__(self, pops, nepol)
    # constructor
    # @param pops number of populations
    # @param nepol number of bins

    def broaden_lim(self, prof, pop, mini, maxi):
        dmin, dmax = self.ranges[prof+str(pop)]
        dmin = min(dmin, mini)
        dmax = max(dmax, maxi)
        self.ranges[prof+str(pop)] = [dmin, dmax]
        # overrule settings for n(r) plot
        if prof == 'nr':
            self.ranges[prof+str(pop)] = [-0.5, 4.0]
        return
    ## \fn broaden_lim(self, prof, pop, mini, maxi)
    # broaden the plotting range
    # @param prof string
    # @param pop int
    # @param mini new minimum
    # @param maxi new maximum

    def add(self, prof):
        if np.isnan(prof.chi2):
            return
        self.chis.append(prof.chi2)
        self.profs.append(prof)
        return
    ## \fn add(self, prof)
    # add a profile
    # @param prof Profile to add

    def merge(self, pc):
        self.chis = np.hstack([self.chis, pc.chis])
        self.profs = np.hstack([self.profs, pc.profs])
    ## \fn merge(self, pc)
    # merge another profile collection into this one

    def cut_subset(self):
        self.chis = np.array(self.chis) # better to work with numpy arrays for subsections
        minchi = min(self.chis)
        if USE_ALL:
            maxchi = max(self.chis)
        else:
            counts, bins = np.histogram(np.log10(self.chis), bins=np.sqrt(len(self.chis)))
            posmax = np.argmin(np.abs(counts-max(counts)))
            for k in range(posmax, len(counts)-1):
                if counts[k+1] > counts[k] or counts[k+1]==0 or counts[k+1]<=max(counts)/10:
                    break
            maxchi = min(30*min(self.chis),10**(bins[k+1]))
        self.subset = [minchi, maxchi]
        print('chi^2 subset: ',minchi, maxchi)
    ## \fn cut_subset(self)
    # set subset to [0, 10*min(chi)] (or 30* minchi, or any value wished)

    def set_x0(self, x0, Binmin, Binmax):
        self.binmin = Binmin
        self.binmax = Binmax
        self.x0 = x0
        self.Mmin.x0  = x0
        self.M99lo.x0 = x0
        self.M95lo.x0 = x0
        self.M68lo.x0 = x0
        self.Mmedi.x0 = x0
        self.M68hi.x0 = x0
        self.M95hi.x0 = x0
        self.M99hi.x0 = x0
        self.Mmax.x0  = x0
    ## \fn set_x0(self, x0)
    # set radii for 1sigma, 2sigma profiles
    # @param x0 radii in pc
    # @param Binmin minimal radii of bins in [pc]
    # @param Binmax maximal radii of bins in [pc]

    def sort_prof(self, prof, pop, gp):
        self.goodprof = []
        self.goodchi = []
        for k in range(len(self.profs)):
            if self.subset[0] <= self.chis[k] <= self.subset[1]:
                self.goodprof.append(self.profs[k].get_prof(prof, pop))
                self.goodchi.append(self.chis[k])
        tmp = gh.sort_profiles_binwise(np.array(self.goodprof).transpose()).transpose()
        ll = len(tmp)
        #norm = 1
        #if prof == 'Sig':
        #    norm = gh.ipol_rhalf_log(gp.xepol, tmp[ll/2], gp.Xscale[0])
        self.Mmin.set_prof(prof,  tmp[0],       pop, gp)
        self.M99lo.set_prof(prof, tmp[ll*0.01], pop, gp)
        self.M95lo.set_prof(prof, tmp[ll*0.05], pop, gp)
        self.M68lo.set_prof(prof, tmp[ll*0.32], pop, gp)
        self.Mmedi.set_prof(prof, tmp[ll/2],    pop, gp)
        self.M68hi.set_prof(prof, tmp[ll*0.68], pop, gp)
        self.M95hi.set_prof(prof, tmp[ll*0.95], pop, gp)
        self.M99hi.set_prof(prof, tmp[ll*0.99], pop, gp)
        self.Mmax.set_prof(prof,  tmp[-1],      pop, gp)
        return tmp
    ## \fn sort_prof(self, prof, pop, gp)
    # sort the list of prof-profiles, and store the {1,2}sigma, min, medi, max in the appropriate place
    # @param prof profile identifier, 'rho', 'beta', ...
    # @param pop population identifier
    # @param gp global parameters
    # @return tmp profiles sorted binwise

    def calculate_J(self, gp):
        if len(self.profs)>0:
            for i in range(len(self.profs)):
                Sigprof = gip.rho_INT_Sig(gp.xepol, self.profs[i].get_prof('rho', 0), gp)
                Jprof = gip.Jpar(gp.xepol, Sigprof, gp)
                # add 3 extension bins
                tck = splrep(np.log(gp.xepol[:-gp.nexp]), np.log(Jprof), k=1, s=0.1)
                Jext = np.exp(splev(np.log(gp.xepol[-gp.nexp:]), tck))
                Jfull = np.hstack([Jprof, Jext])
                for k in range(gp.nepol):
                    Jfull[k] = max(0, Jfull[k])
                self.profs[i].set_prof('J', Jfull, 0, gp)
        else:
            gh.LOG(1, 'len(self.profs) == 0, did not calculate self.profs.J')
    ## \fn calculate_J(self, gp)
    # calculate J from Sig from rho
    # @param gp global parameters

    def sort_profiles(self, gp):
        self.sort_prof('rho', 0, gp)
        self.sort_prof('M', 0, gp)
        self.sort_prof('Sig', 0, gp)
        self.sort_prof('nr', 0, gp)
        self.sort_prof('nrnu', 0, gp)
        self.sort_prof('nu', 0, gp)
        for pop in np.arange(1, gp.pops+1):
            self.sort_prof('betastar', pop, gp)
            self.sort_prof('beta', pop, gp)
            self.sort_prof('nu', pop, gp)
            self.sort_prof('nrnu', pop, gp)
            self.sort_prof('Sig', pop, gp)
            self.sort_prof('sig', pop, gp)
        return
    ## \fn sort_profiles(self, gp)
    # sort all profiles, in a parallel way
    # @param gp global parameters

    def set_analytic(self, x0, gp):
        r0 = x0 # [pc], spherical case
        self.analytic.x0 = r0
        anbeta = []; annu = [];  anSig = []
        if gp.investigate == 'gaia':
            anrho = ga.rho_gaia(r0, gp)[0]
            anM = gip.rho_SUM_Mr(r0, anrho)
            annr = ga.nr3Dtot_gaia(r0, gp)
            tmp_annu = ga.rho_gaia(r0, gp)[1]
            annu.append( tmp_annu )
            anSig.append( gip.rho_INT_Sig(r0, tmp_annu, gp) )
            for pop in np.arange(1, gp.pops+1):
                beta = ga.beta_gaia(r0, gp)[pop]
                anbeta.append(beta)
                nu = ga.rho_gaia(r0,gp)[pop]
                annu.append(nu)
                anSig.append(gip.rho_INT_Sig(r0, nu, gp))
        elif gp.investigate == 'walk':
            anrho = ga.rho_walk(r0, gp)[0]
            anM = gip.rho_SUM_Mr(r0, anrho)
            annr = ga.nr3Dtot_deriv_walk(r0, gp) # TODO too high in case of core
            tmp_annu = ga.rho_walk(r0, gp)[1]
            annu.append( tmp_annu )
            anSig.append( gip.rho_INT_Sig(r0, tmp_annu, gp) )
            for pop in np.arange(1, gp.pops+1):
                beta = ga.beta_walk(r0, gp)[pop]
                anbeta.append(beta)
                nu = ga.rho_walk(r0, gp)[pop]
                dum,dum,dum,nudat,nuerr = np.transpose(np.loadtxt(gp.files.nufiles[pop], unpack=False, skiprows=1))
                locrhalf = np.argmin(abs(gp.xipol-gp.Xscale[pop]))
                nuhalf = nudat[locrhalf]*gp.nu0pc[pop]
                annuhalf = nu[np.argmin(abs(r0-locrhalf))]
                annu.append(nu*nuhalf/annuhalf)
                dum,dum,dum,Sigdat,Sigerr = np.transpose(np.loadtxt(gp.files.Sigfiles[pop], unpack=False, skiprows=1))
                locrhalf = np.argmin(abs(gp.xipol-gp.Xscale[pop]))
                Sighalf = Sigdat[locrhalf]*gp.Sig0pc[pop]
                Sig = gip.rho_INT_Sig(r0, nu, gp)
                anSighalf = Sig[np.argmin(abs(r0-locrhalf))]
                anSig.append(Sig*Sighalf/anSighalf)
        elif gp.investigate == 'triax':
            anrho = ga.rho_triax(r0, gp) # one and only
            anM = gip.rho_SUM_Mr(r0, anrho)
            annr = ga.nr3Dtot_deriv_triax(r0, gp)
            tmp_annu = ga.rho_triax(r0, gp) # TODO, M/L=1 assumed here, wrong
            annu.append(tmp_annu)
            anSig.append( gip.rho_INT_Sig(r0, tmp_annu, gp))
            for pop in np.arange(1, gp.pops+1):
                beta = ga.beta_triax(r0)
                anbeta.append(beta)
                nu = ga.rho_triax(r0, gp) # TODO, assumes M/L=1
                annu.append(nu)
                anSig.append( gip.rho_INT_Sig(r0, nu, gp))
        self.analytic.set_prof('rho', anrho, 0, gp)
        self.analytic.set_prof('M', anM, 0, gp)
        self.analytic.set_prof('nr', annr, 0, gp)
        self.analytic.set_prof('nu', annu[0], 0, gp)
        self.analytic.set_prof('nrnu', -gh.derivipol(np.log(annu[0]), np.log(r0)), 0, gp)
        self.analytic.set_prof('Sig', anSig[0], 0, gp)
        for pop in np.arange(1, gp.pops+1):
            self.analytic.set_prof('beta', anbeta[pop-1], pop, gp)
            self.analytic.set_prof('betastar', anbeta[pop-1]/(2.-anbeta[pop-1]), pop, gp)
            self.analytic.set_prof('nu', annu[pop], pop, gp)
            nrnu = -gh.derivipol(np.log(annu[pop]), np.log(r0))
            self.analytic.set_prof('nrnu', nrnu, pop, gp)
            self.analytic.set_prof('Sig', anSig[pop] , pop, gp)#/ Signorm, pop, gp)
            self.analytic.set_prof('sig', -np.ones(len(r0)), pop, gp)
        return
    ## \fn set_analytic(x0, gp)
    # set analytic curves (later shown in blue)
    # @param x0 radius in [pc]
    # @param gp global parameters

    def write_prof(self, basename, prof, pop, gp):
        output = go.Output()
        uni = unit(prof)
        output.add('radius (models) [pc]', gp.xepol)
        output.add('M 99% CL low ' + uni, self.M99lo.get_prof(prof, pop))
        output.add('M 95% CL low ' + uni, self.M95lo.get_prof(prof, pop))
        output.add('M 68% CL low ' + uni, self.M68lo.get_prof(prof, pop))
        output.add('M median '     + uni, self.Mmedi.get_prof(prof, pop))
        output.add('M 68% CL high '+ uni, self.M68hi.get_prof(prof, pop))
        output.add('M 95% CL high '+ uni, self.M95hi.get_prof(prof, pop))
        output.add('M 99% CL high '+ uni, self.M99hi.get_prof(prof, pop))
        output.write(basename+'output/ascii/prof_'+prof+'_'+str(pop)+'.ascii')
        if (gp.investigate =='walk' or gp.investigate=='gaia') \
           and (prof != 'Sig'):
            out_an = go.Output()
            out_an.add('radius [pc]', self.analytic.x0)
            out_an.add('analytic profile', self.analytic.get_prof(prof, pop))
            out_an.write(basename+'output/analytic/prof_'+prof+'_'+str(pop)+'.analytic')
    ## \fn write_prof(self, basename, prof, pop, gp)
    # write output file for a single profile
    # @param basename directory string
    # @param prof string for profile. rho, nr, betastar, ...
    # @param pop population int
    # @param gp global parameters

    def write_chi2(self, basename, edges, bins):
        output = go.Output()
        output.add('edges', edges[1:])
        output.add('bins', bins)
        output.write(basename+'output/prof_chi2_0.ascii')
        return
    ## \fn write_chi2(self, basename, edges, bins)
    # write ascii file with chi2 information
    # @param basename directory string
    # @param edges array of edges
    # @param bins y-values

    def write_all(self, basename, gp):
        self.write_prof(basename, 'rho', 0, gp)
        self.write_prof(basename, 'J', 0, gp)
        if gp.geom == 'sphere':
            self.write_prof(basename, 'M', 0, gp)
        self.write_prof(basename, 'Sig', 0, gp)
        self.write_prof(basename, 'nu', 0, gp)
        self.write_prof(basename, 'nrnu', 0, gp)
        if gp.geom=='sphere':
            self.write_prof(basename, 'nr', 0, gp)
        for pop in np.arange(1, gp.pops+1):
            self.write_prof(basename, 'betastar', pop, gp)
            self.write_prof(basename, 'beta', pop, gp)
            self.write_prof(basename, 'Sig', pop, gp)
            self.write_prof(basename, 'nu', pop, gp)
            self.write_prof(basename, 'nrnu', pop, gp)
            self.write_prof(basename, 'sig', pop, gp)
    ## \fn write_all(self, basename, gp)
    # write output files for all profiles
    # @param basename directory string
    # @param gp global parameters

    def plot_N_samples(self, ax, prof, pop):
        k=0
        while k<30:
            ind = npr.randint(len(self.profs))
            if self.subset[0] <= self.chis[ind] <= self.subset[1]:
                lp  = self.profs[ind]
                ax.plot(self.Mmedi.x0, lp.get_prof(prof, pop), color='black', alpha=0.1, lw=0.25)
                k += 1
        return
    ## \fn plot_N_samples(self, ax, prof, pop)
    # plot 30 sample profiles on top of the model {1,2} sig boundaries
    # @param ax axis object to draw upon
    # @param prof string
    # @param pop integer for which population

    def plot_data(self, ax, basename, prof, pop, gp):
        output = go.Output()
        r0 = gp.xipol # [pc]
        output.add('radius (data) [pc]', r0)
        if prof == 'Sig':
            # get 2D data here
            # Rbin [Xscale], Binmin [Xscale], Binmax [Xscale], Sig(R)/Sig(0) [1], error [1]
            dum,dum,dum,Sigdat,Sigerr = np.transpose(np.loadtxt(gp.files.Sigfiles[pop], unpack=False, skiprows=1))
            Sigdat *= gp.Sig0pc[pop] # [Msun/pc^2]
            Sigerr *= gp.Sig0pc[pop] # [Msun/pc^2]
            #Signorm = gh.ipol_rhalf_log(gp.xipol, Sigdat, gp.Xscale[pop])
            output.add('data [Msun/pc^2]', Sigdat)
            output.add('error [Msun/pc^2]', Sigerr)
            output.add('data - error [Msun/pc^2]', Sigdat-Sigerr)
            output.add('data + error [Msun/pc^2]', Sigdat+Sigerr)
            ax.fill_between(r0, Sigdat-Sigerr, Sigdat+Sigerr, color='blue', alpha=0.3, lw=1)
            self.broaden_lim('Sig', pop, min(Sigdat-Sigerr)/2, 2*max(Sigdat+Sigerr))
        elif prof == 'nu':
            # get 3D data here
            # Rbin [Xscale], Binmin [Xscale], Binmax [Xscale], nu(R)/nu(0) [1], error [1]
            dum,dum,dum,nudat,nuerr = np.transpose(np.loadtxt(gp.files.nufiles[pop], unpack=False, skiprows=1))
            nudat *= gp.nu0pc[pop] # [Msun/pc^2]
            nuerr *= gp.nu0pc[pop] # [Msun/pc^2]
            output.add('data [Msun/pc^3]', nudat)
            output.add('error [Msun/pc^3]', nuerr)
            output.add('data - error [Msun/pc^2]', nudat-nuerr)
            output.add('data + error [Msun/pc^2]', nudat+nuerr)
            ax.fill_between(r0, nudat-nuerr, nudat+nuerr, color='blue', alpha=0.3, lw=1)
            self.broaden_lim('nu', pop, min(nudat-nuerr)/2., 2*max(nudat+nuerr))
        elif prof == 'sig':
            DATA = np.transpose(np.loadtxt(gp.files.sigfiles[pop], unpack=False, skiprows=1))
            sigdat = DATA[4-1] # [maxsiglosi]
            sigerr = DATA[5-1] # [maxsiglosi]
            sigdat *= gp.maxsiglos[pop]  # [km/s]
            sigerr *= gp.maxsiglos[pop]  # [km/s]
            output.add('data [km/s]', sigdat)
            output.add('error [km/s]', sigerr)
            output.add('data - error [km/s]', sigdat-sigerr)
            output.add('data + error [km/s]', sigdat+sigerr)
            ax.fill_between(r0, sigdat-sigerr, sigdat+sigerr, color='blue', alpha=0.3, lw=1)
            self.broaden_lim('sig', pop, 0., 2*max(sigdat+sigerr))
        output.write(basename+'output/data/prof_'+prof+'_'+str(pop)+'.data')
        return
    ## \fn plot_data(self, ax, basename, prof, pop, gp)
    # plot data as blue shaded region
    # @param ax axis object
    # @param basename string
    # @param prof string of profile
    # @param pop which population to analyze: 0 (all), 1, 2, ...
    # @param gp global parameters

    def plot_labels(self, ax, prof, pop, gp):
        if prof=='chi2':
            ax.set_xlabel('$\\log_{10}\\chi^2$')
            ax.set_ylabel('frequency')
        else:
            ax.set_xlabel('$R\\quad[\\rm{pc}]$')
        if prof == 'rho':
            ax.set_ylabel('$\\rho\\quad[\\rm{M}_\\odot/\\rm{pc}^3]$')
        elif prof == 'M':
            ax.set_ylabel('$M(r)$')
        elif prof == 'nr':
            ax.set_ylabel('$n(r)$')
            self.broaden_lim('nr', 0, -1, 5)
        elif prof == 'beta':
            ax.set_ylabel('$\\beta_'+str(pop)+'$')
            self.broaden_lim('beta', pop, -1.05, 1.05)
        elif prof == 'betastar':
            ax.set_ylabel('$\\beta^*_'+str(pop)+'$')
            self.broaden_lim('betastar', pop, -1.05, 1.05)
        elif prof == 'Sig' and pop > 0:
            ax.set_ylabel('$\\Sigma_'+str(pop)+'\\quad[\\rm{M}_\\odot/\\rm{pc}^2]$')
        elif prof == 'Sig' and pop == 0:
            ax.set_ylabel('$\\Sigma^*\\quad[\\rm{M}_\\odot/\\rm{pc}^2]$')
        elif prof == 'nu' and pop > 0:
            ax.set_ylabel('$\\nu_'+str(pop)+'\\quad[\\rm{M}_\\odot/\\rm{pc}^3]$')
        elif prof == 'nu' and pop == 0:
            ax.set_ylabel('$\\rho^*\\quad[\\rm{M}_\\odot/\\rm{pc}^3]$')
        elif prof == 'nrnu':
            ax.set_ylabel('$n_{\\nu,'+str(pop)+'}(r)$')
        elif prof == 'sig':
            ax.set_ylabel('$\\sigma_{\\rm{LOS},'+str(pop)+'}\\quad[\\rm{km}/\\rm{s}]$')
        return
    ## \fn plot_labels(self, ax, prof, pop, gp)
    # draw x and y labels
    # @param ax axis object to use
    # @param prof string of profile to look at
    # @param pop which population to analyze: 0 (all), 1, 2, ...
    # @param gp global parameters

    def plot_Xscale_3D(self, ax, gp):
        rmin = np.log10(min(gp.xipol))
        rmax = np.log10(max(gp.xipol))
        gp.xfine = np.logspace(rmin, rmax, gp.nfine)
        if gp.investigate == 'walk':
            if gp.pops == 1:
                rhodm, nu1 = ga.rho_walk(gp.xepol, gp)
            else:
                rhodm, nu1, nu2 = ga.rho_walk(gp.xepol, gp)
        elif gp.investigate == 'gaia':
            rhodm, nu1 = ga.rho_gaia(gp.xepol, gp)
        for pop in range(gp.pops):
            # use our models
            nuprof = self.Mmedi.get_prof('nu', pop+1)
            tck = splrep(gp.xepol, nuprof)
            nuproffine = splev(gp.xfine, tck)
            if gp.investigate == 'walk' or gp.investigate == 'gaia':
                # or rather use analytic values, where available
                if pop == 0:
                    nuprof = nu1
                elif pop == 1:
                    nuprof = nu2
            if gp.geom == 'sphere':
                Mprof = gip.rho_SUM_Mr(gp.xfine, nuproffine)
                Mmax = max(Mprof) # Mprof[-1]
                ihalf = -1
                for kk in range(len(Mprof)):
                    # half-light radius (3D) is where mass is more than half
                    # ihalf gives the iindex of where this happens
                    if Mprof[kk] >= Mmax/2 and ihalf < 0:
                        xx = (gp.xfine[kk-1]+gp.xfine[kk])/2
                        print('rhalf = ', xx, ' pc')
                        ax.axvline(xx, color='green', lw=0.5, alpha=0.7)
                        ihalf = kk
    ## \fn plot_Xscale_3D(ax, gp)
    # plot 3D half-light radii, based on median nu model
    # @param ax axis object
    # @param gp global parameters

    def plot_full_distro(self, ax, prof, pop, gp):
        x = self.x0
        y = self.sort_prof(prof, pop, gp) # gives [Nmodels, Nbin] shape matrix
        Nvertbin = np.sqrt(len(y[:,0]))
        xbins = np.hstack([self.binmin, self.binmax[-1]])
        ybins = np.linspace(min(y[:,:]), max(y[:,:]), num=Nvertbin)

        H, xedges, yedges = np.histogram2d(x, y, bins=[xbins, ybins])
        fig, ax = plt.subplots(figsize=(8, 8), dpi=80, facecolor='w', edgecolor='k')

        ax.set_xlim([self.x[0],self.x0[-1]])

        extent = [min(xbins), max(xbins), min(ybins), max(ybins)]
        im = ax.imshow(H.T, extent=extent, aspect='auto', interpolation='nearest', cmap=plt.cm.binary) #plt.cm.Blues)
        fig.colorbar(im)
        if prof == 'Sig' or prof == 'sig':
            for pop in range(gp.pops):
                ax.axvline(gp.Xscale[pop+1], color='blue', lw=0.5) # [pc]
        else:
            self.plot_Xscale_3D(ax, gp)

        return
    ## \fn plot_full_distro(self, ax, prof, pop, gp)
    # plot filled region with grayscale propto # models in log bin
    # @param ax axis object to plot into
    # @param prof string
    # @param pop int
    # @param gp

    def fill_nice(self, ax, prof, pop, gp):
        M99lo = self.M99lo.get_prof(prof, pop)
        M95lo = self.M95lo.get_prof(prof, pop)
        M68lo = self.M68lo.get_prof(prof, pop)
        Mmedi = self.Mmedi.get_prof(prof, pop)
        r0    = gp.xepol
        M68hi = self.M68hi.get_prof(prof, pop)
        M95hi = self.M95hi.get_prof(prof, pop)
        M99hi = self.M99hi.get_prof(prof, pop)
        #print('lengths:', len(r0), len(M95lo), len(M95hi))
        #print('min max M95: ', min(M95lo), max(M95lo), prof)
        ax.fill_between(r0, M99lo, M99hi, color='black', alpha=0.1, lw=0.1)
        ax.plot(r0, M99lo, color='black', lw=0.4)
        ax.plot(r0, M99hi, color='black', lw=0.3)
        ax.fill_between(r0, M95lo, M95hi, color='black', alpha=0.2, lw=0.1)
        ax.plot(r0, M95lo, color='black', lw=0.4)
        ax.plot(r0, M95hi, color='black', lw=0.3)
        ax.fill_between(r0, M68lo, M68hi, color='black', alpha=0.4, lw=0.1)
        ax.plot(r0, M68lo, color='black', lw=0.4)
        ax.plot(r0, M68hi, color='black', lw=0.3)
        ax.plot(r0, Mmedi, 'r', lw=1)
        if prof == 'Sig' or prof == 'sig':
            for pop in range(gp.pops):
                ax.axvline(gp.Xscale[pop+1], color='blue', lw=0.5) # [pc]
        else:
            self.plot_Xscale_3D(ax, gp)
        ax.set_xlim([r0[0], r0[-1]])
        if prof == 'beta' or prof == 'betastar':
            self.broaden_lim('beta', pop, -1, 1)
            self.broaden_lim('betastar', pop, -1, 1)
        elif prof == 'nr':
            self.broaden_lim('nr', pop, -0.5, 5)
        elif prof == 'nrnu':
            self.broaden_lim('nrnu', pop, 0., max(M95hi))
        elif prof == 'M':
            self.broaden_lim('M', 0, min(M68lo), max(M68hi))
        elif prof == 'rho':
            self.broaden_lim('rho', 0, min(M68lo), max(M68hi))
        elif prof == 'nu':
            self.broaden_lim('nu', pop, min(M68lo), max(M68hi))
        return
    ## \fn fill_nice(self, ax, prof, pop, gp)
    # plot filled region for 1sigma and 2sigma confidence interval
    # @param ax axis object to plot into
    # @param prof string
    # @param pop int
    # @param gp

    def plot_profile(self, basename, prof, pop, gp):
        gh.LOG(1, 'prof '+str(prof)+', pop '+str(pop)+', run '+basename)
        fig = plt.figure()
        ax  = fig.add_subplot(111)
        if prof != 'chi2':
            ax.set_xscale('log')
        if prof == 'rho' or prof == 'J' or prof == 'Sig' or\
           prof == 'M' or prof == 'nu':
            ax.set_yscale('log')
        self.plot_labels(ax, prof, pop, gp)
        if len(self.profs)>0:
            if prof == 'chi2':
                goodchi = []
                for k in range(len(self.profs)):
                    # do include all chi^2 values for plot
                    goodchi.append(self.chis[k])
                print('plotting profile chi for '+str(len(goodchi))+' models')
                bins, edges = np.histogram(np.log10(goodchi), range=[-2,6], \
                                           bins=max(6,np.sqrt(len(goodchi))),\
                                           density=True)
                ax.step(edges[1:], bins, where='pre')
                plt.draw()
                self.write_chi2(basename, edges, bins)
                fig.savefig(basename+'output/prof_chi2_0.pdf')
                return

            self.fill_nice(ax, prof, pop, gp)
            # TODO: replace above with full distribution plot (has bugs,
            # File "programs/plotting/gi_collection.py", line 466, in plot_full_distro
            # ybins = np.linspace(min(y[:,:]), max(y[:,:]), num=Nvertbin)
            # ValueError: The truth value of an array with more than one element is ambiguous. Use a.any() or a.all()
            #self.plot_full_distro(ax, prof, pop, gp)

            self.plot_N_samples(ax, prof, pop)
            if prof == 'Sig' or prof == 'sig':
                self.plot_data(ax, basename, prof, pop, gp)

            if (gp.investigate == 'gaia' or gp.investigate=='triax') and prof != 'sig' or (gp.investigate=='walk' and prof!='sig' and prof!='nu' and prof!='Sig'):
                r0 = self.analytic.x0
                y0 = self.analytic.get_prof(prof, pop)
                self.broaden_lim(prof, pop, min(y0), max(y0))
                ax.plot(r0, y0, 'b--', lw=2)
            ax.set_ylim(self.ranges[prof+str(pop)])
            plt.draw()
        else:
            gh.LOG(1, 'empty self.profs')
        fig.savefig(basename+'output/pdf/prof_'+prof+'_'+str(pop)+'.pdf')
        return 1
    ## \fn plot_profile(self, basename, prof, pop, gp)
    # plot single profile
    # @param basename
    # @param prof string of profile to look at
    # @param pop population number
    # @param gp global parameters

    def __repr__(self):
        return "Profile Collection with "+str(len(self.profs))+" Profiles"
    ## \fn __repr__(self)
    # string representation for ipython
