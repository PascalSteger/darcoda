#!/usr/bin/env python3

## @file
# collect profiles and perform actions on them

# (c) 2014 ETHZ Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import numpy.random as npr
import pdb

import matplotlib.pyplot as plt
plt.ioff()

from gl_class_profiles import Profiles
import gl_output as go
import gl_helper as gh
import gl_analytic as ga
import gl_project as glp
import gl_physics as phys

USE_ALL = False

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
    elif prof == 'sig' or prof == 'sig_vec':
        unit = '[km/s]'
    elif prof == 'sigz2_vec':
        unit = '[km^2/s^2]'
    elif prof == 'nu_vec':
        unit = '[stars/kpc^3]'
    elif prof == 'rho_DM_vec':
        unit = '[Msum/kpc^3]'
    else:
        unit = 'a.u.'
    return unit
## \fn unit(prof)
# return string with units for profile
# @param prof string for profile


class ProfileCollection():
    def __init__(self, pops, nipol):
        self.chis = []
        self.goodchi = []  # will contain all chi^2 for models in self.subset
        self.goodprof = [] # will contain the corresponding profiles
        self.profs = []
        self.subset = [0., np.inf]
        self.x0 = np.array([])
        self.binmin = np.array([])
        self.binmax = np.array([])
        self.Mmin  = Profiles(pops, nipol)
        self.M95lo = Profiles(pops, nipol)
        self.M68lo = Profiles(pops, nipol)
        self.Mmedi = Profiles(pops, nipol)
        self.M68hi = Profiles(pops, nipol)
        self.M95hi = Profiles(pops, nipol)
        self.Mmax  = Profiles(pops, nipol)
        self.analytic = Profiles(pops, 100)
    ## \fn __init__(self, pops, nipol)
    # constructor
    # @param pops number of populations
    # @param nipol number of bins

    def add(self, prof):
        if np.isnan(prof.chi2):
            return
        self.chis.append(prof.chi2)
        self.profs.append(prof)
        return
    ## \fn add(self, prof)
    # add a profile
    # @param prof Profile to add



    def cut_subset(self):
        self.chis = np.array(self.chis) # better to work with numpy arrays for subsections
        minchi = min(self.chis)
        if USE_ALL:
            maxchi = max(self.chis)
        else:
            counts, bins = np.histogram(np.log10(self.chis), bins=np.sqrt(len(self.chis)))
            mincounts = -1
            found_max = False; found_min = False
            for k in range(len(counts)):
                # flag if we found the first max peak in the chi2 distribution
                if counts[k] < mincounts and not found_max:
                    found_max = True

                # if still on rising part, store actual entry as new maximum
                # and set starting minimum chi2 at highest value, to be changed later on
                if counts[k] >= mincounts and not found_max:
                    mincounts = counts[k]
                    maxcounts = counts[k]
                # if on decreasing branch, continue to set down border
                if counts[k] < maxcounts and found_max and not found_min:
                    mincounts = counts[k]
                if counts[k] >= maxcounts and found_max:
                    found_min = True
                    break
            maxchi = 10**(bins[k+1])
        self.subset = [minchi, 10*minchi]# maxchi]
    ## \fn cut_subset(self)
    # set subset to [0, 10*min(chi)] (or 30* minchi, or any value wished)


    def set_x0(self, x0):
        self.x0 = x0
        self.Mmin.x0  = x0
        self.M95lo.x0 = x0
        self.M68lo.x0 = x0
        self.Mmedi.x0 = x0
        self.M68hi.x0 = x0
        self.M95hi.x0 = x0
        self.Mmax.x0  = x0
    ## \fn set_x0(self, x0)
    # set radii for 1sigma, 2sigma profiles
    # @param x0 radii in pc


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
        self.M95lo.set_prof(prof, tmp[ll*0.05], pop, gp)
        self.M68lo.set_prof(prof, tmp[ll*0.32], pop, gp)
        self.Mmedi.set_prof(prof, tmp[ll/2],    pop, gp)
        self.M68hi.set_prof(prof, tmp[ll*0.68], pop, gp)
        self.M95hi.set_prof(prof, tmp[ll*0.95], pop, gp)
        self.Mmax.set_prof(prof,  tmp[-1],      pop, gp)
        return tmp
    ## \fn sort_prof(self, prof, pop, gp)
    # sort the list of prof-profiles, and store the {1,2}sigma, min, medi, max in the appropriate place
    # @param prof profile identifier, 'rho', 'beta', ...
    # @param pop population identifier
    # @param gp global parameters
    # @return tmp profiles sorted binwise

    def calculate_M(self, gp):
        if len(self.profs)>0:
            for i in range(len(self.profs)):
                Mprof = glp.rho_SUM_Mr(gp.xipol, self.profs[i].get_prof('rho', 1))
                self.profs[i].set_prof('M', Mprof, 1, gp)
                if gp.pops == 2:
                    Mprof = glp.rho_SUM_Mr(gp.xipol, self.profs[i].get_prof('rho', 2))
                    self.profs[i].set_prof('M', Mprof, 2, gp)
        else:
            gh.LOG(1, 'len(self.profs)==0, did not calculate self.profs.M')
        return
    ## \fn calculate_M(self, gp)
    # calculate M profiles from rho, as this has not been done prior to pc2.save
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
    ## \fn sort_profiles(self, gp)
    # sort all profiles
    # @param gp global parameters

    def sort_profiles_disc(self, gp):
        self.sort_prof('nu_vec', 0, gp)
        #self.sort_prof('sig_vec', 0, gp) #Old plots
        self.sort_prof('sigz2_vec', 0, gp)

        self.sort_prof('kz_rho_DM_vec', 0, gp)
        self.sort_prof('kz_nu_vec', 0, gp)

        self.sort_prof('rho_DM_vec', 0, gp)
        self.sort_prof('Sig_DM_vec', 0, gp)

        if gp.baryonmodel not in ['none', 'simplenu_baryon']:
            return

        self.sort_prof('rho_baryon_vec', 0, gp)
        self.sort_prof('Sig_baryon_vec', 0, gp)

        self.sort_prof('rho_total_vec', 0, gp)
        self.sort_prof('Sig_total_vec', 0, gp)


    def set_analytic(self, x0, gp):
        r0 = x0 # [pc], spherical case
        self.analytic.x0 = r0
        anbeta = []; annu = [];  anSig = []
        if gp.investigate == 'gaia':
            anrho = ga.rho_gaia(r0, gp)[0]
            anM = glp.rho_SUM_Mr(r0, anrho)
            annr = ga.nr3Dtot_gaia(r0, gp)
            tmp_annu = ga.rho_gaia(r0, gp)[1]
            annu.append( tmp_annu )
            anSig.append( glp.rho_INT_Sig(r0, tmp_annu, gp) )
            for pop in np.arange(1, gp.pops+1):
                beta = ga.beta_gaia(r0, gp)[pop]
                anbeta.append(beta)
                nu = ga.rho_gaia(r0,gp)[pop]
                annu.append(nu)
                anSig.append(glp.rho_INT_Sig(r0, nu, gp))
        elif gp.investigate == 'walk':
            anrho = ga.rho_walk(r0, gp)[0]
            anM = glp.rho_SUM_Mr(r0, anrho)
            annr = ga.nr3Dtot_deriv_walk(r0, gp) # TODO too high in case of core
            tmp_annu = ga.rho_walk(r0, gp)[1]
            annu.append( tmp_annu )
            anSig.append( glp.rho_INT_Sig(r0, tmp_annu, gp) )
            for pop in np.arange(1, gp.pops+1):
                beta = ga.beta_walk(r0, gp)[pop]
                anbeta.append(beta)
                nu = ga.rho_walk(r0, gp)[pop]
                annu.append(nu)
                anSig.append(glp.rho_INT_Sig(r0, nu, gp))
        elif gp.investigate == 'triax':
            anrho = ga.rho_triax(r0, gp)[0]
            anM = glp.rho_SUM_Mr(r0, anrho)
            annr = ga.nr3Dtot_deriv_triax(r0, gp)
            tmp_annu = ga.rho_triax(r0, gp)[1]
            annu.append(tmp_annu)
            anSig.append( glp.rho_INT_Sig(r0, tmp_annu, gp))
            for pop in np.arange(1, gp.pops+1):
                beta = ga.beta_triax(r0, gp)[pop]
                anbeta.append(beta)
                nu = ga.rho_triax(r0, gp)[pop]
                annu.append(nu)
                anSig.append( glp.rho_INT_Sig(r0, nu, gp))
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
            self.analytic.set_prof('sig', -np.ones(len(r0)), pop, gp) # TODO: find analytic profile
        return
    ## \fn set_analytic(x0, gp)
    # set analytic curves (later shown in blue)
    # @param x0 radius in [pc]
    # @param gp global parameters


    def write_prof(self, basename, prof, pop, gp):
        output = go.Output()
        uni = unit(prof)
        output.add('radius (models) [pc]', self.Mmedi.x0)
        output.add('M 95% CL low ' + uni, self.M95lo.get_prof(prof, pop))
        output.add('M 68% CL low ' + uni, self.M68lo.get_prof(prof, pop))
        output.add('M median '     + uni, self.Mmedi.get_prof(prof, pop))
        output.add('M 68% CL high '+ uni, self.M68hi.get_prof(prof, pop))
        output.add('M 95% CL high '+ uni, self.M95hi.get_prof(prof, pop))
        output.write(basename+'/output/ascii/prof_'+prof+'_'+str(pop)+'.ascii')
        if (gp.investigate =='walk' or gp.investigate=='gaia') \
           and (prof != 'Sig'):
            out_an = go.Output()
            out_an.add('radius [pc]', self.analytic.x0)
            out_an.add('analytic profile', self.analytic.get_prof(prof, pop))
            out_an.write(basename+'/output/analytic/prof_'+prof+'_'+str(pop)+'.analytic')
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
        output.write(basename+'/output/prof_chi2_0.ascii')
        return
    ## \fn write_chi2(self, basename, edges, bins)
    # write ascii file with chi2 information
    # @param basename directory string
    # @param edges array of edges
    # @param bins y-values


    def write_all(self, basename, gp):
        self.write_prof(basename, 'rho', 0, gp)
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

    def write_all_disc(self, basename, gp):
        self.write_prof(basename, 'nu_vec', 0, gp)
        #self.write_prof(basename, 'sig_vec', 0, gp)
        self.write_prof(basename, 'sigz2_vec', 0, gp)

        self.write_prof(basename, 'kz_rho_DM_vec', 0, gp)
        self.write_prof(basename, 'kz_nu_vec', 0, gp)

        self.write_prof(basename, 'rho_DM_vec', 0, gp)
        self.write_prof(basename, 'Sig_DM_vec', 0, gp)

        self.write_prof(basename, 'rho_baryon_vec', 0, gp)
        self.write_prof(basename, 'Sig_baryon_vec', 0, gp)

        self.write_prof(basename, 'rho_total_vec', 0, gp)
        self.write_prof(basename, 'Sig_total_vec', 0, gp)

    def plot_N_samples(self, ax, prof, pop):
        k=0
        while k<100:
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
        print('Plotting data for prof = ', prof)
        output = go.Output()
        r0 = self.Mmedi.x0 # [pc]
        output.add('radius (data) [pc]', r0)
        if prof == 'Sig':
            # get 2D data here
            # Rbin [Xscale], Binmin [Xscale], Binmax [Xscale], Sig(R)/Sig(0) [1], error [1]
            dum,dum,dum,Sigdat,Sigerr = np.transpose(np.loadtxt(gp.files.Sigfiles[pop], \
                                                              unpack=False, skiprows=1))
            #Sigdat *= gp.Sig0pc[pop] # [Msun/pc^2]
            #Sigerr *= gp.Sig0pc[pop] # [Msun/pc^2]
            #Signorm = gh.ipol_rhalf_log(gp.xipol, Sigdat, gp.Xscale[pop])
            output.add('data [Msun/pc^2]', Sigdat)
            output.add('error [Msun/pc^2]', Sigerr)
            output.add('data - error [Msun/pc^2]', Sigdat-Sigerr)
            output.add('data + error [Msun/pc^2]', Sigdat+Sigerr)
            ax.fill_between(r0, Sigdat-Sigerr, Sigdat+Sigerr, \
                            color='blue', alpha=0.3, lw=1)
            ax.set_ylim([min(Sigdat-Sigerr)/2., 2.*max(Sigdat+Sigerr)])
        elif prof == 'nu' or prof == 'nu_vec':
            # get 3D data here
            # Rbin [Xscale], Binmin [Xscale], Binmax [Xscale], nu(R)/nu(0) [1], error [1]
            dum,dum,dum,nudat,nuerr = np.transpose(np.loadtxt(gp.files.nufiles[pop], \
                                                              unpack=False, skiprows=1))
            #nudat *= gp.nu0pc[pop] # [Msun/pc^2]
            #nuerr *= gp.nu0pc[pop] # [Msun/pc^2]
            output.add('data [Msun/pc^3]', nudat)
            output.add('error [Msun/pc^3]', nuerr)
            output.add('data - error [Msun/pc^2]', nudat-nuerr)
            output.add('data + error [Msun/pc^2]', nudat+nuerr)
            ax.fill_between(r0, nudat-nuerr, nudat+nuerr, \
                            color='blue', alpha=0.3, lw=1)
            #ax.set_ylim([min(nudat-nuerr)/2., 2.*max(nudat+nuerr)]) #bodge
        elif prof == 'sig' or prof == 'sig_vec':
            DATA = np.transpose(np.loadtxt(gp.files.sigfiles[pop], \
                                           unpack=False, skiprows=1))
            sigdat = DATA[4-1] # [maxsiglosi]
            sigerr = DATA[5-1] # [maxsiglosi]
            #sigdat *= gp.maxsiglos[pop]  # [km/s]
            #sigerr *= gp.maxsiglos[pop]  # [km/s]
            output.add('data [km/s]', sigdat)
            output.add('error [km/s]', sigerr)
            output.add('data - error [km/s]', sigdat-sigerr)
            output.add('data + error [km/s]', sigdat+sigerr)
            ax.fill_between(r0, sigdat-sigerr, sigdat+sigerr, \
                            color='blue', alpha=0.3, lw=1)
            #ax.set_ylim([0., 2.*max(sigdat+sigerr)]) #bodge
        elif prof == 'sigz2_vec':
            DATA = np.transpose(np.loadtxt(gp.files.sigfiles[pop], \
                                           unpack=False, skiprows=1))
            sigz2dat = DATA[4-1] # [maxsiglosi]
            sigz2err = DATA[5-1] # [maxsiglosi]
            #sigdat *= gp.maxsiglos[pop]  # [km/s]
            #sigerr *= gp.maxsiglos[pop]  # [km/s]
            output.add('data [km^2/s^2]', sigz2dat)
            output.add('error [km^2/s^2]', sigz2err)
            output.add('data - error [km/s]', sigz2dat-sigz2err)
            output.add('data + error [km/s]', sigz2dat+sigz2err)
            ax.fill_between(r0, sigz2dat-sigz2err, sigz2dat+sigz2err, \
                            color='blue', alpha=0.3, lw=1)
            #ax.set_ylim([0., 2.*max(sigdat+sigerr)]) #bodge
        output.write(basename+'/output/data/prof_'+prof+'_'+str(pop)+'.data')
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
        elif gp.geom == 'disc' and prof != 'chi2':
            ax.set_xlabel('$z\\quad[\\rm{kpc}]$')
        elif gp.geom == 'sphere' and prof != 'chi2':
            ax.set_xlabel('$R\\quad[\\rm{pc}]$')

        if prof == 'rho':
            ax.set_ylabel('$\\rho\\quad[\\rm{M}_\\odot/\\rm{pc}^3]$')
        elif prof == 'M':
            ax.set_ylabel('$M(r)$')
            ax.set_ylim([1e5, 1e11])
        elif prof == 'nr':
            ax.set_ylabel('$n(r)$')
            ax.set_ylim([-1,5])
        elif prof == 'beta':
            ax.set_ylabel('$\\beta_'+str(pop)+'$')
            ax.set_ylim([-1.05,1.05])
        elif prof == 'betastar':
            ax.set_ylabel('$\\beta^*_'+str(pop)+'$')
            ax.set_ylim([-1.05,1.05])
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
        # Disc cases
        elif prof == 'nu_vec':
            ax.set_ylabel('$\\nu_{\\rm{Tr},'+str(pop)+'}\\quad[\\rm{stars}/\\rm{kpc}^3]$')
            ax.set_ylim(1.0E3, 4.0E4)
        elif prof == 'sig_vec':
            ax.set_ylabel('$\\sigma_{z}\\quad[\\rm{km}/\\rm{s}]$')
            ax.set_ylim(10., 40.)
        elif prof == 'sigz2_vec':
            ax.set_ylabel('$\\sigma_{z}^2\\quad[\\rm{km}^2/\\rm{s}^2]$')
            #ax.set_ylim(100., 1600.)

        elif prof == 'kz_nu_vec':
            ax.set_ylabel('$k(z)_{\\nu}$')
            ax.set_ylim(0., 6.)
        elif prof == 'kz_rho_DM_vec':
            ax.set_ylabel('$k(z)_{\\rho, \\rm{DM}}$')
            ax.set_ylim(-10., 10.)

        elif prof == 'rho_DM_vec':
            ax.set_ylabel('$\\rho_{\\rm{DM}}\\quad[\\rm{M}_\\odot/\\rm{kpc}^3]$')
            ax.set_ylim(1E3, 1E9)
        elif prof == 'Sig_DM_vec':
            ax.set_ylabel('$\\Sigma_{\\rm{DM}} \\quad[\\rm{M}_\\odot/\\rm{kpc}^2]$')
            ax.set_ylim(0,1.0E8)

        elif prof == 'rho_baryon_vec':
            ax.set_ylabel('$\\rho_{\\rm{baryon}}\\quad[\\rm{M}_\\odot/\\rm{kpc}^3]$')
            ax.set_ylim(1E3, 1E9)
        elif prof == 'Sig_baryon_vec':
            ax.set_ylabel('$\\Sigma_{\\rm{baryon}} \\quad[\\rm{M}_\\odot/\\rm{kpc}^2]$')
            ax.set_ylim(0,1.0E8)

        elif prof == 'rho_total_vec':
            ax.set_ylabel('$\\rho_{\\rm{total}}\\quad[\\rm{M}_\\odot/\\rm{kpc}^3]$')
            ax.set_ylim(1E3, 1E9)
        elif prof == 'Sig_total_vec':
            ax.set_ylabel('$\\Sigma_{\\rm{total}} \\quad[\\rm{M}_\\odot/\\rm{kpc}^2]$')
            ax.set_ylim(0,1.0E8)



#
#
#        elif prof == 'sig' or prof == 'sig_vec':
#        unit = '[km/s]'
#    elif prof == 'nu_vec':
#        unit = '[stars/kpc^3]'
#    elif prof == 'rho_DM_vec':
#        unit = '[Msum/kpc^3]'

        return
    ## \fn plot_labels(self, ax, prof, pop, gp)
    # draw x and y labels
    # @param ax axis object to use
    # @param prof string of profile to look at
    # @param pop which population to analyze: 0 (all), 1, 2, ...
    # @param gp global parameters


    def plot_Xscale_3D(self, ax, gp):
        if gp.investigate == 'walk':
            if gp.pops == 1:
                rhodm, nu1 = ga.rho_walk(gp.xipol, gp)
            else:
                rhodm, nu1, nu2 = ga.rho_walk(gp.xipol, gp)
        elif gp.investigate == 'gaia':
            rhodm, nu1 = ga.rho_gaia(gp.xipol, gp)
        for pop in range(gp.pops):
            # use our models
            nuprof = self.Mmedi.get_prof('nu', pop+1)
            if gp.investigate == 'walk' or gp.investigate == 'gaia':
                # or rather use analytic values, where available
                if pop == 0:
                    nuprof = nu1
                elif pop == 1:
                    nuprof = nu2
            if gp.geom == 'sphere':
                Mprof = glp.rho_SUM_Mr(gp.xipol, nuprof)
                Mmax = max(Mprof) # Mprof[-1]
                ihalf = -1
                for kk in range(len(Mprof)):
                    # half-light radius (3D) is where mass is more than half
                    # ihalf gives the iindex of where this happens
                    if Mprof[kk] >= Mmax/2 and ihalf <0:
                        xx = (gp.xipol[kk-1]+gp.xipol[kk])/2
                        ax.axvline(xx, color='green', lw=0.5, alpha=0.7)
                        ihalf = kk
    ## \fn plot_Xscale_3D(ax, gp)
    # plot 3D half-light radii, based on median nu model
    # @param ax axis object
    # @param gp global parameters



    def fill_nice(self, ax, prof, pop, gp):
        M95lo = self.M95lo.get_prof(prof, pop)
        M68lo = self.M68lo.get_prof(prof, pop)
        Mmedi = self.Mmedi.get_prof(prof, pop)
        r0    = self.Mmedi.x0
        M68hi = self.M68hi.get_prof(prof, pop)
        M95hi = self.M95hi.get_prof(prof, pop)
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
        ax.set_xlim([r0[0]/2,r0[-1]*1.0]) #HS r0[-1] previously multiplied by 1.5
        #if prof == 'kz_nu_vec':
        #    ax.set_ylim(-1., 6.)
        # if gp.pops==1:
        #     ax.set_xlim([r0[0], 5*gp.Xscale[1]])
        # else:
        #     ax.set_xlim([r0[0], 5*max(gp.Xscale[1], gp.Xscale[2])])
        return
    ## \fn fill_nice(self, ax, prof, pop, gp)
    # plot filled region for 1sigma and 2sigma confidence interval
    # @param ax axis object to plot into
    # @param prof string
    # @param pop int
    # @param gp


    def plot_full_distro(self, ax, prof, pop, gp):
        x = self.x0
        y = self.sort_prof(prof, pop, gp) # gives [Nmodels, Nbin] shape matrix
        Nvertbin = np.sqrt(len(y[:,0]))
        xbins = np.hstack([self.binmin, self.binmax[-1]])
        ybins = np.linspace(min(y[:,:]), max(y[:,:]), num=Nvertbin)

        H, xedges, yedges = np.histogram2d(x, y, bins=[xbins, ybins])
        fig, ax = plt.subplots(figsize=(8, 8), dpi=80, facecolor='w', edgecolor='k')

        ax.set_xlim([self.x[0]/2,self.x0[-1]*1.5])

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


    def plot_profile(self, basename, prof, pop, gp):
        gh.LOG(1, 'plotting profile '+str(prof)+' for pop '+str(pop)+' in run '+basename)

        fig = plt.figure()
        ax  = fig.add_subplot(111)


        if prof != 'chi2':
            ax.set_xscale('log')

        if prof == 'rho' or prof == 'Sig' or\
           prof == 'M' or prof == 'nu':
            ax.set_yscale('log')

        if prof == 'nu_vec' or prof == 'rho_DM_vec' or prof == 'rho_baryon_vec' or prof == 'rho_total_vec':
            ax.set_xscale('linear')
            ax.set_yscale('log')

        if prof == 'kz_rho_DM_vec' or prof == 'kz_nu_vec' or prof == 'sig_vec' or prof == 'sigz2_vec' or prof == 'Sig_DM_vec'  or prof == 'Sig_baryon_vec'  or prof == 'Sig_total_vec':
            ax.set_xscale('linear')
            ax.set_yscale('linear')

        self.plot_labels(ax, prof, pop, gp)

        if len(self.profs)>0:
            if prof == 'chi2':
                goodchi = []
                for k in range(len(self.profs)):
                    # do include all chi^2 values for plot
                    goodchi.append(self.chis[k])
                print('plotting profile chi for '+str(len(goodchi))+' models')
                print('min, max, maxsubset found chi2: ', min(self.chis), max(self.chis), self.subset[1])
                bins, edges = np.histogram(np.log10(goodchi), range=[-2,6], \
                                           bins=max(6,np.sqrt(len(goodchi))),\
                                           density=True)
                ax.step(edges[1:], bins, where='pre')
                plt.draw()
                self.write_chi2(basename, edges, bins)
                fig.savefig(basename+'/output/prof_chi2_0.pdf')
                return

            self.fill_nice(ax, prof, pop, gp)
            self.plot_N_samples(ax, prof, pop)
            self.plot_bins(ax, prof, gp)
            if prof == 'Sig' or prof=='nu' or prof == 'sig' or prof == 'nu_vec' or prof == 'sig_vec' or prof == 'sigz2_vec':
                self.plot_data(ax, basename, prof, pop, gp)

            if gp.investigate == 'simplenu':
                self.plot_model_simplenu(ax, basename, prof, gp)

            if (gp.investigate == 'walk' or gp.investigate == 'gaia') \
               and (prof != 'Sig'):
                r0 = self.analytic.x0
                y0 = self.analytic.get_prof(prof, pop)
                ax.plot(r0, y0, 'b--', lw=2)

            #ax.set_xlim(0, gp.z_bincenters[-1]) #bodge
            ax.set_xlim(0, np.round(gp.z_binmaxs[-1], 1))

            plt.draw()
        else:
            gh.LOG(1, 'empty self.profs')
        fig.savefig(basename+'/output/pdf/prof_'+prof+'_'+str(pop)+'.pdf')
        return
    ## \fn plot_profile(self, basename, prof, pop, gp)
    # plot single profile
    # @param basename
    # @param prof string of profile to look at
    # @param pop population number
    # @param gp global parameters


    def plot_from_ascii_chi(self, ax, chidata):
        edges = chidata[:,0]
        bins  = chidata[:,1]
        ax.step(edges, bins, where='pre')
    ## \fn plot_from_ascii_chi(self, ax, chidata)
    # plot chi ascii files
    # @param ax axis profile
    # @param chidata data


    def plot_from_ascii(self, ax, data, data_analytic, gp):
        r0 = data[:,0]
        M95lo = data[:,1]
        M68lo = data[:,2]
        Mmedi = data[:,3]
        M68hi = data[:,4]
        M95hi = data[:,5]
        ax.fill_between(r0, M95lo, M95hi, color='black', alpha=0.2, lw=0.1)
        ax.plot(r0, M95lo, color='black', lw=0.4)
        ax.plot(r0, M95hi, color='black', lw=0.3)
        ax.fill_between(r0, M68lo, M68hi, color='black', alpha=0.4, lw=0.1)
        ax.plot(r0, M68lo, color='black', lw=0.4)
        ax.plot(r0, M68hi, color='black', lw=0.3)
        ax.plot(r0, Mmedi, 'r', lw=1)
        for pop in range(gp.pops):
            ax.axvline(gp.Xscale[pop+1], color='gray', lw=0.5) # [pc]

        ranalytic = data_analytic[:,0]
        analytic  = data_analytic[:,1]
        ax.plot(ranalytic, analytic, lw=2, color='blue')
        return
    ## \fn plot_from_ascii(self, ax, data, data_analytic, gp)
    # plot a single plot from previously calculated data
    # @param ax axis object to plot into
    # @param data array of (M95lo, M68lo, Mmedi, M68hi, M95hi) for the profile in question
    # @param data_analytic analytic profile
    # @param gp global parameters


    def plot_model_simplenu(self, ax, basename, prof, gp):
        G1 = 4.299e-6   # Newton's constant in (km)^2*kpc/(Msun*s^2)
        K=1500.
        F=267.65
        D=0.18
        z0=0.4
        normC = 22.85
        nu0 = 1.
        Ntr = 10000.
        zmax= 1.2

        nu0 = Ntr/(z0*(1-np.exp(-zmax/z0)))

        zvec=np.linspace(0, zmax, 100)

        nuvec = nu0*np.exp(-zvec/z0)

        Kzvec_total = -((K*zvec)/(np.sqrt(zvec**2 + D**2)) + 2.*F*zvec)
        Kzvec_baryon = -((K*zvec)/(np.sqrt(zvec**2 + D**2)))
        Kzvec_DM = -(2.*F*zvec)


        Sigma_z_total = (1000.**2)*abs(Kzvec_total)/(2*np.pi*4.299) #Msun kpc^-1
        Sigma_z_baryon = (1000.**2)*abs(Kzvec_baryon)/(2*np.pi*4.299) #Msun kpc^-1
        Sigma_z_DM = (1000.**2)*abs(Kzvec_DM)/(2*np.pi*4.299) #Msun kpc^-1

        rho_z_total = (1/(4*np.pi*G1)) * abs((K*(D**2)/((D**2 + zvec**2)**(1.5))) + 2.*F)
        rho_z_baryon = (1/(4*np.pi*G1)) * abs((K*(D**2)/((D**2 + zvec**2)**(1.5))))
        rho_z_DM = (1/(4*np.pi*G1)) * abs(2.*F) * np.ones(len(zvec))

        k_z_rho_total = (3*(D**2)*K*zvec) / ( ((D**2 + zvec**2)**2.5) * ((K*(D**2))/((D**2 + zvec**2)**1.5) + 2*F))
        k_z_rho_baryon = (3*(D**2)*K*zvec) / ( ((D**2 + zvec**2)**2.5) * ((K*(D**2))/((D**2 + zvec**2)**1.5))) # not currently used

        # Backwards compatibility: if using old data, then all mass is in DM
        if gp.baryonmodel not in ['simplenu_baryon']:
            gh.LOG(1, 'Simplenu Analytic: No baryon model, all mass is in DM.')
            Sigma_z_DM = Sigma_z_total
            rho_z_DM = rho_z_total

        if prof == 'nu_vec':
            ax.plot(zvec, nuvec, 'g-', alpha=0.5)

        elif prof == 'Sig_total_vec':
            ax.plot(zvec, Sigma_z_total, 'g-', alpha=0.5)
        elif prof == 'Sig_baryon_vec':
            ax.plot(zvec, Sigma_z_baryon, 'g-', alpha=0.5)
        elif prof == 'Sig_DM_vec':
            ax.plot(zvec, Sigma_z_DM, 'g-', alpha=0.5)

        elif prof == 'rho_total_vec':
            ax.plot(zvec, rho_z_total, 'g-', alpha = 0.5)
        elif prof == 'rho_baryon_vec':

            #Prior range on baryons
            #pdb.set_trace()
            def rho_baryon(z, K, D):
                return (1/(4*np.pi*G1)) * abs((K*(D**2)/((D**2 + zval**2)**(1.5))))

            rho_z_baryon_prior_max=[]
            rho_z_baryon_prior_min=[]

            K_vec = np.linspace(gp.simplenu_baryon_K_min, gp.simplenu_baryon_K_max, 1000)
            D_vec = np.linspace(gp.simplenu_baryon_D_min+1.E-6, gp.simplenu_baryon_D_max, 1000)
            K_vec, D_vec = np.meshgrid(K_vec, D_vec)

            for zval in zvec:
                rho_grid = rho_baryon(zval, K_vec, D_vec)
                rho_z_baryon_prior_max.append(rho_grid.max())
                rho_z_baryon_prior_min.append(rho_grid.min())

            ax.plot(zvec, rho_z_baryon, 'g-', alpha = 0.5)
            ax.fill_between(zvec, rho_z_baryon_prior_min, rho_z_baryon_prior_max, color='r', alpha=0.1, lw=1)
            #ax.plot(zvec, rho_z_baryon_prior_max,'g-', alpha = 0.5, linewidth=1)
            #ax.plot(zvec, rho_z_baryon_prior_min,'g-', alpha = 0.5, linewidth=1)


        elif prof == 'rho_DM_vec':
            ax.plot(zvec, rho_z_DM, 'g-', alpha = 0.5)

        # if all mass is described by DM, then plot kz_rho_DM_vec
        # if the simplenu_baryon model is used, then kz_rho_DM_vec should equal zero
        elif prof == 'kz_rho_DM_vec' and gp.baryonmodel in ['none', 'sim']: #i.e. if all mass is described by DM
            ax.plot(zvec, k_z_rho_total, 'g-', alpha = 0.5)


        elif prof == 'sigz2_vec':
            sigz2_analytic = phys.sigz2(zvec, Sigma_z_total, nuvec, (normC**2)*nuvec[0])
            ax.plot(zvec, sigz2_analytic, 'g-', alpha = 0.5)


        return


    def plot_bins(self, ax, prof, gp):
        [ax.axvline(x, color='c', linewidth=0.1) for x in gp.z_bincenters]
        return





    def __repr__(self):
        return "Profile Collection with "+str(len(self.profs))+" Profiles"
    ## \fn __repr__(self)
    # string representation for ipython
