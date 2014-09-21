#!/usr/bin/env ipython

## @file
# collect profiles and perform actions on them

# (c) 2014 ETHZ Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import numpy.random as npr
import pdb, sys
import matplotlib
#matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

plt.ioff()

from gl_class_profiles import Profiles
import gl_output as go
import gl_helper as gh
import gl_analytic as ga
import gl_project as glp


def unit(prof):
    if prof == 'rho':
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
    # @param subset [chi_min, chi_max] for models to be included in analysis


    def add(self, prof):
        self.chis.append(prof.chi2)
        self.profs.append(prof)
    ## \fn add(self, prof)
    # add a profile
    # @param prof Profile to add


    def cut_subset(self):
        minchi = min(self.chis)
        maxchi = max(self.chis)
        self.subset = [0., 30.*minchi]
    ## \fn cut_subset(self)
    # set subset to [0, 10*min(chi)] (or 30* minchi, or any value wished)


    def set_x0(self, x0):
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
        #    norm = gh.ipol_rhalf_log(gp.xepol, tmp[ll/2], gp.Rscale[0])

        self.Mmin.set_prof(prof,  tmp[0],       pop, gp)
        self.M95lo.set_prof(prof, tmp[ll*0.05], pop, gp)
        self.M68lo.set_prof(prof, tmp[ll*0.32], pop, gp)
        self.Mmedi.set_prof(prof, tmp[ll/2],    pop, gp)
        self.M68hi.set_prof(prof, tmp[ll*0.68], pop, gp)
        self.M95hi.set_prof(prof, tmp[ll*0.95], pop, gp)
        self.Mmax.set_prof(prof,  tmp[-1],      pop, gp)

    ## \fn sort_prof(self, prof, pop, gp)
    # sort the list of prof-profiles, and store the {1,2}sigma, min, medi, max in the appropriate place
    # @param prof profile identifier, 'rho', 'beta', ...
    # @param pop population identifier
    # @param gp global parameters


    def calculate_M(self, gp):
        for i in range(len(self.profs)):
            Mprof = glp.rho_SUM_Mr(gp.xipol, self.profs[i].get_prof('rho', 1))
            self.profs[i].set_prof('M', Mprof, 1, gp)
            if gp.pops == 2:
                Mprof = glp.rho_SUM_Mr(gp.xipol, self.profs[i].get_prof('rho', 2))
                self.profs[i].set_prof('M', Mprof, 2, gp)
        return
    ## \fn calculate_M(self, gp)
    # calculate M profiles from rho, as this has not been done prior to pc2.save
    # @param gp global parameters


    def sort_profiles(self, gp):
        self.sort_prof('rho', 0, gp)
        self.sort_prof('M', 0, gp)
        if gp.geom=='sphere':
            self.sort_prof('nr', 0, gp)
        self.sort_prof('rhostar', 0, gp)
        self.sort_prof('betastar', 1, gp)
        self.sort_prof('beta', 1, gp)
        self.sort_prof('Sig', 1, gp)
        self.sort_prof('sig', 1, gp)
        if gp.pops == 2:
            self.sort_prof('betastar', 2, gp)
            self.sort_prof('beta', 2, gp)
            self.sort_prof('Sig', 2, gp)
            self.sort_prof('sig', 2, gp)
    ## \fn sort_profiles(self, gp)
    # sort all profiles
    # @param gp global parameters


    def set_analytic(self, x0, gp):
        from gl_project import rho_INT_Rho, rho_param_INT_Rho, rho_SUM_Mr
        if gp.investigate == 'gaia':
            r0 = x0 # [pc], spherical case
            self.analytic.set_x0(r0)
            rhoanalytic = ga.rhogaiatot_3D(r0)
            self.analytic.set_prof('rho', rhoanalytic, 0, gp)
            Manalytic = rho_SUM_Mr(r0, rhoanalytic)
            self.analytic.set_prof('M', Manalytic, 0, gp)            
            self.analytic.set_prof('nr', ga.nrgaiatot_3D(r0), 0, gp)
            for pop in range(1, gp.pops+1):
                # r0 scaling in pc, fine
                self.analytic.set_prof('beta', ga.betagaia(r0), pop, gp)
                be = ga.betagaia(r0)
                self.analytic.set_prof('betastar', be/(2.-be), pop, gp)
                analytic_3D = ga.nugaiatot_3D(r0)
                Sig_an = rho_INT_Rho(r0, analytic_3D)
                Signorm = gh.ipol_rhalf_log(r0, analytic_3D, gp.Rscale[pop])
                self.analytic.set_prof('Sig', Sig_an / Signorm, pop, gp)
                self.analytic.set_sig(0.*r0) # TODO: find analytical profile
        elif gp.investigate == 'walk':
            r0 = x0 # [pc], spherical case
            self.analytic.x0 = r0
            rhoanalytic = ga.rhowalk_3D(r0, gp)[0]
            self.analytic.set_prof('rho', rhoanalytic, 0, gp)
            Manalytic = rho_SUM_Mr(r0, rhoanalytic)
            self.analytic.set_prof('M', Manalytic, 0, gp)
            self.analytic.set_prof('nr', ga.nrwalktot_3D_deriv(r0, gp), 0, gp)
            # TODO: rho^* and Sigma^* profiles from rhowalk
            for pop in range(1, gp.pops+1):
                be = ga.betawalker(r0, gp)[pop-1]
                self.analytic.set_prof('beta', be, pop, gp)
                self.analytic.set_prof('betastar', be/(2.-be), pop, gp)
                analytic_3D = ga.rhowalk_3D(r0, gp)[pop]
                Sig_an = rho_INT_Rho(r0, analytic_3D)
                #Signorm = gh.ipol_rhalf_log(r0, Sig_an, gp.Rscale[pop])
                self.analytic.set_prof('Sig', Sig_an, pop, gp)
                self.analytic.set_prof('sig', 0.*r0-1, pop, gp) # TODO
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
        output.write(basename+'/prof_'+prof+'_'+str(pop)+'.ascii')

        if gp.investigate =='walk' or gp.investigate=='gaia':
            out_an = go.Output()
            out_an.add('radius [pc]', self.analytic.x0)
            out_an.add('analytic profile', self.analytic.get_prof(prof, pop))
            out_an.write(basename+'/prof_'+prof+'_'+str(pop)+'.analytic')
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
        output.write(basename+'/prof_chi2_0.ascii')
        return
    ## \fn write_chi2(self, basename, edges, bins)
    # write ascii file with chi2 information
    # @param basename directory string
    # @param edges array of edges
    # @param bins y-values


    def write_all(self, basename, gp):
        self.write_prof(basename, 'rho', 0, gp)
        self.write_prof(basename, 'M', 0, gp)
        if gp.geom=='sphere':
            self.write_prof(basename, 'nr', 0, gp)
        self.write_prof(basename, 'betastar', 1, gp)
        self.write_prof(basename, 'beta', 1, gp)
        self.write_prof(basename, 'Sig', 1, gp)
        self.write_prof(basename, 'sig', 1, gp)
        if gp.pops == 2:
            self.write_prof(basename, 'betastar', 2, gp)
            self.write_prof(basename, 'beta', 2, gp)
            self.write_prof(basename, 'Sig', 2, gp)
            self.write_prof(basename, 'sig', 2, gp)
    ## \fn write_all(self, basename, gp)
    # write output files for all profiles
    # @param basename directory string
    # @param gp global parameters


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
        output = go.Output()
        r0 = self.Mmedi.x0 # [pc]
        output.add('radius (data) [pc]', r0)
        if prof == 'Sig' or prof == 'Sigstar':
            # get 2D data here
            # Rbin [Rscale], Binmin [Rscale], Binmax [Rscale], Sig(R)/Sig(0) [1], error [1]
            dum,dum,dum,Sigdat,Sigerr = np.transpose(np.loadtxt(gp.files.Sigfiles[pop], \
                                                              unpack=False, skiprows=1))
            Sigdat *= gp.Sig0pc[pop] # [Msun/pc^2]
            Sigerr *= gp.Sig0pc[pop] # [Msun/pc^2]
            #Signorm = gh.ipol_rhalf_log(gp.xipol, Sigdat, gp.Rscale[pop])
            output.add('data [Msun/pc^2]', Sigdat)
            output.add('error [Msun/pc^2]', Sigerr)
            output.add('data - error [Msun/pc^2]', Sigdat-Sigerr)
            output.add('data + error [Msun/pc^2]', Sigdat+Sigerr)
            ax.fill_between(r0, Sigdat-Sigerr, Sigdat+Sigerr, \
                            color='blue', alpha=0.3, lw=1)
            ax.set_ylim([min(Sigdat-Sigerr)/2., 2.*max(Sigdat+Sigerr)])
        elif prof == 'sig':
            DATA = np.transpose(np.loadtxt(gp.files.sigfiles[pop], \
                                           unpack=False, skiprows=1))
            sigdat = DATA[4-1] # [maxsiglosi]
            sigerr = DATA[5-1] # [maxsiglosi]
            sigdat *= gp.maxsiglos[pop]  # [km/s]
            sigerr *= gp.maxsiglos[pop]  # [km/s]
            output.add('data [km/s]', sigdat)
            output.add('error [km/s]', sigerr)
            output.add('data - error [km/s]', sigdat-sigerr)
            output.add('data + error [km/s]', sigdat+sigerr)
            ax.fill_between(r0, sigdat-sigerr, sigdat+sigerr, \
                            color='blue', alpha=0.3, lw=1)
            ax.set_ylim([0., 2.*max(sigdat+sigerr)])
        output.write(basename+'/prof_'+prof+'_'+str(pop)+'.data')
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
        elif prof == 'rhostar':
            ax.set_ylabel('$\\rho^*\\quad[\\rm{M}_\\odot/\\rm{pc}^3]$')
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
        elif prof == 'Sig':
            ax.set_ylabel('$\\Sigma_'+str(pop)+'\\quad[\\rm{M}_\\odot/\\rm{pc}^2]$')
        elif prof == 'Sigstar':
            ax.set_ylabel('$\\Sigma^*\\quad[\\rm{M}_\\odot/\\rm{pc}^2]$')
        elif prof == 'sig':
            ax.set_ylabel('$\\sigma_{\\rm{LOS},'+str(pop)+'}\quad[\\rm{km}/\\rm{s}]$')
        return
    ## \fn plot_labels(self, ax, prof, pop, gp)
    # draw x and y labels
    # @param ax axis object to use
    # @param prof string of profile to look at
    # @param pop which population to analyze: 0 (all), 1, 2, ...
    # @param gp global parameters


    def plot_Rscale_3D(self, ax, gp):
        if gp.investigate == 'walk' or gp.investigate == 'gaia':
            if gp.pops == 1:
                rhodm, nu1 = ga.rhowalk_3D(gp.xipol, gp)
            else:
                rhodm, nu1, nu2 = ga.rhowalk_3D(gp.xipol, gp)
        for pop in range(gp.pops):
            # use our models
            nuprof = self.Mmedi.get_prof('nu', pop+1)
            if gp.investigate == 'walk' or gp.investigate == 'gaia':
                # or rather use analytic values, where available
                if pop == 0:
                    nuprof = nu1
                elif pop == 1:
                    nuprof = nu2

            Mprof = glp.rho_SUM_Mr(gp.xipol, nuprof)
            Mmax = max(Mprof) # Mprof[-1]
            ihalf = -1
            for kk in range(len(Mprof)):
                # half-light radius (3D) is where mass is more than half
                # ihalf gives the index of where this happens
                if Mprof[kk] >= Mmax/2 and ihalf <0:
                    xx = (gp.xipol[kk-1]+gp.xipol[kk])/2
                    ax.axvline(xx, color='green', lw=0.5, alpha=0.7)
                    ihalf = kk
    ## \fn plot_Rscale_3D(ax, gp)
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
                ax.axvline(gp.Rscale[pop+1], color='blue', lw=0.5) # [pc]
        else:
            self.plot_Rscale_3D(ax, gp)
        ax.set_xlim([r0[0],r0[-1]])
        # TODO: get better range
        #ax.set_xlim([gp.xipol[0], gp.xipol[-1]]) # gp.dat.binmax[-1]*gp.Rscale[0]])
        # if gp.pops==1:
        #     ax.set_xlim([r0[0], 5*gp.Rscale[1]])
        # else:
        #     ax.set_xlim([r0[0], 5*max(gp.Rscale[1], gp.Rscale[2])])
        return
    ## \fn fill_nice(self, ax, prof, pop, gp)
    # plot filled region for 1sigma and 2sigma confidence interval
    # @param ax axis object to plot into
    # @param prof string
    # @param pop int
    # @param gp


    def plot_profile(self, basename, prof, pop, gp):
        fig = plt.figure()
        ax  = fig.add_subplot(111)

        if prof != 'chi2':
            ax.set_xscale('log')

        if prof == 'rho' or prof == 'rhostar' or prof == 'Sig' or prof == 'M' or prof=='Sigstar':
            ax.set_yscale('log')

        self.plot_labels(ax, prof, pop, gp)

        if prof == 'chi2':
            goodchi = []
            for k in range(len(self.profs)):
                if self.subset[0] <= self.chis[k] <= self.subset[1]:
                    goodchi.append(self.chis[k])
            print('plotting profile '+prof+' for '+str(len(goodchi))+' models')
            bins, edges = np.histogram(np.log10(goodchi), range=[-2,6], \
                                       bins=max(6,np.sqrt(len(goodchi))),\
                                       density=True)
            ax.step(edges[1:], bins, where='pre')
            plt.draw()
            self.write_chi2(basename, edges, bins)
            fig.savefig(basename+'/prof_chi2_0.png')
            pp = PdfPages(basename+'/prof_chi2_0.pdf')
            pp.savefig(fig)
            pp.close()

            return

        self.fill_nice(ax, prof, pop, gp)
        self.plot_N_samples(ax, prof, pop)
        if prof == 'Sigstar' or prof == 'Sig' or prof == 'sig':
            self.plot_data(ax, basename, prof, pop, gp)
        else:
            if gp.investigate == 'walk' or gp.investigate == 'gaia':
                ax.plot(self.analytic.x0, self.analytic.get_prof(prof, pop), 'b--', lw=2)

        plt.draw()
        fig.savefig(basename+'/prof_'+prof+'_'+str(pop)+'.png')
        pp = PdfPages(basename+'/prof_'+prof+'_'+str(pop)+'.pdf')
        pp.savefig(fig)
        pp.close()
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
            ax.axvline(gp.Rscale[pop+1], color='gray', lw=0.5) # [pc]

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


    def __repr__(self):
        return "Profile Collection with "+str(len(self.profs))+" Profiles"
    ## \fn __repr__(self)
    # string representation for ipython

    def plot_overview(self, basename, gp):
        rho_data       = np.loadtxt(basename+'/prof_rho_0.ascii', skiprows=1, delimiter=',')
        nr_data        = np.loadtxt(basename+'/prof_nr_0.ascii', skiprows=1, delimiter=',')
        betastar1_data = np.loadtxt(basename+'/prof_betastar_1.ascii', skiprows=1, delimiter=',')
        beta1_data     = np.loadtxt(basename+'/prof_beta_1.ascii', skiprows=1, delimiter=',')
        Sig1_data      = np.loadtxt(basename+'/prof_Sig_1.ascii', skiprows=1, delimiter=',')
        sig1_data      = np.loadtxt(basename+'/prof_sig_1.ascii', skiprows=1, delimiter=',')
        if gp.pops == 2:
            betastar2_data = np.loadtxt(basename+'/prof_betastar_2.ascii', skiprows=1, delimiter=',')
            beta2_data     = np.loadtxt(basename+'/prof_beta_2.ascii', skiprows=1, delimiter=',')
            Sig2_data      = np.loadtxt(basename+'/prof_Sig_2.ascii', skiprows=1, delimiter=',')
            sig2_data      = np.loadtxt(basename+'/prof_sig_2.ascii', skiprows=1, delimiter=',')
        if gp.investigate == 'walk' or gp.investigate=='gaia':
            rho_analytic_data       = np.loadtxt(basename+'/prof_rho_0.analytic', skiprows=1, delimiter=',')
            nr_analytic_data        = np.loadtxt(basename+'/prof_nr_0.analytic', skiprows=1, delimiter=',')
            betastar1_analytic_data = np.loadtxt(basename+'/prof_betastar_1.analytic', skiprows=1, delimiter=',')
            beta1_analytic_data     = np.loadtxt(basename+'/prof_beta_1.analytic', skiprows=1, delimiter=',')
            Sig1_analytic_data      = np.loadtxt(basename+'/prof_Sig_1.analytic', skiprows=1, delimiter=',')
            sig1_analytic_data      = np.loadtxt(basename+'/prof_sig_1.analytic', skiprows=1, delimiter=',')
            if gp.pops == 2:
                betastar2_analytic_data = np.loadtxt(basename+'/prof_betastar_2.analytic', skiprows=1, delimiter=',')
                beta2_analytic_data     = np.loadtxt(basename+'/prof_beta_2.analytic', skiprows=1, delimiter=',')
                Sig2_analytic_data      = np.loadtxt(basename+'/prof_Sig_2.analytic', skiprows=1, delimiter=',')
                sig2_analytic_data      = np.loadtxt(basename+'/prof_sig_2.analytic', skiprows=1, delimiter=',')

        # produce figure with plots for "GOOD" and "ALL" model sets
        # +------+------+
        # | rho  | nr   |
        # +------+------+
        # | b*1  | beta1|
        # +------+------+
        # | Sig1 | sig1 |
        # +------+------+
        # | b*2  | beta2| if pops > 1
        # +------+------+
        # | Sig2 | sig2 |
        # +------+------+

        if gp.pops == 1:
            fig, ((ax_rho, ax_nr), (ax_bs1, ax_b1), (ax_Sig1, ax_sig1)) = plt.subplots(nrows=3, ncols=2, sharex=True)
        elif gp.pops == 2:
            fig, ((ax_rho, ax_nr), (ax_bs1, ax_b1), (ax_Sig1, ax_sig1), (ax_bs2, ax_b2), (ax_Sig2, ax_sig2)) = plt.subplots(nrows=5, ncols=2, sharex=True)

        # ax3 = subplot(313,  sharex=ax1, sharey=ax1)
        ax_rho.set_xscale('log')
        ax_nr.set_xscale('log')
        ax_bs1.set_xscale('log');
        ax_b1.set_xscale('log')
        ax_Sig1.set_xscale('log')
        ax_sig1.set_xscale('log')
        if gp.pops == 2:
            ax_bs2.set_xscale('log')
            ax_b2.set_xscale('log')
            ax_Sig2.set_xscale('log')
            ax_sig2.set_xscale('log')

        ax_rho.set_yscale('log')
        ax_nr.set_yscale('linear')
        ax_bs1.set_yscale('linear')
        ax_b1.set_yscale('linear')
        ax_Sig1.set_yscale('log')
        ax_sig1.set_yscale('linear')
        if gp.pops == 2:
            ax_bs2.set_yscale('linear')
            ax_b2.set_yscale('linear')
            ax_Sig2.set_yscale('log')
            ax_sig2.set_yscale('linear')
            ax_bs2.set_ylim([-1,1])
            ax_b2.set_ylim([-1,1])

        ax_nr.set_ylim([-1,6])
        ax_bs1.set_ylim([-1,1])
        ax_b1.set_ylim([-1,1])

        if gp.investigate=='walk' or gp.investigate=='gaia':
            self.plot_from_ascii(ax_rho,  rho_data,     rho_analytic_data, gp)
            self.plot_from_ascii(ax_nr,   nr_data,      nr_analytic_data,  gp)
            self.plot_from_ascii(ax_bs1,  betastar1_data, betastar1_analytic_data,gp)
            self.plot_from_ascii(ax_b1,   beta1_data,    beta1_analytic_data,gp)
            self.plot_from_ascii(ax_Sig1, Sig1_data,     Sig1_analytic_data, gp)
            self.plot_from_ascii(ax_sig1, sig1_data,     sig1_analytic_data, gp)
            if gp.pops == 2:
                self.plot_from_ascii(ax_bs2,  betastar2_data, betastar2_analytic_data,  gp)
                self.plot_from_ascii(ax_b2,   beta2_data,    beta2_analytic_data, gp)
                self.plot_from_ascii(ax_Sig2, Sig2_data,     Sig2_analytic_data, gp)
                self.plot_from_ascii(ax_sig2, sig2_data,     sig2_analytic_data, gp)

        fig.tight_layout()

        plt.draw()
        fig.savefig(basename+'/ov.png')
        pp = PdfPages(basename+'/ov.pdf')
        pp.savefig(fig)
        pp.close()
    ## \fn plot_overview(self, basename, gp)
    # plot a panel with all profiles
    # NOT USED ANYMORE
    # @param basename directory name for output files
    # @param gp global parameters


