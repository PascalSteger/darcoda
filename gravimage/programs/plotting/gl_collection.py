#!/usr/bin/env python3

## @file
# collect profiles and perform actions on them

# (c) 2014 ETHZ Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import numpy.random as npr
import math
import pdb

import matplotlib.pyplot as plt
plt.ioff()

from gl_class_profiles import Profiles
import plotting.gl_output as go
import gl_helper as gh
import gl_analytic as ga
import gl_project as glp
import gl_physics as phys

#USE_ALL = False
USE_ALL = True  # SS: Did not seem to affect things, what is this?

# The labels given in this function do not seem to correlate with what is actually plotted
def unit(prof):
    if prof == 'rho' or prof == 'nu':
        print ('prof:',prof)
        pdb.set_trace()
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
        print ('prof:',prof)
        unit = '[Msum/kpc^3]' 
    elif prof == 'sigmaRz_vec':
        unit = '[kpc]'
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
        self.BestFit = Profiles(pops, nipol)
        self.analytic = Profiles(pops, 100)
    ## \fn __init__(self, pops, nipol)
    # constructor
    # @param pops number of populations
    # @param nipol number of bins

    def add(self, prof):
        if np.isnan(prof.chi2):
            return
        self.chis.append(prof.chi2)
        #pdb.set_trace() 
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
        #self.subset = [minchi, maxchi]
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


    def sort_prof(self, prof, pop, gp): #SS: no longer in use(?)
        self.goodprof = []
        self.goodchi = []
        for k in range(len(self.profs)):
            if self.subset[0] <= self.chis[k] <= self.subset[1]:
                self.goodprof.append(self.profs[k].get_prof(prof, pop))
                self.goodchi.append(self.chis[k])
        
        if prof == 'nu_vec' or prof == 'sigz2_vec': # Printing best fit model:
            goodchi_arr = np.array(self.goodchi)
            minchi_index = np.argmin(goodchi_arr)
            if prof == 'nu_vec':
                print ('minchi:',goodchi_arr[minchi_index],' of',goodchi_arr.shape,'models')
            goodprof_arr = np.array(self.goodprof)
            #print ('goodprof shape:',goodprof_arr.shape)
            print ('best fit',prof,'profile:',self.goodprof[minchi_index])
        
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


    def weighted_sort_prof(self, prof, pop, gp):
        self.goodprof = []
        self.goodchi = []
        self.prof_weights=[]
        for k in range(len(self.profs)):
            if self.subset[0] <= self.chis[k] <= self.subset[1]:
                self.goodprof.append(self.profs[k].get_prof(prof, pop))
                self.goodchi.append(self.chis[k])
                self.prof_weights.append(self.profs[k].mn_weight)
        #pdb.set_trace()
        ## BODGE Renormalizing so that weights sum to 1  TODO  FIXME !!!
        #self.prof_weights = self.prof_weights*(10./np.sum(self.prof_weights))
        print ('Sum of weights:',np.sum(self.prof_weights))
        #pdb.set_trace()
        bin_prof_values = np.array(self.goodprof).transpose() # values for this bin from all the profiles
        bin_prof_weights = []
        bin_prof_credreg_bounds=[] # credibility region boundaries for this bin

        #print (self.goodprof[np.argmin(self.goodchi)])

        for lter in range(len(bin_prof_values)):
            sort_indices = np.argsort(bin_prof_values[lter])
            bin_prof_values[lter] = bin_prof_values[lter][sort_indices]
            bin_prof_weights.append(np.array(self.prof_weights)[sort_indices])

            cred_reg_bounds = [0.025, 0.16, 0.5, 0.84, 0.975]
            bin_i_prof_credreg_bounds = [bin_prof_values[lter][0]]

            for mter in range(len(cred_reg_bounds)):
                bin_weight_sum = 0.0

                for nter in range(len(bin_prof_values[lter])):
                    bin_weight_sum += bin_prof_weights[lter][nter]

                    if bin_weight_sum >= cred_reg_bounds[mter]:
                        bin_i_prof_credreg_bounds.append(bin_prof_values[lter][nter])
                        #print('found boundary, bin_weight_sum:',bin_weight_sum,'credreg:',cred_reg_bounds[mter])
                        break
               # pdb.set_trace()

            bin_i_prof_credreg_bounds.append(bin_prof_values[lter][-1])
            #print(bin_i_prof_credreg_bounds)
            bin_prof_credreg_bounds.append(bin_i_prof_credreg_bounds)


        credreg_bound_profs = np.array(bin_prof_credreg_bounds).transpose()

        if prof == 'rho_DM_vec': # Possibility to rescale plots
            rescale_prof = 1e6
        else:
            rescale_prof = 1.
            
        print ('Sorting prof:',prof)
        #print ('credreg_bound_profs:',credreg_bound_profs)
        #pdb.set_trace()
        self.Mmin.set_prof(prof,  credreg_bound_profs[0]/rescale_prof, pop, gp)
        self.M95lo.set_prof(prof, credreg_bound_profs[1]/rescale_prof, pop, gp)
        self.M68lo.set_prof(prof, credreg_bound_profs[2]/rescale_prof, pop, gp)
        self.Mmedi.set_prof(prof, credreg_bound_profs[3]/rescale_prof, pop, gp)
        self.M68hi.set_prof(prof, credreg_bound_profs[4]/rescale_prof, pop, gp)
        self.M95hi.set_prof(prof, credreg_bound_profs[5]/rescale_prof, pop, gp)
        self.Mmax.set_prof(prof,  credreg_bound_profs[6]/rescale_prof, pop, gp)
        self.BestFit.set_prof(prof, self.goodprof[np.argmin(self.goodchi)]/rescale_prof, pop, gp)
        print ('BestFit profile for',prof,':',self.goodprof[np.argmin(self.goodchi)])
        #pdb.set_trace()
        return credreg_bound_profs
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

    def weighted_sort_profiles_disc(self, gp):
        self.weighted_sort_prof('nu_vec', 0, gp)
        #self.sort_prof('sig_vec', 0, gp) #Old plots
        self.weighted_sort_prof('sigz2_vec', 0, gp)

        self.weighted_sort_prof('kz_rho_DM_vec', 0, gp)
        self.weighted_sort_prof('kz_nu_vec', 0, gp)

        self.weighted_sort_prof('rho_DM_vec', 0, gp)
        self.weighted_sort_prof('Sig_DM_vec', 0, gp)

        if gp.baryonmodel not in ['none', 'simplenu_baryon']:
            return

        self.weighted_sort_prof('rho_baryon_vec', 0, gp)
        self.weighted_sort_prof('Sig_baryon_vec', 0, gp)

        self.weighted_sort_prof('rho_total_vec', 0, gp)
        self.weighted_sort_prof('Sig_total_vec', 0, gp)

        if gp.tilt:
            self.weighted_sort_prof('sigmaRz_vec', 0, gp)
            self.weighted_sort_prof('tilt_vec', 0, gp)
            self.weighted_sort_prof('R_param',0,gp)


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
        #output.add('Best fit model '+uni, 
        print ('******************************************')
        profile_to_plot = self.Mmedi.get_prof(prof, pop)
        if prof == 'R_param':
            print ('Low 95% for',prof,':',self.M95lo.get_prof(prof, pop)[0])
            print ('Low 68% for',prof,':',self.M68lo.get_prof(prof, pop)[0])
            print ('Median for',prof,':',self.Mmedi.get_prof(prof, pop)[0])
            print ('High 68% for',prof,':',self.M68hi.get_prof(prof, pop)[0])
            print ('High 95% for',prof,':',self.M95hi.get_prof(prof, pop)[0])
        else:
            print ('Median for',prof,':',profile_to_plot)

        if prof == 'nu_vec':
            dum,dum,dum,nudat,nuerr = np.transpose(np.loadtxt(gp.files.nufiles[pop], \
                                                              unpack=False, skiprows=1))
            print ('nudat:',nudat)
            print ('nuerr:',nuerr)
            chi2_nu = np.sum(np.square((nudat-profile_to_plot)/nuerr))/np.size(nudat)
            print ('chi2 for nu part:',chi2_nu)
        elif prof == 'sigz2_vec':
            DATA = np.transpose(np.loadtxt(gp.files.sigfiles[pop], \
                                           unpack=False, skiprows=1))
            sigz2dat = DATA[4-1] # [maxsiglosi]
            sigz2err = DATA[5-1] # [maxsiglosi]
            print ('sigz2dat:',sigz2dat)
            print ('sigz2err:',sigz2err)
            chi2_sigz2 = np.sum(np.square((sigz2dat-profile_to_plot)/sigz2err))/np.size(sigz2dat)
            print ('chi2 for sigz2 part:',chi2_sigz2)
            
        #print (self.M95hi.get_prof(prof, pop))
        #pdb.set_trace()
        if prof != 'tilt_vec' and prof != 'R_param':
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
        
        if gp.tilt:
            self.write_prof(basename, 'sigmaRz_vec', 0, gp)
            self.write_prof(basename, 'tilt_vec', 0, gp)
            self.write_prof(basename, 'R_param', 0, gp)

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

        elif prof == 'sigmaRz_vec':
            DATA = np.transpose(np.loadtxt(gp.files.tiltfiles[pop], \
                                           unpack=False, skiprows=1))
            tiltdat = DATA[4-1] # [maxsiglosi]
            tilterr = DATA[5-1] # [maxsiglosi]
            #output.add('data [km^2/s^2]', np.sqrt(tiltdat))  # SS: are these used ???
            #output.add('error [km^2/s^2]', tilterr)
            #output.add('data - error [km/s]', tiltdat-tilterr)
            #output.add('data + error [km/s]', tiltdat+tilterr)
            #pdb.set_trace()
            ax.fill_between(r0, tiltdat-tilterr, tiltdat+tilterr, \
                            color='blue', alpha=0.3, lw=1)
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
            #ax.set_ylim(1.0E4, 4.0E5)
        elif prof == 'sig_vec':
            ax.set_ylabel('$\\sigma_{z}\\quad[\\rm{km}/\\rm{s}]$')
            ax.set_ylim(10., 40.)
        elif prof == 'sigz2_vec':
            ax.set_ylabel('$\\sigma_{z}^2\\quad[\\rm{km}^2/\\rm{s}^2]$')
            #ax.set_ylim(500., 1100.)

        elif prof == 'kz_nu_vec':
            ax.set_ylabel('$k(z)_{\\nu}$')
            #ax.set_ylim(0., 6.)
        elif prof == 'kz_rho_DM_vec':
            ax.set_ylabel('$k(z)_{\\rho, \\rm{DM}}$')
            ax.set_ylim(-10., 10.)

        elif prof == 'rho_DM_vec':
            ax.set_ylabel('$\\rho_{\\rm{DM}}\\quad[10^6\\rm{M}_\\odot/\\rm{kpc}^3]$')
            #ax.set_ylim(1E6, 1E8)
        elif prof == 'Sig_DM_vec':
            ax.set_ylabel('$\\Sigma_{\\rm{DM}} \\quad[\\rm{M}_\\odot/\\rm{kpc}^2]$')
            #ax.set_ylim(0,1.0E8)

        elif prof == 'rho_baryon_vec':
            ax.set_ylabel('$\\rho_{\\rm{baryon}}\\quad[\\rm{M}_\\odot/\\rm{kpc}^3]$')
            #ax.set_ylim(1E3, 1E9)
        elif prof == 'Sig_baryon_vec':
            ax.set_ylabel('$\\Sigma_{\\rm{baryon}} \\quad[\\rm{M}_\\odot/\\rm{kpc}^2]$')
            #ax.set_ylim(0,1.0E8)

        elif prof == 'rho_total_vec':
            ax.set_ylabel('$\\rho_{\\rm{total}}\\quad[\\rm{M}_\\odot/\\rm{kpc}^3]$')
            #ax.set_ylim(1E3, 1E9)
        elif prof == 'Sig_total_vec':
            ax.set_ylabel('$\\Sigma_{\\rm{total}} \\quad[\\rm{M}_\\odot/\\rm{kpc}^2]$')
            #ax.set_ylim(0,1.0E8)

        elif prof == 'sigmaRz_vec':
            ax.set_ylabel('$\\sigma_{Rz}\\quad[\\rm{km}^2/\\rm{s}^2]$')



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
        BestFit = self.BestFit.get_prof(prof, pop)
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
        #ax.plot(r0, BestFit, 'm', lw=.5)  # 'm'
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

        if prof == 'chi2':  # SS: Could not find where x-axis is set to logscale
            fig = plt.figure()
            ax  = fig.add_subplot(111)
            ax.set_yscale('log') 

        if prof == 'rho' or prof == 'Sig' or\
           prof == 'M' or prof == 'nu':
            fig = plt.figure()
            ax  = fig.add_subplot(111)
            ax.set_yscale('log')

        if prof == 'nu_vec' or  prof == 'rho_baryon_vec' or prof == 'rho_total_vec' or prof == 'sigz2_vec':
            fig = plt.figure()
            ax  = fig.add_subplot(111)
            ax.set_xscale('linear')  # SS: was 'log' before
            ax.set_yscale('log')

        if prof == 'rho_DM_vec':
            #fig = plt.figure()
            right_hand_side_label = False  # GeV/cm3 on right hand sider or not
            if (right_hand_side_label):
                fig = plt.figure(figsize=(2.3,2.3))
                ax  = fig.add_subplot(111)
                plt.subplots_adjust(left=0.21, right=0.79, bottom=0.22, top=0.8)
                # Above fits left and right labels but shrinks plot surface
                #plt.subplots_adjust(left=0.1, right=0.75)
                # Above fits only right label but preseves plot surface square
                # Python3 default: left=0.125, right=0.9, bottom=0.1, top=0.9
                ax2 = ax.twinx()
                ax2.set_ylabel('$\\rho_{\\rm{DM}}\\quad[\\rm{GeV}/\\rm{cm}^3]$')
                ax_ylim = np.array([6,16])    # TAG 
                kpc = 3.0857E19   # kpc in m
                Msun = 1.9891E30  # Sun's mass in kg
                GeV = 1.78266E-27 # GeV in kg
                ax_ylim2 = 1E6*Msun/(GeV*(100.*kpc)**3)*ax_ylim
                ax.set_ylim(ax_ylim) 
                ax2.set_ylim(ax_ylim2)
            else:
                fig = plt.figure()
                ax  = fig.add_subplot(111)
                #ax.set_ylim([0,14])    # TAG 
            ax.set_xscale('linear')  # SS: changed rhoDM plots to linear scale
            ax.set_yscale('linear')  # plot this in log or linear scale...?



        if prof == 'kz_rho_DM_vec' or prof == 'kz_nu_vec' or prof == 'sig_vec' or prof == 'Sig_DM_vec'  or prof == 'Sig_baryon_vec'  or prof == 'Sig_total_vec' or prof == 'sigmaRz_vec':
            fig = plt.figure()
            ax  = fig.add_subplot(111)
            ax.set_xscale('linear')  # SS: was 'log' before
            ax.set_yscale('linear')

        #if prof != 'chi2':
        #    ax.set_xscale('log')

        self.plot_labels(ax, prof, pop, gp)

        if len(self.profs)>0:
            if prof == 'chi2':
                goodchi = []
                for k in range(len(self.profs)):
                    # do include all chi^2 values for plot
                    goodchi.append(self.chis[k])
                print('plotting profile chi for '+str(len(goodchi))+' models')
                print('min, max, maxsubset found chi2: ', min(self.chis), max(self.chis), self.subset[1])
                #pdb.set_trace()
                #chi_arr = np.array(self.chis) 
                #print(chi_arr[np.where(chi_arr < 1.87)])
                bins, edges = np.histogram(np.log10(goodchi), range=[-3,3], \
                                           bins=max(6,np.sqrt(len(goodchi))),\
                                           density=True)
                # BODGE above to get the code running. Does it make any difference?
                # Exchanged range = [-3,3] to range = [-3,4]   TODO FIXME

                ax.step(edges[1:], bins, where='pre')
                ax.xaxis.get_majorticklabels()[0].set_horizontalalignment("left")
                plt.draw()
                self.write_chi2(basename, edges, bins)
                fig.savefig(basename+'/output/prof_chi2_0.pdf')
                return

            self.fill_nice(ax, prof, pop, gp)
            #self.plot_N_samples(ax, prof, pop) #SS: Don't plot thin gray lines
            self.plot_bins(ax, prof, gp)
            if prof == 'Sig' or prof=='nu' or prof == 'sig' or prof == 'nu_vec' or prof == 'sig_vec' or prof == 'sigz2_vec' or prof == 'sigmaRz_vec':
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

    def true_sigz2_func(self,binmin,binmax,nsmallbin,nbin,z0,Ntr,nu0,K,D,F,gp):
        #C = 22.85**2
        C = 40.0**2
        truesig2_arr = np.zeros(nbin)
        for kbin in range(nbin):
            zvec = np.zeros(nsmallbin+1) # one extra point since lists both binmin & binmax
            zvec[0] = binmin[kbin]
            zvec[-1] = binmax[kbin]
            for k in range(1,nsmallbin+1):
                exptemp = math.exp(-zvec[k-1]/z0)-(math.exp(-zvec[0]/z0)-math.exp(-zvec[-1]/z0))/nsmallbin
                zvec[k] = -z0*math.log(exptemp)
            nuvec = nu0*np.exp(-zvec/z0)
            if gp.baryonmodel not in ['simplenu_baryon']:  # assume DM only
                Kzvec_total = -2.*F*zvec
            else:
                Kzvec_total = -((K*zvec)/(np.sqrt(zvec**2 + D**2)) + 2.*F*zvec)
            Sigma_z_total = (1000.**2)*abs(Kzvec_total)/(2*np.pi*4.299) #Msun kpc^-1
            #smallsig2_vec = sigz2(zvec,nuvec,C)
            tiltvec = np.zeros(len(nuvec))   # Tilt not implemented here
            smallsig2_vec = phys.sigz2(zvec, Sigma_z_total, tiltvec, nuvec, C*nuvec[0])
            sig2 = np.sum(smallsig2_vec)/nsmallbin
            C = smallsig2_vec[-1]
            truesig2_arr[kbin] = sig2
        return truesig2_arr
    ## \fn true_sigz2_func
    # Calculates the true sigz2 value for each bin (with bin borders binmin & binmax), by
    # dividing the bin into nsmallbin subbins, with theoret. equally many tracer stars per bin.
    # The true sigz2 value is then calculated as the average over the subbin results.


    def plot_model_simplenu(self, ax, basename, prof, gp):
        Ntr = 930811  # N.o. tracers in simple2_1e6 with z < 2.4 kpc
        #Ntr = 9516. 
        G1 = 4.299e-6   # Newton's constant in (km)^2*kpc/(Msun*s^2)
        K = 1500.
        F = 267.65
        D = 0.18
        K_dd = 300.  # For 'normal' dark disk
        #K_dd = 900.   # For bdd, i.e. big dark disk
        D_dd = 2.5
        z0 = 0.9  # Thick disk
        #z0 = 0.4  # Thin disk
        normC = 40.0
        #normC = 22.85

        zvec=np.linspace(0, gp.data_z_cut, 100)

        Kzvec_total = -((K*zvec)/(np.sqrt(zvec**2 + D**2)) + 2.*F*zvec)
        Kzvec_baryon = -((K*zvec)/(np.sqrt(zvec**2 + D**2)))
        Kzvec_DD = -((K_dd*zvec)/(np.sqrt(zvec**2 + D_dd**2)))
        Kzvec_const_DM = -(2.*F*zvec)

        rho_z_DM_const = (1/(4*np.pi*G1)) * abs(2.*F) * np.ones(len(zvec))
        if prof == 'nu_vec':
            print ('gp.dd_data:',gp.dd_data,' K_dd:',K_dd)
        dd_data = False
        #if gp.dd_data:
        if dd_data:
            rho_z_DM =  rho_z_DM_const +  (1/(4*np.pi*G1)) * abs((K_dd*(D_dd**2)/((D_dd**2 + zvec**2)**(1.5))))
            Kzvec_DM = Kzvec_const_DM + Kzvec_DD
        else:    
            rho_z_DM = rho_z_DM_const
            Kzvec_DM = Kzvec_const_DM
        #rho_z_total = (1/(4*np.pi*G1)) * abs((K*(D**2)/((D**2 + zvec**2)**(1.5))) + 2.*F)
        rho_z_baryon = (1/(4*np.pi*G1)) * abs((K*(D**2)/((D**2 + zvec**2)**(1.5))))
        rho_z_total = rho_z_baryon + rho_z_DM
        #Sigma_z_total = (1000.**2)*abs(Kzvec_total)/(2*np.pi*4.299)
        Sigma_z_baryon = (1000.**2)*abs(Kzvec_baryon)/(2*np.pi*4.299)
        Sigma_z_DM = (1000.**2)*abs(Kzvec_DM)/(2*np.pi*4.299)
        Sigma_z_total = Sigma_z_baryon + Sigma_z_DM

        k_z_rho_total = (3*(D**2)*K*zvec) / ( ((D**2 + zvec**2)**2.5) * ((K*(D**2))/((D**2 + zvec**2)**1.5) + 2*F))
        k_z_rho_baryon = (3*(D**2)*K*zvec) / ( ((D**2 + zvec**2)**2.5) * ((K*(D**2))/((D**2 + zvec**2)**1.5))) # not currently used

        # Backwards compatibility: if using old data, then all mass is in DM
        if gp.baryonmodel not in ['simplenu_baryon']:
            gh.LOG(1, 'Simplenu Analytic: No baryon model, all mass is in DM.')
            Sigma_z_DM = Sigma_z_total
            rho_z_DM = rho_z_total

        nu0 = Ntr/(z0*(1-np.exp(-gp.data_z_cut/z0)))
        nuvec = nu0*np.exp(-zvec/z0)

        #bincentermed, binmin, binmax, nudat, nuerr = gh.readcol5(gp.files.nufiles[0])
        bincentermed, binmin, binmax, dum, dum = gh.readcol5(gp.files.nufiles[0])
        nbin = np.size(bincentermed)
        truen_arr = nu0*z0*(np.exp(-1.*binmin/z0)-np.exp(-1.*binmax/z0))  # true n.o. stars in bins
        truenu_arr = truen_arr/(binmax-binmin)

        if prof == 'nu_vec':
            ##ax.plot(zvec, nuvec, 'g-', alpha=0.5)
            #ax.plot(bincentermed, truenu_arr, 'g-', alpha=0.5)
            # Uncomment above line to plot theory curve for nu
            dfflnsdlcn = 1.   # needs to have an indented block...

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
            ax.plot(zvec, rho_z_DM/1e6, 'g-', alpha = 0.5) # Plot in units of 10^6 Msun/kpc3

        # if all mass is described by DM, then plot kz_rho_DM_vec
        # if the simplenu_baryon model is used, then kz_rho_DM_vec should equal zero
        elif prof == 'kz_rho_DM_vec' and gp.baryonmodel in ['none', 'sim']: #i.e. if all mass is described by DM
            ax.plot(zvec, k_z_rho_total, 'g-', alpha = 0.5)


        elif prof == 'sigz2_vec':
            ##sigz2_analytic = phys.sigz2(zvec, Sigma_z_total, nuvec, (normC**2)*nuvec[0])
            ##ax.plot(zvec, sigz2_analytic, 'g-', alpha = 0.5)
            true_sigz2_arr = self.true_sigz2_func(binmin,binmax,1000,nbin,z0,Ntr,nu0,K,D,F,gp)
            #ax.plot(bincentermed, true_sigz2_arr, 'g-', alpha = 0.5)
            # Uncomment above line to plot theory curve for sigz2
            
        elif prof == 'sigmaRz_vec':
            A = -0.0087
            n = 1.44
            sigmaRz = A*(1000.*zvec)**n 
            ax.plot(zvec, sigmaRz, 'g-', alpha = 0.5)

        return

    def plot_bins(self, ax, prof, gp):
        [ax.axvline(x, color='c', linewidth=0.1) for x in gp.z_bincenters]
        return





    def __repr__(self):
        return "Profile Collection with "+str(len(self.profs))+" Profiles"
    ## \fn __repr__(self)
    # string representation for ipython
