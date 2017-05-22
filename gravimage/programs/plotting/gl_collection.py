#!/usr/bin/env python3

## @file
# collect profiles and perform actions on them

# (c) 2014 ETHZ Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import numpy.random as npr
import math
import pdb
import matplotlib.ticker as ticker

import matplotlib.pyplot as plt
plt.ioff()

from gl_class_profiles import Profiles
import plotting.gl_output as go
import gl_helper as gh
import gl_analytic as ga
import gl_project as glp
import gl_physics as phys

import barrett.posterior as posterior
import barrett.data as data

USE_ALL = True
#USE_ALL = True  # SS: Did not seem to affect things, what is this?

import numpy as np
from matplotlib import rcParams, rc
 
# common setup for matplotlib

# Font code form Jon: 
#fs = 14
#fs = 10
fs = 9
params = {'backend': 'pdf',
          'savefig.dpi': 300, # save figures to 300 dpi
          'axes.labelsize': fs,
          'font.size': fs,
          'legend.fontsize': fs,
          'xtick.labelsize': fs,
          'ytick.major.pad': 6,
          'xtick.major.pad': 6,
          'ytick.labelsize': fs,
          'text.usetex': True,
          'font.family':'sans-serif',
          # free font similar to Helvetica
          'font.sans-serif':'FreeSans'}
 
# use of Sans Serif also in math mode
rc('text.latex', preamble=r'\usepackage{sfmath}')
 
rcParams.update(params)
# ... end ... (font code from Jon)


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
    elif prof == 'sigmaRz_vec':
        unit = '[kpc]'
    else:
        unit = 'a.u.'
    return unit
## \fn unit(prof)
# return string with units for profile
# @param prof string for profile


class ProfileCollection():  # (reffered to as 'pc in plot_prof_disc)
    def __init__(self, pops, nipol, no_Sigrho_bins): # nipol = gp.nbins
        self.chis = []
        self.goodchi = []  # will contain all chi^2 for models in self.subset
        self.goodprof = [] # will contain the corresponding profiles
        self.profs = []
        self.subset = [0., np.inf]
        self.x0 = np.array([])
        self.binmin = np.array([])
        self.binmax = np.array([])
        self.Mmin  = Profiles(pops, nipol, no_Sigrho_bins)
        self.M997lo = Profiles(pops, nipol, no_Sigrho_bins)
        self.M95lo = Profiles(pops, nipol, no_Sigrho_bins)
        self.M68lo = Profiles(pops, nipol, no_Sigrho_bins)
        self.Mmedi = Profiles(pops, nipol, no_Sigrho_bins)
        self.M68hi = Profiles(pops, nipol, no_Sigrho_bins)
        self.M95hi = Profiles(pops, nipol, no_Sigrho_bins)
        self.M997hi = Profiles(pops, nipol, no_Sigrho_bins)
        self.Mmax  = Profiles(pops, nipol, no_Sigrho_bins)
        self.BestFit = Profiles(pops, nipol, no_Sigrho_bins)
        #self.analytic = Profiles(pops, 100)
    ## \fn __init__(self, pops, nipol)
    # constructor
    # @param pops number of populations
    # @param nipol number of bins # now nipol is an array ...[pop][0 or 1]

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
    # (SS: exporting the profiles in phys_MNoutput_profiles.save to self.profs)



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


    def set_z_values(self, z_vecs, nu_z_vecs, z_vecs_comb_w0):
        # Most of these seems to not be used. (probably only Mmedi.z_vecs used):
        self.Mmin.z_vecs = z_vecs
        self.M997lo.z_vecs = z_vecs
        self.M95lo.z_vecs = z_vecs
        self.M68lo.z_vecs = z_vecs
        self.Mmedi.z_vecs = z_vecs
        self.M68hi.z_vecs = z_vecs
        self.M95hi.z_vecs = z_vecs
        self.M997hi.z_vecs = z_vecs
        self.Mmax.z_vecs = z_vecs

        self.Mmin.nu_z_vecs = nu_z_vecs
        self.M997lo.nu_z_vecs = nu_z_vecs
        self.M95lo.nu_z_vecs = nu_z_vecs
        self.M68lo.nu_z_vecs = nu_z_vecs
        self.Mmedi.nu_z_vecs = nu_z_vecs
        self.M68hi.nu_z_vecs = nu_z_vecs
        self.M95hi.nu_z_vecs = nu_z_vecs
        self.M997hi.nu_z_vecs = nu_z_vecs
        self.Mmax.nu_z_vecs = nu_z_vecs

        self.Mmin.z_vecs_comb_w0 = z_vecs_comb_w0
        self.M997lo.z_vecs_comb_w0 = z_vecs_comb_w0
        self.M95lo.z_vecs_comb_w0 = z_vecs_comb_w0
        self.M68lo.z_vecs_comb_w0 = z_vecs_comb_w0
        self.Mmedi.z_vecs_comb_w0 = z_vecs_comb_w0
        self.M68hi.z_vecs_comb_w0 = z_vecs_comb_w0
        self.M95lo.z_vecs_comb_w0 = z_vecs_comb_w0
        self.M997lo.z_vecs_comb_w0 = z_vecs_comb_w0
        self.Mmax.z_vecs_comb_w0 = z_vecs_comb_w0

        print ('z values for rho and Sigma plots:',z_vecs_comb_w0)
        print ('z values for nu plots:',nu_z_vecs)
        print ('z values for sigz and sigRz plots:',z_vecs)


    def sort_prof(self, prof, pop, gp): #SS: no longer in use(?)
        print ('In sort_prof')
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

        if prof == 'rho_DM_vec': # Possibility to rescale plots
            #rescale_prof = 1e6
            rescale_prof = 1.
        else:
            rescale_prof = 1.


        self.Mmin.set_prof(prof,  tmp[0]/rescale_prof,       pop, gp)
        self.M95lo.set_prof(prof, tmp[ll*0.05]/rescale_prof, pop, gp)
        self.M68lo.set_prof(prof, tmp[ll*0.32]/rescale_prof, pop, gp)
        self.Mmedi.set_prof(prof, tmp[ll/2]/rescale_prof,    pop, gp)
        self.M68hi.set_prof(prof, tmp[ll*0.68]/rescale_prof, pop, gp)
        self.M95hi.set_prof(prof, tmp[ll*0.95]/rescale_prof, pop, gp)
        self.Mmax.set_prof(prof,  tmp[-1]/rescale_prof,      pop, gp)
        return tmp
    ## \fn sort_prof(self, prof, pop, gp)
    # sort the list of prof-profiles, and store the {1,2}sigma, min, medi, max in the appropriate place
    # @param prof profile identifier, 'rho', 'beta', ...
    # @param pop population identifier
    # @param gp global parameters
    # @return tmp profiles sorted binwise


    def weighted_sort_prof(self, prof, pop, gp):
        print ('In weighted_sort_prof: prof=',prof,' pop:',pop)
        # All pop indep quantities (rho_DM_vec etc) are only called for pop=0
        # (while the pop dep (like nu_vecs) are called for both pop=0 and 1).
        self.goodprof = []
        self.goodchi = []
        self.prof_weights=[]

        for k in range(len(self.profs)):
            if self.subset[0] <= self.chis[k] <= self.subset[1]:
                # Currently no chi cut in use  (i.e. ok if >= 0 & <= inf)
                self.goodprof.append(self.profs[k].get_prof(prof, pop))
                # Guess it would be possible to call eg rho_DM twice here
                # (for pop 1 & 2) and add it to one larger plot
                # (need to be coordinated with z_vec).
                self.goodchi.append(self.chis[k])
                self.prof_weights.append(self.profs[k].mn_weight)

        #Check weight distribution
        max_weight = max(self.prof_weights)
        sum_weight = sum(self.prof_weights)

        if 0.5 < max_weight < 0.9:
            print('WARNING: more than 0.5 of the weight is in one point')
            print('max weight = ', max_weight)
        elif 0.9 < max_weight:
            print('WARNING: more than 0.9 of the weight is in one point')
            print('max weight = ', max_weight)

        if sum_weight < 0.99:
            print('WARNING: weight sum is less than 0.99: sum = ', sum_weight)

        bin_prof_values = np.array(self.goodprof).transpose() # values for this bin from all the profiles
        #print ('shape(bin_prof_values):',np.shape(bin_prof_values))
        #print ('bpv[0]:',bin_prof_values[0])
        #print ('bpv[1]:',bin_prof_values[1])
        bin_prof_weights = []
        bin_prof_credreg_bounds=[] # credibility region boundaries for this bin

        #print (self.goodprof[np.argmin(self.goodchi)])

        for lter in range(len(bin_prof_values)):
            sort_indices = np.argsort(bin_prof_values[lter])
            bin_prof_values[lter] = bin_prof_values[lter][sort_indices]
            bin_prof_weights.append(np.array(self.prof_weights)[sort_indices])

            #cred_reg_bounds = [0.025, 0.16, 0.5, 0.84, 0.975] #0.0015, , 0.9985
            cred_reg_bounds = [0.0015, 0.025, 0.16, 0.5, 0.84, 0.975, 0.9985]
            bin_i_prof_credreg_bounds = [bin_prof_values[lter][0]]

            #self.prof_weights = np.ones(len(self.prof_weights))/(len(self.prof_weights))

            for mter in range(len(cred_reg_bounds)):
                bin_weight_sum = 0.0

                for nter in range(len(bin_prof_values[lter])):

                    bin_weight_sum += bin_prof_weights[lter][nter]

                    if bin_weight_sum >= cred_reg_bounds[mter]:
                        bin_i_prof_credreg_bounds.append(bin_prof_values[lter][nter])
                        #print('found boundary ', mter)
                        break

            bin_i_prof_credreg_bounds.append(bin_prof_values[lter][-1])
            #print(bin_i_prof_credreg_bounds)
            bin_prof_credreg_bounds.append(bin_i_prof_credreg_bounds)


        credreg_bound_profs = np.array(bin_prof_credreg_bounds).transpose()

        if prof == 'rho_DM_vec': # Possibility to rescale plots
            #rescale_prof = 1e6
            rescale_prof = 1.
        else:
            rescale_prof = 1.

        #print ('prof:',prof, 'pop: ', pop)


        #95-68 bounds
        #self.Mmin.set_prof(prof,  credreg_bound_profs[0]/rescale_prof, pop, gp)
        #self.M95lo.set_prof(prof, credreg_bound_profs[1]/rescale_prof, pop, gp)
        #self.M68lo.set_prof(prof, credreg_bound_profs[2]/rescale_prof, pop, gp)
        #self.Mmedi.set_prof(prof, credreg_bound_profs[3]/rescale_prof, pop, gp)
        #self.M68hi.set_prof(prof, credreg_bound_profs[4]/rescale_prof, pop, gp)
        #self.M95hi.set_prof(prof, credreg_bound_profs[5]/rescale_prof, pop, gp)
        #self.Mmax.set_prof(prof,  credreg_bound_profs[6]/rescale_prof, pop, gp)

        #99.7 95 68 bounds
        #print ('Before cbp: prof:',prof) 
        #print ('credreg_bound_profs:',credreg_bound_profs)
        self.Mmin.set_prof(prof,  credreg_bound_profs[0]/rescale_prof, pop, gp)
        self.M997lo.set_prof(prof, credreg_bound_profs[1]/rescale_prof, pop, gp)
        self.M95lo.set_prof(prof, credreg_bound_profs[2]/rescale_prof, pop, gp)
        self.M68lo.set_prof(prof, credreg_bound_profs[3]/rescale_prof, pop, gp)
        self.Mmedi.set_prof(prof, credreg_bound_profs[4]/rescale_prof, pop, gp)
        self.M68hi.set_prof(prof, credreg_bound_profs[5]/rescale_prof, pop, gp)
        self.M95hi.set_prof(prof, credreg_bound_profs[6]/rescale_prof, pop, gp)
        self.M997hi.set_prof(prof, credreg_bound_profs[7]/rescale_prof, pop, gp)
        self.Mmax.set_prof(prof,  credreg_bound_profs[8]/rescale_prof, pop, gp)





        self.BestFit.set_prof(prof, self.goodprof[np.argmin(self.goodchi)]/rescale_prof, pop, gp)
        #print ('BestFit profile for',prof,':',self.goodprof[np.argmin(self.goodchi)])

        return credreg_bound_profs
    ## \fn sort_prof(self, prof, pop, gp)
    # sort the list of prof-profiles, and store the {1,2}sigma, min, medi, max in the appropriate place
    # @param prof profile identifier, 'rho', 'beta', ...
    # @param pop population identifier
    # @param gp global parameters
    # @return tmp profiles sorted binwise


    def weighted_hist_heatmap(self,ax, prof, pop, gp):
        print('In weighted_hist_heatmap, prof = ', prof, 'pop = ', pop)
        self.goodprof = []
        self.goodchi = []
        self.prof_weights=[]
        for k in range(len(self.profs)):
            if self.subset[0] <= self.chis[k] <= self.subset[1]:
                self.goodprof.append(self.profs[k].get_prof(prof, pop))
                self.goodchi.append(self.chis[k])
                self.prof_weights.append(self.profs[k].mn_weight)

        bin_prof_values = np.array(self.goodprof).transpose() # values for this bin from all the profiles

        N_models = len(bin_prof_values[0])

        #Compile bin centers for all points

        bin_prof_values_bincents = [np.ones(N_models)*gp.z_bincenter_vecs[pop][ii] for ii in range(0, len(gp.z_bincenter_vecs[pop]))]

        #Construct list of all points
        y_data = np.concatenate(bin_prof_values)
        z_data = np.concatenate(bin_prof_values_bincents)
        weight_data = np.tile(self.prof_weights, len(gp.z_bincenter_vecs[pop]))

        if weight_data.all() ==0:
            weight_data = np.ones(len(weight_data)) * (1./N_models)

        #Set bin edges
        z_edges = np.append(gp.z_binmin_vecs[pop], gp.z_binmax_vecs[pop][-1])
        ax_yscale = ax.get_yscale()
        if ax_yscale == 'linear':
            y_edges = np.linspace(min(y_data), max(y_data), 100)
        elif ax_yscale =='log':
            y_edges = np.logspace(np.log10(min(y_data)), np.log10(max(y_data)), 100)

        #Construct histogram
        Hist, z_edges2, y_edges2 = np.histogram2d(z_data, y_data, weights = weight_data, bins=(z_edges, y_edges))

        #Mask zero values
        Hist_masked = np.ma.masked_where(Hist==0, Hist)

        #Plot color histogram
        if prof == 'rho_DM_vec': # Possibility to rescale plots
            #rescale_prof = 1e6
            rescale_prof = 1.
        else:
            rescale_prof = 1.

        ax.pcolormesh(z_edges, y_edges/rescale_prof, np.transpose(Hist_masked), cmap='BuPu')

        #Pull in Credible Region intervals

        print('weighted_hist_heatmap, prof = ', prof, 'pop = ', pop)
        M95lo = self.M95lo.get_prof(prof, pop)
        M68lo = self.M68lo.get_prof(prof, pop)
        Mmedi = self.Mmedi.get_prof(prof, pop)
        M68hi = self.M68hi.get_prof(prof, pop)
        M95hi = self.M95hi.get_prof(prof, pop)

        # TAGTAG !!

        [ax.hlines(M95lo[ii], z_edges[ii], z_edges[ii+1], lw=0.1, linestyle=':') for ii in range(0, gp.nbins[pop])]
        [ax.hlines(M68lo[ii], z_edges[ii], z_edges[ii+1], lw=0.1) for ii in range(0, gp.nbins[pop])]
        ax.plot(gp.z_bincenter_vecs[pop], Mmedi, 'r.', markersize=0.8)
        [ax.hlines(M68hi[ii], z_edges[ii], z_edges[ii+1], lw=0.1) for ii in range(0, gp.nbins[pop])]
        [ax.hlines(M95hi[ii], z_edges[ii], z_edges[ii+1], lw=0.1, linestyle=':') for ii in range(0, gp.nbins[pop])]

        return




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
        print ('In sort_profiles')
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
        print ('In sort_profiles_disc')
        for t_pop in range(0, gp.ntracer_pops):
            #self.sort_prof('kz_nu_vecs', t_pop, gp) # kz_nu_vecs not viable opt
            self.sort_prof('nu_vecs', t_pop, gp)
            self.sort_prof('sigz2_vecs', t_pop, gp)
            self.sort_prof('chi2_nu_vecs', t_pop, gp)

            if gp.tilt:
                self.sort_prof('sigmaRz2_vecs', t_pop, gp)
                self.sort_prof('tilt_vecs', t_pop, gp)

            #Legacy Support
            #self.sort_prof('kz_nu_vec', t_pop, gp)
            #self.sort_prof('nu_vec', t_pop, gp)
            #self.sort_prof('sigz2_vec', t_pop, gp)
            #if gp.tilt:
            #    self.sort_prof('sigmaRz_vec', t_pop, gp)
            #    self.sort_prof('tilt_vec', t_pop, gp)

        #self.sort_prof('kz_rho_DM_vec', 0, gp)  # stop entertaining kz_rho_DM

        self.sort_prof('rho_DM_vec', 0, gp)
        self.sort_prof('Sig_DM_vec', 0, gp)

        if gp.baryonmodel not in ['none', 'simplenu_baryon', 'trivial_baryon', 'obs_baryon', 'simplenu_baryon_gausian']:
            return

        self.sort_prof('rho_baryon_vec', 0, gp)
        self.sort_prof('Sig_baryon_vec', 0, gp)

        self.sort_prof('rho_total_vec', 0, gp)
        self.sort_prof('Sig_total_vec', 0, gp)


    def weighted_sort_profiles_disc(self, gp):
        print ('In weighted_sort_profiles_disc')
        for t_pop in range(0, gp.ntracer_pops):
            #self.weighted_sort_prof('kz_nu_vecs', t_pop, gp) # SS trying
            print ('Before w_s_p nu_vecs, t_pop:',t_pop)
            self.weighted_sort_prof('nu_vecs', t_pop, gp)
            print ('Before w_s_p sigz2_vecs, t_pop:',t_pop)
            self.weighted_sort_prof('sigz2_vecs', t_pop, gp)
            self.weighted_sort_prof('chi2_sigz2_vecs', t_pop, gp)
            self.weighted_sort_prof('chi2_sigRz2_vecs', t_pop, gp)

            if gp.tilt:
                self.weighted_sort_prof('sigmaRz2_vecs', t_pop, gp)
                self.weighted_sort_prof('tilt_vecs', t_pop, gp)

            #Legacy Support
            #self.weighted_sort_prof('kz_nu_vec', t_pop, gp)
            #self.weighted_sort_prof('nu_vec', t_pop, gp)
            #self.weighted_sort_prof('sigz2_vec', t_pop, gp)

            #if gp.tilt:
            #    self.weighted_sort_prof('sigmaRz_vec', t_pop, gp)
            #    self.weighted_sort_prof('tilt_vec', t_pop, gp)

        #self.weighted_sort_prof('kz_rho_DM_vec', 0, gp) # SS trying

        self.weighted_sort_prof('rho_DM_vec', 0, gp) 
        self.weighted_sort_prof('Sig_DM_vec', 0, gp) 

        if gp.baryonmodel not in ['none', 'simplenu_baryon', 'trivial_baryon', 'obs_baryon', 'simplenu_baryon_gaussian']:
            return

        self.weighted_sort_prof('rho_baryon_vec', 0, gp)
        self.weighted_sort_prof('Sig_baryon_vec', 0, gp)

        self.weighted_sort_prof('rho_total_vec', 0, gp)
        self.weighted_sort_prof('Sig_total_vec', 0, gp)




    def set_analytic(self, x0, gp):
        print ('In set_analytic')
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


    def write_prof(self, basename, prof, pop, gp):  # Used
        print ('******************************************')
        print ('In write_prof, prof=',prof)
        output = go.Output()
        uni = unit(prof)

        if prof in ['kz_rho_DM_vec', 'rho_DM_vec', 'Sig_DM_vec', 'rho_baryon_vec', 'Sig_baryon_vec',
            'rho_total_vec', 'Sig_total_vec']:
            output.add('radius (models) [pc]', self.Mmedi.z_vecs_comb_w0)
            #print ('z_vec:',self.Mmedi.z_vecs_comb_w0)
        elif prof in ['nu_vecs']:
            output.add('radius (models) [pc]', self.Mmedi.nu_z_vecs[pop])
        elif prof in ['sigz2_vecs', 'tilt_vecs', 'sigmaRz2_vecs']: # , 'kz_nu_vecs', 'z_vecs']
            output.add('radius (models) [pc]', self.Mmedi.z_vecs[pop])

        output.add('M 99.7% CL low ' + uni, self.M997lo.get_prof(prof, pop))
        output.add('M 95% CL low ' + uni, self.M95lo.get_prof(prof, pop))
        output.add('M 68% CL low ' + uni, self.M68lo.get_prof(prof, pop))
        output.add('M median '     + uni, self.Mmedi.get_prof(prof, pop))
        output.add('M 68% CL high '+ uni, self.M68hi.get_prof(prof, pop))
        output.add('M 95% CL high '+ uni, self.M95hi.get_prof(prof, pop))
        output.add('M 99.7% CL high '+ uni, self.M997hi.get_prof(prof, pop))
        #output.add('Best fit model '+uni,
        profile_to_plot = self.Mmedi.get_prof(prof, pop)
        print ('pop:',pop)
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

        if prof != 'tilt_vec':
            #print ('Before: prof=',prof,' pop=',pop,' *************')
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
        print ('In write_all')
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
        #self.write_prof(basename, 'kz_rho_DM_vec', 0, gp) #DELETE ME
        for t_pop in range(0, gp.ntracer_pops):
            #self.write_prof(basename, 'kz_nu_vecs', t_pop, gp)
            self.write_prof(basename, 'nu_vecs', t_pop, gp)
            self.write_prof(basename, 'sigz2_vecs', t_pop, gp)
            if gp.tilt:
                self.write_prof(basename, 'sigmaRz2_vecs', t_pop, gp)
                self.write_prof(basename, 'tilt_vecs', t_pop, gp)
            #self.write_prof(basename, 'kz_nu_vec', t_pop, gp)
            #self.write_prof(basename, 'nu_vec', t_pop, gp)
            #self.write_prof(basename, 'sigz2_vec', t_pop, gp)
            #if gp.tilt:
            #    self.write_prof(basename, 'sigmaRz_vec', t_pop, gp)
            #    self.write_prof(basename, 'tilt_vec', t_pop, gp)

            self.write_prof(basename, 'rho_DM_vec', t_pop, gp)
            self.write_prof(basename, 'Sig_DM_vec', t_pop, gp)

            self.write_prof(basename, 'rho_baryon_vec', t_pop, gp)
            self.write_prof(basename, 'Sig_baryon_vec', t_pop, gp)

            self.write_prof(basename, 'rho_total_vec', t_pop, gp)
            self.write_prof(basename, 'Sig_total_vec', t_pop, gp)
         # SS trying, above was like kz_... before
        #self.write_prof(basename, 'kz_rho_DM_vec', 0, gp)


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
        if prof == 'nu_vecs':
            r0 = self.Mmedi.nu_z_vecs[pop] # [pc]
        else:
            r0 = self.Mmedi.z_vecs[pop] # [pc]
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
        elif prof == 'nu' or prof == 'nu_vec' or prof == 'nu_vecs':
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
            #ax.fill_between(r0, nudat-nuerr, nudat+nuerr, \
            #                color='blue', alpha=0.3, lw=1)

            #ax.errorbar(r0,nudat,yerr=nuerr, linestyle="None",capsize = 1.,color='blue')    # TAG nu
            ax.errorbar(r0,nudat,yerr=nuerr, fmt='.',capsize=1.5,color='blue',zorder=10)
            #ax.errorbar(r0,nudat,yerr=nuerr, capsize = 0.05 ,color='blue') 
            print ('r0:',r0)
            print ('nudat:',nudat)
            #pdb.set_trace()
            #ax.set_ylim([min(nudat-nuerr)/2., 2.*max(nudat+nuerr)]) #bodge
        elif prof == 'sig' or prof == 'sig_vec' or prof == 'sig_vec':
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
            #ax.fill_between(r0, sigdat-sigerr, sigdat+sigerr, \
            #                color='blue', alpha=0.3, lw=1)
            ax.errorbar(r0,sigdat,yerr=sigerr, fmt='.',capsize = 1.5,color='blue',zorder=10) 
            #ax.errorbar(r0,sigdat,yerr=sigerr, linestyle="None",capsize = 1.,color='blue') 
            #ax.set_ylim([0., 2.*max(sigdat+sigerr)]) #bodge
        elif prof == 'sigz2_vec' or prof == 'sigz2_vecs':
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
            #ax.fill_between(r0, sigz2dat-sigz2err, sigz2dat+sigz2err, \
            #                color='blue', alpha=0.3, lw=1)

            ax.errorbar(r0,sigz2dat,yerr=sigz2err, fmt='.',capsize = 1.5,color='blue',zorder=10)
            #ax.errorbar(r0,sigz2dat,yerr=sigz2err, linestyle="None",capsize = 1.,color='blue') 
            #ax.set_ylim([0., 2.*max(sigdat+sigerr)]) #bodge

        elif prof == 'sigmaRz_vec' or prof == 'sigmaRz2_vecs':
            DATA = np.transpose(np.loadtxt(gp.files.sigRz2_files[pop], \
                                           unpack=False, skiprows=1))
            tiltdat = DATA[4-1] # [maxsiglosi]
            tilterr = DATA[5-1] # [maxsiglosi]
            #output.add('data [km^2/s^2]', np.sqrt(tiltdat))  # SS: are these used ???
            #output.add('error [km^2/s^2]', tilterr)
            #output.add('data - error [km/s]', tiltdat-tilterr)
            #output.add('data + error [km/s]', tiltdat+tilterr)

            print ('r0:',r0)
            print ('tiltdat-tilterr=',tiltdat-tilterr)
            print ('tiltdat+tilterr=',tiltdat+tilterr)
            # Plot the blue band showing the data:
            #ax.fill_between(r0, tiltdat-tilterr, tiltdat+tilterr, \
            #                    color='blue', alpha=0.3, lw=1)
            ax.errorbar(r0,tiltdat,yerr=tilterr, fmt='.',capsize = 1.5,color='blue',zorder=10) 
            #ax.errorbar(r0,tiltdat,yerr=tilterr, linestyle="None",capsize = 1.,color='blue') 
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
        elif prof == 'nu_vec' or prof == 'nu_vecs':
            ax.set_ylabel('$\\nu_{\\rm{Tr},'+str(pop)+'}\\quad[\\rm{stars}/\\rm{kpc}^3]$')
            #ax.set_ylim(1.0E4, 4.0E5)
        elif prof == 'sig_vec' or prof == 'sig_vecs':
            ax.set_ylabel('$\\sigma_{z}\\quad[\\rm{km}/\\rm{s}]$')
            ax.set_ylim(10., 40.)
        elif prof == 'sigz2_vec' or prof == 'sigz2_vecs':
            if gp.fit_to_sigz2:
                ax.set_ylabel('$\\sigma_{z,'+str(pop)+'}^2\\quad[\\rm{km}^2/\\rm{s}^2]$')
            else:
                ax.set_ylabel('$\\sigma_{z,'+str(pop)+'}\\quad[\\rm{km}/\\rm{s}]$')
                #ax.set_ylim(36., 48.)  # Useful for plotting old pop. TAG

        elif prof == 'kz_nu_vec' or prof == 'kz_nu_vecs':
            ax.set_ylabel('$k(z)_{\\nu,'+str(pop)+'} $')
            #ax.set_ylim(0., 6.)
        elif prof == 'kz_rho_DM_vec':
            ax.set_ylabel('$k(z)_{\\rho, \\rm{DM}}$')
            ax.set_ylim(-10., 10.)

        elif prof == 'rho_DM_vec':
            #ax.set_ylabel('$\\rho_{\\rm{DM}}\\quad[10^{-3}\\rm{M}_\\odot/\\rm{pc}^3]$')
            ax.set_ylabel('$\\rho_{\\rm{DM}}\\quad[\\rm{M}_\\odot/\\rm{pc}^3]$')
            #ax.set_ylabel('$\\rho_{\\rm{DM}}\\quad[\\rm{M}_\\odot/\\rm{kpc}^3]$')
            # Do not set ylim here!!  Set rho_ylim instead
        elif prof == 'Sig_DM_vec':
            ax.set_ylabel('$\\Sigma_{\\rm{DM}} \\quad[\\rm{M}_\\odot/\\rm{pc}^2]$')
            #ax.set_ylim(0,1.0E8)

        elif prof == 'rho_baryon_vec':
            ax.set_ylabel('$\\rho_{\\rm{baryon}}\\quad[\\rm{M}_\\odot/\\rm{pc}^3]$')
            #ax.set_ylim(1E3, 1E9)
        elif prof == 'Sig_baryon_vec':
            ax.set_ylabel('$\\Sigma_{\\rm{baryon}} \\quad[\\rm{M}_\\odot/\\rm{pc}^2]$')
            #ax.set_ylim(0,1.0E8)

        elif prof == 'rho_total_vec':
            #ax.set_ylabel('$\\rho_{\\rm{total}}\\quad[\\rm{M}_\\odot/\\rm{kpc}^3]$')
            ax.set_ylabel('$\\rho\\quad[\\rm{M}_\\odot/\\rm{pc}^3]$')
            #ax.set_ylim(1E3, 1E9)
            ax.set_ylim(0., 0.3E8)

        elif prof == 'Sig_total_vec':
            #ax.set_ylabel('$\\Sigma_{\\rm{total}} \\quad[\\rm{M}_\\odot/\\rm{kpc}^2]$')
            ax.set_ylabel('$\\Sigma \\quad[\\rm{M}_\\odot/\\rm{pc}^2]$')
            #ax.set_ylim(0,1.0E8)

        elif prof == 'sigmaRz_vec' or prof == 'sigmaRz2_vecs':
            ax.set_ylabel('$\\sigma_{Rz,'+str(pop)+'}\\quad[\\rm{km}^2/\\rm{s}^2]$')
            #ax.set_ylim(-10,100)




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
        print('fill_nice, prof = ', prof, 'pop = ', pop)
        #pdb.set_trace()
        if prof in ['kz_rho_DM_vec', 'rho_DM_vec', 'Sig_DM_vec', 'rho_baryon_vec',
            'Sig_baryon_vec', 'rho_total_vec', 'Sig_total_vec']:
            r0 = self.Mmedi.z_vecs_comb_w0
        elif prof in ['nu_vecs']:
            r0 = self.Mmedi.nu_z_vecs[pop]
        elif prof in ['z_vecs', 'kz_nu_vecs', 'nu_vecs', 'sigz2_vecs', 'tilt_vecs','sigmaRz2_vecs', 'chi2_sigz2_vecs', 'chi2_sigRz2_vecs']:
            r0 = self.Mmedi.z_vecs[pop]

        M95lo = self.M95lo.get_prof(prof, pop)
        M68lo = self.M68lo.get_prof(prof, pop)
        Mmedi = self.Mmedi.get_prof(prof, pop)
        BestFit = self.BestFit.get_prof(prof, pop)

        M68hi = self.M68hi.get_prof(prof, pop)
        M95hi = self.M95hi.get_prof(prof, pop)
        ax.fill_between(r0, M95lo, M95hi, color='black', alpha=0.2, lw=0.1)
        ax.plot(r0, M95lo, color='black', lw=0.4)
        ax.plot(r0, M95hi, color='black', lw=0.3)
        ax.fill_between(r0, M68lo, M68hi, color='black', alpha=0.4, lw=0.1)
        ax.plot(r0, M68lo, color='black', lw=0.4)
        ax.plot(r0, M68hi, color='black', lw=0.3)
        ax.plot(r0, Mmedi, 'r', lw=1) # Drawing the red median line
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

        if prof == 'Sig_baryon_vec':
            young_bary_Sig_min = np.array([29038624.,29556174.,30029287.,30467336.,30874921.,
                                           31255877.,31613414.,31950222.,32268565.,32570353.,
                                           32857204.,33130496.,33391406.,33640946.,33879988.,
                                           34109291.,34329518.,34541249.,34745001.,34941229.])
            young_bary_Sig_max = np.array([49387938.,49677359.,49929360.,50149362.,50341979.,
                                           50518227.,50674976.,50813817.,50956433.,51092113.,
                                           51213229.,51321730.,51419281.,51507309.,51587036.,
                                           51659509.,51725628.,51786165.,51841785.,51893060.])

            old_bary_Sig_min = np.array([30638698.,31610401.,32438533.,33149405.,33779845.,
                                         34315761.,34790192.,35223447.,35620826.,35986456.,
                                         36323654.,36635168.,36923329.,37190157.,37434538.,
                                         37658100.,37865959.,38059362.,38239427.,38407164.])
            old_bary_Sig_max = np.array([49989901.,50476024.,50841798.,51122426.,51342212.,
                                         51518001.,51668249.,51799780.,51909713.,52003296.,
                                         52084249.,52155234.,52218182.,52274515.,52325297.,
                                         52371342.,52413282.,52451622.,52486769.,52519062.])

            both_bary_Sig_min = np.array([29295203.,30530919.,31540943.,32389533.,33126345.,
                                          33776613.,34322334.,34803158.,35240437.,35640319.,
                                          36007474.,36345574.,36657595.,36946019.,37212967.,
                                          37460279.,37689582.,37902325.,38099811.,38283223.])
            both_bary_Sig_max = np.array([49607396.,50316804.,50820972.,51185272.,51453725.,
                                          51656025.,51812239.,51935969.,52036464.,52120047.,
                                          52191068.,52252540.,52306578.,52354680.,52397928.,
                                          52437115.,52472837.,52505551.,52535617.,52563325.])

        # Profiles generated by taking the min and max of 1.000.000 random parameter conf.
            if r0[-1] < 2.:  # Bodgy way to tell if pop under study is young
                ax.plot(r0,young_bary_Sig_min,color='green',lw=0.75,ls='dashed',dashes=(3,1))
                ax.plot(r0,young_bary_Sig_max,color='green',lw=0.75,ls='dashed',dashes=(3,1))
            else: # run is for old tracers, or both young and old tracers
                if r0[0] > 0.6:  # old pop only
                    ax.plot(r0,old_bary_Sig_min,color='green',lw=0.75,ls='dashed',dashes=(3,1))
                    ax.plot(r0,old_bary_Sig_max,color='green',lw=0.75,ls='dashed',dashes=(3,1))
                else:  # run has both old and young pop
                    ax.plot(r0,both_bary_Sig_min,color='green',lw=0.75,ls='dashed',dashes=(3,1))
                    ax.plot(r0,both_bary_Sig_max,color='green',lw=0.75,ls='dashed',dashes=(3,1))

            print ('Sig_baryon_model: r0=',r0)
            print ('*************************************************************')

        if prof == 'rho_baryon_vec':
            young_bary_rho_min = np.array([4188383.,3639628.,3188091.,2816620.,2510592.,
                                           2257762.,1981270.,1723154.,1500925.,1301532.,
                                           1132532.,991569.,873890.,775520.,693141.,
                                           623996.,565797.,516652.,474994.,439530.])
            young_bary_rho_max = np.array([9427400.,8322270.,7410273.,6645553.,6102458.,
                                           5627087.,5208404.,4837275.,4520961.,4244411.,
                                           3995199.,3793673.,3610437.,3442034.,3286276.,
                                           3154888.,3033907.,2920469.,2813602.,2712520.])

            old_bary_rho_min = np.array([3144525.,2364568.,1740064.,1266311.,945920.,
                                         731169.,583254.,480137.,406917.,353253.,
                                         312541.,280554.,254191.,231283.,211508.,
                                         194133.,178662.,164747.,152142.,140662.])
            old_bary_rho_max = np.array([7208423.,5779680.,4844702.,4193084.,3683998.,
                                         3278733.,2979109.,2720941.,2493819.,2291909.,
                                         2110820.,1946074.,1797384.,1664597.,1542381.,
                                         1429719.,1326909.,1232101.,1144345.,1063079.])

            both_bary_rho_min = np.array([4238843.,3073546.,2259546.,1638758.,1172239.,
                                          864858.,664870.,529982.,437570.,372272.,
                                          324471.,288141.,259490.,236116.,216484.,
                                          199063.,182910.,168533.,155617.,143932.])

            both_bary_rho_max = np.array([9415725.,7035723.,5717103.,4752183.,4066671.,
                                          3566511.,3161179.,2857129.,2615685.,2402181.,
                                          2210752.,2037551.,1879909.,1735866.,1603896.,
                                          1482762.,1371420.,1268970.,1174621.,1087668.])


            if r0[-1] < 2.:  # Bodgy way to tell if pop under study is young
                ax.plot(r0,young_bary_rho_min,color='green',lw=0.75,ls='dashed',dashes=(3,1))
                ax.plot(r0,young_bary_rho_max,color='green',lw=0.75,ls='dashed',dashes=(3,1))
            else: # run is for old tracers, or both young and old tracers
                if r0[0] > 0.6:  # old pop only
                    ax.plot(r0,old_bary_rho_min,color='green',lw=0.75,ls='dashed',dashes=(3,1))
                    ax.plot(r0,old_bary_rho_max,color='green',lw=0.75,ls='dashed',dashes=(3,1))
                else:  # run has both old and young pop
                    ax.plot(r0,both_bary_rho_min,color='green',lw=0.75,ls='dashed',dashes=(3,1))
                    ax.plot(r0,both_bary_rho_max,color='green',lw=0.75,ls='dashed',dashes=(3,1))
        return
    ## \fn fill_nice(self, ax, prof, pop, gp)
    # plot filled region for 1sigma and 2sigma confidence interval
    # @param ax axis object to plot into
    # @param prof string
    # @param pop int
    # @param gp


    def plot_full_distro(self, ax, prof, pop, gp): # Not used
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


    def plot_profile(self, basename, prof, pop, profile_source, gp):
        print ('In plot_profile: prof:',prof)
        plt.ioff()
        gh.LOG(1, 'plotting profile '+str(prof)+' for pop '+str(pop)+' in run '+basename)


        if prof != 'chi2':
            fig = plt.figure()
            ax  = fig.add_subplot(111)
            ax.set_xscale('log')

        if prof == 'chi2':  # SS: Could not find where x-axis is set to logscale
            fig = plt.figure()
            ax  = fig.add_subplot(111)
            ax.set_yscale('log')

        if prof == 'rho' or prof == 'Sig' or\
           prof == 'M' or prof == 'nu':
            fig = plt.figure()
            ax  = fig.add_subplot(111)
            ax.set_yscale('log')

        if prof in ['nu_vecs']:   # Log plots
            fig = plt.figure()
            ax  = fig.add_subplot(111)
            ax.set_xscale('linear') #Test usually linear
            ax.set_yscale('log') #SS: was linear before

        if gp.z_vec_Sigrho[-1] < 2.: # Bodgy way to determine if young pop
            rho_ylim = np.array([0,2.7e7]) # Sets ylim for all density components # TAGTAG
        else: # old pop, or both pops
            if gp.z_vec_Sigrho[0] > 0.6:  # old pop
                rho_ylim = np.array([0,2.9e7])
            else:  # both pops
                rho_ylim = np.array([0,2.2e7])

        if prof == 'rho_DM_vec':
            #fig = plt.figure()
            right_hand_side_label = True  # GeV/cm3 on right hand side or not
            if (right_hand_side_label):
                #fig = plt.figure(figsize=(2.3,2.3)) # Removed as it messes up font size
                fig = plt.figure()
                ax  = fig.add_subplot(111)
                #plt.subplots_adjust(left=0.21, right=0.79, bottom=0.22, top=0.8)
                # Above fits left and right labels but shrinks plot surface
                # left=0.3, right=0.95, bottom=0.25, top=0.9 makes DM plot surface & location 
                         # same as others, i.e. cuts away GeV/cm3 axis
                #plt.subplots_adjust(left=0.3, right=0.95, bottom=0.25, top=0.9)
                plt.subplots_adjust(left=0.1, right=0.75, bottom=0.25, top=0.9)
                ax2 = ax.twinx()
                ax2.set_ylabel('$\\rho_{\\rm{DM}}\\quad[\\rm{GeV}/\\rm{cm}^3]$')
                #ax_ylim = np.array([6,20])    # TAG
                ax_ylim = rho_ylim
                kpc = 3.0857E19   # kpc in m
                pc = 3.0857E16   # pc in m 
                Msun = 1.9891E30  # Sun's mass in kg
                GeV = 1.78266E-27 # GeV in kg
                ax_ylim2 = Msun/(GeV*(100.*kpc)**3)*ax_ylim
                #ax_ylim2 = 1E6*Msun/(GeV*(100.*kpc)**3)*ax_ylim
                #ax.set_ylim(ax_ylim)  # Not needed, set elsewhere
                ax2.set_ylim(ax_ylim2)
            else:
                fig = plt.figure()
                ax  = fig.add_subplot(111)
                #ax.set_ylim([0,14])    # TAG
            ax.set_xscale('linear')  # SS: changed rhoDM plots to linear scale
            ax.set_yscale('linear')  # plot this in log or linear scale...?


        if prof in ['kz_rho_DM_vec','kz_nu_vecs','sig_vecs', 'Sig_DM_vec',
                    'Sig_baryon_vec','Sig_total_vec', 'sigmaRz2_vecs',
                    'chi2_sigz2_vecs','chi2_sigRz2_vecs','rho_total_vec',
                    'sigz2_vecs','rho_baryon_vec']:  # Linear plots
            fig = plt.figure()
            ax  = fig.add_subplot(111)
            ax.set_xscale('linear')  # SS: was 'log' before
            ax.set_yscale('linear') #HS TEST 6 JAN

        self.plot_labels(ax, prof, pop, gp)

        if len(self.profs)>0:
            if prof == 'chi2':
                goodchi = []
                for k in range(len(self.profs)):
                    # do include all chi^2 values for plot
                    goodchi.append(self.chis[k])

                print('plotting profile chi for '+str(len(goodchi))+' models')
                print('min, max, maxsubset found chi2: ', min(self.chis), max(self.chis), self.subset[1])
                #chi_arr = np.array(self.chis)
                #print(chi_arr[np.where(chi_arr < 1.87)])

                bins, edges = np.histogram(np.log10(goodchi),\
                                           bins=max(6,np.sqrt(len(goodchi))),\
                                           density=True)
                ax.step(edges[1:], bins, where='pre')
                plt.draw()
                self.write_chi2(basename, edges, bins)

                fig.savefig(basename+'/output/' + profile_source + '_prof_chi2_0.pdf')
                plt.close(fig)
                return

            #Choice of plotting style
            #self.weighted_hist_heatmap(ax, prof, pop, gp)
            self.fill_nice(ax, prof, pop, gp)   # TAG

            #self.plot_N_samples(ax, prof, pop) #SS: Don't plot thin gray lines
            #self.plot_bins(ax, prof, pop, gp)  # Remove to not plot vertical bin lines
            if prof == 'Sig' or prof=='nu' or prof == 'sig' or prof == 'nu_vecs' or prof == 'sig_vecs' or prof == 'sigz2_vecs' or prof == 'sigmaRz2_vecs':
                #self.plot_data(ax, basename, prof, pop, gp)  # Plotting the data
                kfnnskdjncjnfksjdn = 1
            if gp.investigate == 'simplenu':
                self.plot_model_simplenu(ax, basename, prof, gp)

            if (gp.investigate == 'walk' or gp.investigate == 'gaia') \
               and (prof != 'Sig'):
                r0 = self.analytic.x0
                y0 = self.analytic.get_prof(prof, pop)
                ax.plot(r0, y0, 'b--', lw=2)

            #ax.set_xlim(0, gp.z_bincenters[-1]) #bodge
            #ax.set_xlim(0, np.round(gp.z_binmaxs[-1], 1))
            #ax.set_xlim(0., gp.z_all_pts[-1])

            #ax.set_xlim(0., max(np.concatenate(gp.z_binmax_vecs))) #Standard
            if prof == 'nu_vecs':   # Set nu x-range the same as other
                #ax.set_xlim(0., gp.nu_z_bincenter_vecs[0][-1]*1.3)
                ax.set_xlim(0.3, 2.7)
            else:
            #ax.set_xlim(0., gp.z_vec_Sigrho[-1]*1.2) # Tightening the plot TESTING !!
                if gp.z_vec_Sigrho[-1] < 2.: # Bodgy way to determine if young or old pop
                    #ax.set_xlim(0.5, 1.3)
                    #plt.xticks(np.array([0.6,0.8,1.0,1.2]))
                    ax.set_xlim(0.4, 1.4)
                    plt.xticks(np.array([0.4,0.6,0.8,1.0,1.2,1.4]))
                else: # old pop (or both pop)
                    ax.set_xlim(0.3, 2.7)
                    #ax.set_xlim(0.3, 1.4)  # TAGTAGTAG
            #plt.xticks(np.arange(0., gp.z_vec_Sigrho[-1]*1.3, 0.5)) # with 0.0 tick
            #plt.xticks(np.arange(0.5, gp.z_vec_Sigrho[-1]*1.3, 0.5)) # without 0.0 tick
            #plt.xticks(np.arange(0.5, gp.z_vec_Sigrho[-1]*1.2, 0.5)) # without 0.0 tick
            if prof in ['Sig_total_vec','Sig_baryon_vec','Sig_DM_vec']:
                if gp.z_vec_Sigrho[-1] < 2: # Bodgy way to determine if young pop
                    ax.set_ylim(0.,0.9e8 )   # TAGTAG  
                else: # old pop, or both pop
                    if gp.z_vec_Sigrho[0] > 0.6:  # old pop
                        ax.set_ylim(0.,1.5e8 )   
                    else:  # both pops in run
                        ax.set_ylim(0.,1.2e8 )   

                # Rescale y axis to get units in Msun/pc2 instead
                scale_y=1E-6
                ticks = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y*scale_y))
                ax.yaxis.set_major_formatter(ticks)

                yticks = ax.yaxis.get_major_ticks()  # TESTING !!!!
                yticks[0].label1.set_visible(False)
         
            #ax.set_xlim(0., 0.5)
            if prof in ['rho_total_vec','rho_baryon_vec','rho_DM_vec']:
                ax.set_ylim(rho_ylim ) # all rho have same ylim values

                # Rescale y axis to get units in Msun/pc3 instead
                scale_y=1E-9
                ticks = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y*scale_y))
                ax.yaxis.set_major_formatter(ticks)

            plt.minorticks_on()

            if prof == 'sigmaRz_vec' or prof == 'sigmaRz2_vecs': #TEST 6 Jan 16
                lknlfdnrkjfnkej = 1
                #ax.set_xlim(0., 0.25)

            #for tick in ax.xaxis.get_majorticklabels():
            #    tick.set_horizontalalignment("left")

            if prof == 'Sig' or prof=='nu' or prof == 'sig' or prof == 'nu_vecs' or prof == 'sig_vecs' or prof == 'sigz2_vecs' or prof == 'sigmaRz2_vecs':
                self.plot_data(ax, basename, prof, pop, gp)  # Plotting the data
           
            # Below line moves the leftmost element (0.0) of z axis to the right:
            #ax.xaxis.get_majorticklabels()[0].set_horizontalalignment("left")

            plt.draw()
        else:
            gh.LOG(1, 'empty self.profs')

        fig.savefig(basename+'/output/pdf/' + 'band_' + profile_source + '_prof_'+prof+'_'+str(pop)+'.pdf')
        plt.close(fig)
        #print ('End of plot_profile')
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

    def true_sigz2_func(self,binmin,binmax,nsmallbin,nbin,z0,Ntr,nu0,K,D,F,gp): # Not used
        print ('In true_sigz2_func. Think the sigz calculation below uses C in wrong way!!')
        pdb.set_trace()
        C = 22.85**2
        #C = 40.0**2
        truesig2_arr = np.zeros(nbin)
        for kbin in range(nbin):
            zvec = np.zeros(nsmallbin+1) # one extra point since lists both binmin & binmax
            zvec[0] = binmin[kbin]
            zvec[-1] = binmax[kbin]
            for k in range(1,nsmallbin+1):
                exptemp = math.exp(-zvec[k-1]/z0)-(math.exp(-zvec[0]/z0)-math.exp(-zvec[-1]/z0))/nsmallbin
                zvec[k] = -z0*math.log(exptemp)
            nuvec = nu0*np.exp(-zvec/z0)
            if gp.baryonmodel not in ['simplenu_baryon', 'trivial_baryon', 'obs_baryon', 'simplenu_baryon_gaussian']:  # assume DM only
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
        #Reconstruct number of tracers stars
        #Required later for some theory calculations.

        # BODGE: !!! TODO FIXME !!!
        #Ntr = np.round(sum(gp.dat.nu[0]*(gp.z_binmax_vecs[0] - gp.z_binmin_vecs[0])))

        #Baryon Model
        G1 = 4.299e-6   # Newton's constant in (km)^2*kpc/(Msun*s^2)
        K = 1500.
        F = 267.65
        D = 0.18

        #Set dark disc model
        plot_dd_data = False  # SS: Disabled DD option below
        D_dd = 2.5
        #if 'bdd' in gp.external_data_file[0]:
        #    K_dd = 900.
        #    plot_dd_data = True
        #elif 'dd' in gp.external_data_file[0]:
        #    K_dd = 300
        #    plot_dd_data = True

        #Set thick or thin disk
        #if 'simple2' in gp.external_data_file[0]:
        #    z0 = 0.9
        #elif 'simple' in gp.external_data_file[0]:
        #    z0 = 0.4
        z0 = 0.9  # SS: Set scale height of nu tracers explicitly instead

        #normC = 40.0
        normC = 22.85

        #Define z range over which to plot
        zvec=np.linspace(0, max(np.concatenate(gp.z_binmax_vecs)), 100)


        #Kzvec_total = -((K*zvec)/(np.sqrt(zvec**2 + D**2)) + 2.*F*zvec)

        #Define baryon density and force
        rho_z_baryon = (1/(4*np.pi*G1)) * abs((K*(D**2)/((D**2 + zvec**2)**(1.5))))
        Kzvec_baryon = -((K*zvec)/(np.sqrt(zvec**2 + D**2)))

        #Set DM density and force, adding dd if necessary
        rho_z_DM_const = (1/(4*np.pi*G1)) * abs(2.*F) * np.ones(len(zvec))
        Kzvec_const_DM = -(2.*F*zvec)

        if plot_dd_data:
            rho_z_dd = (1/(4*np.pi*G1)) * abs((K_dd*(D_dd**2)/((D_dd**2 + zvec**2)**(1.5))))
            rho_z_DM =  rho_z_DM_const + rho_z_dd
            Kzvec_DD = -((K_dd*zvec)/(np.sqrt(zvec**2 + D_dd**2)))
            Kzvec_DM = Kzvec_const_DM + Kzvec_DD
        else:
            rho_z_DM = rho_z_DM_const
            Kzvec_DM = Kzvec_const_DM

        #Calculate total density
        rho_z_total = rho_z_baryon + rho_z_DM
        #Sigma_z_total = (1000.**2)*abs(Kzvec_total)/(2*np.pi*4.299)
        Sigma_z_baryon = (1000.**2)*abs(Kzvec_baryon)/(2*np.pi*4.299)
        Sigma_z_DM = (1000.**2)*abs(Kzvec_DM)/(2*np.pi*4.299)
        Sigma_z_total = Sigma_z_baryon + Sigma_z_DM

        #k_z parameterization of rho total and baryon (unknown purpose HS04/12/15)
        k_z_rho_total = (3*(D**2)*K*zvec) / ( ((D**2 + zvec**2)**2.5) * ((K*(D**2))/((D**2 + zvec**2)**1.5) + 2*F))
        k_z_rho_baryon = (3*(D**2)*K*zvec) / ( ((D**2 + zvec**2)**2.5) * ((K*(D**2))/((D**2 + zvec**2)**1.5))) # not currently used


        # Backwards compatibility: if using old data, then all mass is in DM
        if gp.baryonmodel not in ['simplenu_baryon', 'trivial_baryon', 'obs_baryon', 'simplenu_baryon_gaussian']:
            gh.LOG(1, 'Simplenu Analytic: No baryon model, all mass is in DM.')
            Sigma_z_DM = Sigma_z_total
            rho_z_DM = rho_z_total


        pop=0
        ############
        ## Theory Curve for nu calculation
        #try:
        #   data_z_cut = gp.data_z_cut[pop]
        #except ValueError:
        #   data_z_cut = gp.data_z_cut

        #nu0 = Ntr/(z0*(1-np.exp(-data_z_cut/z0)))
        #nuvec = nu0*np.exp(-zvec/z0)

        #bincentermed, binmin, binmax, dum, dum = gh.readcol5(gp.files.nufiles[0])
        #nbin = np.size(bincentermed)
        #truen_arr = nu0*z0*(np.exp(-1.*binmin/z0)-np.exp(-1.*binmax/z0))  # true n.o. stars in bins
        #truenu_arr = truen_arr/(binmax-binmin)
        #if prof == 'nu_vec':
        #    ##ax.plot(zvec, nuvec, 'g-', alpha=0.5)
        #    #ax.plot(bincentermed, truenu_arr, 'g-', alpha=0.5)
        ############

        if prof == 'Sig_total_vec':
            lsdnsdlvjnaosfjsaz = 1 
            #ax.plot(zvec, Sigma_z_total, 'g-', alpha=0.5)
        #elif prof == 'Sig_baryon_vec':
            #print ('***GREEN LINE**********************************************')
            #ax.plot(zvec, Sigma_z_baryon, 'g-', alpha=0.5)
        #elif prof == 'Sig_DM_vec':
            #ax.plot(zvec, Sigma_z_DM, 'g-', alpha=0.5)

        #elif prof == 'rho_total_vec':
            #ax.plot(zvec, rho_z_total, 'g-', alpha = 0.5)
        elif prof == 'rho_baryon_vec':
            #Plot the prior range on baryons
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

            #ax.plot(zvec, rho_z_baryon, 'g-', alpha = 0.5)
            #ax.fill_between(zvec, rho_z_baryon_prior_min, rho_z_baryon_prior_max, color='r', alpha=0.1, lw=1)
            ##ax.plot(zvec, rho_z_baryon_prior_max,'g-', alpha = 0.5, linewidth=1)
            ##ax.plot(zvec, rho_z_baryon_prior_min,'g-', alpha = 0.5, linewidth=1)


        elif prof == 'rho_DM_vec':
            #rescale_prof = 10.**6
            rescale_prof = 1.
            #ax.plot(zvec, rho_z_DM/rescale_prof, 'g--', dash_joinstyle='round', dashes=(6,2), alpha = 1.0, lw=2)

        # if all mass is described by DM, then plot kz_rho_DM_vec
        # if the simplenu_baryon model is used, then kz_rho_DM_vec should equal zero
        elif prof == 'kz_rho_DM_vec' and gp.baryonmodel in ['none', 'sim']: #i.e. if all mass is described by DM
            ax.plot(zvec, k_z_rho_total, 'g.=', alpha = 0.5)


        elif prof == 'sigz2_vec':
            ##sigz2_analytic = phys.sigz2(zvec, Sigma_z_total, nuvec, (normC**2)*nuvec[0])
            ##ax.plot(zvec, sigz2_analytic, 'g-', alpha = 0.5)

            # BODGE: !!! TODO FIXME !!!
            #true_sigz2_arr = self.true_sigz2_func(binmin,binmax,1000,nbin,z0,Ntr,nu0,K,D,F,gp)
            sjnjnsdkjn = 1

            #ax.plot(bincentermed, true_sigz2_arr, 'g-', alpha = 0.5)
            # Uncomment above line to plot theory curve for sigz2

        elif prof == 'sigmaRz_vec' or prof == 'sigmaRz2_vecs':
            try:
                positive_sigRz_data_sign = gp.positive_sigRz_data_sign
            except ValueError:
                positive_sigRz_data_sign = False

            if positive_sigRz_data_sign:
                A = 180.08
            else:
                A = -180.08
            n0 = 1.44 #0 to avoid conflict with PDB command

            sigmaRz2 = A*(zvec)**n0 #z [kpc]
            #ax.plot(zvec, sigmaRz2, 'g-', alpha = 0.5)

            #print ('zvec:',zvec)
            #print ('sigmaRz2 green line:',sigmaRz2)
            # Above: plot of green line with true values

        return

    def plot_bins(self, ax, prof, pop, gp):  # Used (not anymore)
        #print ('In plot_bins, prof=',prof)
        if prof in ['kz_rho_DM_vec', 'rho_DM_vec', 'Sig_DM_vec', 'rho_baryon_vec',
            'Sig_baryon_vec', 'rho_total_vec', 'Sig_total_vec']:
            ##r0 = gp.z_all_pts_sorted
            r0 = gp.z_vec_Sigrho  # REMOVE these vertical lines !?!?
        elif prof in ['nu_vecs', 'kz_nu_vecs']:
            r0 = gp.nu_z_bincenter_vecs[pop]
        elif prof in ['z_vecs', 'sigz2_vecs', 'tilt_vecs',
            'sigmaRz2_vecs', 'chi2_sigz2_vecs', 'chi2_sigRz2_vecs']:
            r0 = gp.z_bincenter_vecs[pop]

        #[ax.axvline(x, color='c', linewidth=0.1) for x in r0]
        [ax.axvline(x, color='gray', linewidth=0.1) for x in r0]
        return


    def __repr__(self):
        return "Profile Collection with "+str(len(self.profs))+" Profiles"
    ## \fn __repr__(self)
    # string representation for ipython
