#!/usr/bin/env python3

## @file
# plot all marginals

# (c) 2013 ETHZ Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import pdb
import gl_params as gp
import gl_helper as gh
import gl_plot   as gpl
import gl_analytic as ga
import gl_physics as phys
from gl_project import rho_INT_Rho

def read_scale():
    gp.Rscale = []
    gp.Nu0rscale = [];  gp.Nu0pc = []
    gp.totmass = [];    gp.maxvlos = []
    for i in range(gp.pops+1):
        A = np.loadtxt(gp.files.get_scale_file(i), unpack=False, skiprows=1)
        gp.Rscale.append(A[0])
        gp.Nu0rscale.append(A[1])
        gp.Nu0pc.append(A[2])
        gp.totmass.append(A[3])
        gp.maxvlos.append(A[4])
## \fn read_scale()
# read scale file, store into gp.*scale


def read_models():
    # read in all accepted models
    print(basename+'/{ev.dat, phys_live.points}')
    REJECTED = np.loadtxt(basename+'/ev.dat', skiprows=0, unpack=False)
    LIVE = np.loadtxt(basename+'/phys_live.points', skiprows=0, unpack=False)
    # ALL = np.vstack([REJECTED[:,:-3], LIVE[:,:-2]])
    # for debugging. TODO: comment out
    ALL = np.vstack([REJECTED[-4:,:-3], LIVE[-2:,:-2]])
    return ALL
## \fn read_models
# read in all models, concatenate them


def fill_nice(r0, M95lo, M68lo, Mmedi, M68hi, M95hi):
    gpl.fill_between(r0, M95lo, M95hi, color='black', alpha=0.2, lw=0.1)
    gpl.plot(r0, M95lo, color='black', lw=0.4)
    gpl.plot(r0, M95hi, color='black', lw=0.3)
    gpl.fill_between(r0, M68lo, M68hi, color='black', alpha=0.4, lw=0.1)
    gpl.plot(r0, M68lo, color='black', lw=0.4)
    gpl.plot(r0, M68hi, color='black', lw=0.3)
    gpl.plot(r0, Mmedi, 'r', lw=1)
    gpl.axvline(gp.rstarhalf, color='gray', lw=0.5) # [pc]
    gpl.xlim([r0[0], r0[-1]])
    return
## \fn fill_nice(r0, M95log, M68lo, Mmedi, M68hi, M95hi)
# plot filled region for 1sigma and 2sigma confidence interval
# @param r0 radius in [pc]
# @param M95lo
# @param M68lo
# @param Mmedi
# @param M68hi
# @param M95hi


def physical(r0, prof, pop, tmp_rho, tmp_nu, tmp_beta):
    if prof == "rho":
        tmp_prof = phys.rho(r0, tmp_rho)
    elif prof == 'nr':
        tmp_prof = tmp_rho[1:]
    elif prof == "nu":
        tmp_prof = rho_INT_Rho(r0, phys.rho(r0, tmp_nu))
    elif prof == "betastar":
        tmp_prof = phys.mapping_beta_star_poly(r0, tmp_beta)
    elif prof == "beta":
        tmp_prof = phys.beta(r0, tmp_beta)
    elif prof == "sig":
        tmp_sig, tmp_kap = phys.sig_kap_los(r0, pop, tmp_rho, tmp_nu, tmp_beta)
        tmp_prof = tmp_sig
    elif prof == "kap":
        tmp_sig, tmp_kap = phys.sig_kap_los(r0, pop, tmp_rho, tmp_nu, tmp_beta)
        tmp_prof = tmp_kap
    return tmp_prof
## \fn physical(r0, prof, pop, tmp_rho, tmp_nu, tmp_beta)
# calculate physical representations of profiles
# @param r0 radius in [pc]
# @param prof string of which profile to calculate
# @param pop  component (0,1,2) for (overall, 1st stellar tracer, 2nd stellar tracer)
# @param tmp_rho internal representation of rho
# @param tmp_nu
# @param tmp_beta


def plot_labels(prof, pop):
    gpl.xlabel('$R\\quad[\\rm{pc}]$')
    if prof=="rho":
        gpl.ylabel('$\\rho\\quad[\\rm{M}_\\odot/\\rm{pc}^3]$')
    elif prof=="nr":
        gpl.ylabel('$n(r)$')
    elif prof=="beta":
        gpl.ylabel('$\\beta_'+str(pop)+'$')
        gpl.ylim([-5,1])
    elif prof=="betastar":
        gpl.ylabel('$\\beta^*_'+str(pop)+'$')
        gpl.ylim([-1,1])
    elif prof=="nu":
        gpl.ylabel('$\\Sigma_'+str(pop)+'\\quad[\\rm{M}_\\odot/\\rm{pc}^2]$')
    elif prof=="sig":
        analytic = 0.*gp.xipol
        gpl.ylabel('$\\sigma_{\\rm{LOS},'+str(pop)+'}\quad[\\rm{km}/\\rm{s}]$')
    return
## \fn plot_labels(prof)
# draw x and y labels
# @param prof string of profile to look at
# @param pop which population to analyze: 0 (all), 1, 2, ...


def plot_data(r0, prof, pop):
    if prof=='nu':
        # get 2D data here
        # Rbin [Rscale], Binmin [Rscale], Binmax [Rscale], Nu(R)/Nu(0) [1], error [1]
        DATA_2D = np.transpose(np.loadtxt(gp.files.nufiles[pop], unpack=False, skiprows=1))
        Nudat = DATA_2D[4-1]; Nuerr = DATA_2D[5-1]
        Nudat *= gp.Nu0pc[pop]; Nuerr *= gp.Nu0pc[pop]
        dat = Nudat; err = Nuerr
    elif prof=='sig':
        DATA = np.transpose(np.loadtxt(gp.files.sigfiles[pop], unpack=False, skiprows=1))
        sigdat = DATA[4-1];        sigerr = DATA[5-1]
        sigdat *= gp.maxvlos[pop]; sigerr *= gp.maxvlos[pop]
        dat = sigdat; err = sigerr
    gpl.fill_between(r0, dat-err, dat+err, color='blue', alpha=0.2, lw=1)
    return
## \fn plot_data(r0, prof, pop)
# plot data as blue shaded region
# @param r0 radii [pc]
# @param prof string of profile
# @param pop which population to analyze: 0 (all), 1, 2, ...


def plot_analytic(r0, prof, pop=0):
    if gp.investigate == 'gaia':
        if prof=="rho":
            analytic = ga.rhogaiatot_3D(r0)
        elif prof=="nr":
            analytic = ga.nrgaiatot_3D(r0)
        elif prof=="beta":
            analytic = ga.betagaia(r0)
        elif prof=="betastar":
            b = ga.betagaia(r0)
            analytic = b/(2.-b)
        elif prof=="nu":
            analytic_3D = ga.nugaiatot_3D(r0)
            analytic = rho_INT_Rho(r0, analytic_3D)
        elif prof=="sig":
            analytic = 0.*r0 # TODO
    elif gp.investigate == 'walk':
        if prof=='rho':
            analytic = ga.rhowalk_3D(r0)[pop]
        elif prof=='nr':
            analytic = ga.rhowalk_3D(r0)[pop] # TODO
        elif prof=='beta':
            analytic = ga.walker_delta(pop, r0) # TODO
        elif prof=='betastar':                  # TODO
            b = ga.walker_delta(pop, r0)
            analytic = b/(2.-b)
        elif prof=='nu':                # TODO: mass ratio M_tracer/M_DM
            pdb.set_trace()
            analytic_3D = ga.rhowalk_3D(r0)[pop]
            analytic = rho_INT_Rho(r0, analytic_3D)
        elif prof=='sig':
            analytic = 0.*r0            # TODO
    gpl.plot(r0, analytic, 'b')
    return analytic
## \fn plot_analytic(r0, prof, pop=0)
# plot analytic curves in blue
# @param r0 radius in [pc]
# @param prof string of profile
# @param pop which population to analyze (attention, default is 0 here!)


def run(basename, prof):
    read_scale()

    # read in half light radius, Nu0[munit/Rscale^2], Nu0[munit/pc^2], munit
    print(gp.files.get_scale_file(0))
    A = np.loadtxt(gp.files.get_scale_file(0), unpack=False, skiprows=1)
    gp.rstarhalf = A[0] # halflight radius (all tracer part.s), stored first ("0")
    # read in radial bins
    radii, binmin, binmax, nudat1, nuerr1 = gh.readcol5(gp.files.nufiles[0])
    gp.xipol = radii*gp.rstarhalf       # [pc]
    maxr = max(radii)
    radii = np.hstack([radii, 2*maxr, 4*maxr, 8*maxr])
    gp.xepol = radii*gp.rstarhalf
    
    models = read_models()
    
    # use physical representation for profiles
    profs = []
    for i in range(len(models)):
        M = models[i]
        tmp_rho = M[0:gp.nepol]
        off = gp.nipol
        tmp_nu = M[off:off+gp.nepol]
        off += gp.nipol
        tmp_beta = M[off:off+gp.nbeta]
        # TODO: include handling of N populations

        try:
            tmp_prof = physical(gp.xepol, prof, pop, tmp_rho, tmp_nu, tmp_beta)
            profs.append(tmp_prof)
        except Exception as detail:
            print('handling error in profile', detail)
            continue

    sortedprofs = gh.sort_profiles_binwise(np.array(profs).transpose())

    Mmin, M95lo, M68lo, Mmedi, \
      M68hi, M95hi, Mmax = gh.get_median_1_2_sig(sortedprofs)

    # plot ranges
    if prof=='rho' or prof=='nu':
        gpl.startlog(); gpl.xscale('log')
    else:
        gpl.start()
    plot_labels(prof, pop)
            
    if prof=='nr':
        gpl.xscale('log')
        # gpl.xlim([0.,2.*gp.xipol[-1]])
        gp.xipol = np.hstack([gp.xipol[0]/2., gp.xipol, gp.rinfty*gp.xipol[-1]])
    if prof=='nu' or prof=='sig':
        plot_data(gp.xipol, prof, pop)
    plot_analytic(gp.xepol, prof, pop)

    radfill = gp.xepol
    if prof == 'sig':
        radfill = gp.xipol
    fill_nice(radfill, M95lo, M68lo, Mmedi, M68hi, M95hi)

    gpl.savefig(basename+'/prof_'+prof+'_'+str(pop)+'.png')
    gpl.ion(); gpl.show(); gpl.ioff()
    return
## \fn run()
# main functionality

if __name__=="__main__":
    import select_run
    basename, prof, pop = select_run.run()
    run(basename, prof)
    gpl.show()
