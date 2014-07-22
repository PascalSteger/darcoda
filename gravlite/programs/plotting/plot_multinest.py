#!/usr/bin/env python3

## @file
# plot all marginals

# (c) 2013 ETHZ Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import pdb


def read_scale(basename, gp):
    gp.Rscale = []
    gp.Nu0pc = []
    gp.maxvlos = []
    for i in range(gp.pops+1):
        # basename + 'scale_' + str(i) + '.txt'
        A = np.loadtxt(gp.files.get_scale_file(i),\
                       unpack=False, skiprows=1)
        gp.Rscale.append(A[0])
        gp.Nu0pc.append(A[2])
        gp.maxvlos.append(A[4])
## \fn read_scale(basename, gp)
# read scale file, store into gp.*scale
# @param basename string
# @param gp


def correct_E_error(filename):
    import os
    os.system("sed -i 's/\\([0-9]\\)-\\([0-9]\\)/\\1E-\\2/g' " + filename)
    return
## \fn correct_E_error(filename)
# replace 3.4230210-301 with 3.4230210E-301 the sed way
# @param filename string


def read_models(basename):
    # read in all accepted models
    print(basename+'/{ev.dat, phys_live.points}')
    correct_E_error(basename + '/ev.dat')
    correct_E_error(basename + '/phys_live.points')
    REJECTED = np.loadtxt(basename+'/ev.dat', skiprows=0, unpack=False)
    LIVE = np.loadtxt(basename+'/phys_live.points', skiprows=0, unpack=False)
    ALL = np.vstack([REJECTED[:,:-3], LIVE[:,:-2]])
    # for debugging, based on live points
    # ALL = np.vstack([LIVE[:,:-2]])
    # for debugging, fixed 100 models after N=1000 iterations
    # ALL = np.vstack([REJECTED[1000:1100,:-2]])
    return ALL
## \fn read_models(basename)
# read in all models, concatenate them
# @param basename string


def fill_nice(r0, M95lo, M68lo, Mmedi, M68hi, M95hi, gp):
    import gl_plot as gpl
    gpl.fill_between(r0, M95lo, M95hi, color='black', alpha=0.2, lw=0.1)
    gpl.plot(r0, M95lo, color='black', lw=0.4)
    gpl.plot(r0, M95hi, color='black', lw=0.3)
    gpl.fill_between(r0, M68lo, M68hi, color='black', alpha=0.4, lw=0.1)
    gpl.plot(r0, M68lo, color='black', lw=0.4)
    gpl.plot(r0, M68hi, color='black', lw=0.3)
    gpl.plot(r0, Mmedi, 'r', lw=1)
    for i in range(gp.pops):
        gpl.axvline(gp.Rscale[i+1], color='gray', lw=0.5) # [pc]
    gpl.xlim([r0[0], r0[-1]])
    return
## \fn fill_nice(r0, M95lo, M68lo, Mmedi, M68hi, M95hi, gp)
# plot filled region for 1sigma and 2sigma confidence interval
# @param r0 radius in [pc]
# @param M95lo lower 2sigma bound
# @param M68lo lower 1sigma bound
# @param Mmedi median
# @param M68hi upper 1sigma bound
# @param M95hi upper 2sigma bound
# @param gp


def physical(r0, prof, pop, tmp_rho, tmp_nu, tmp_beta, gp):
    import gl_physics as phys
    from gl_project import rho_INT_Rho, rho_param_INT_Rho
    if prof == 'rho':
        tmp_prof = phys.rho(r0, tmp_rho, gp)
    elif prof == 'nr':
        tmp_prof = tmp_rho[1:]
    elif prof == 'betastar':
        tmp_prof = phys.mapping_beta_star_poly(r0, tmp_beta, gp)
    elif prof == 'beta':
        tmp_prof = phys.beta(r0, tmp_beta, gp)
    elif prof == 'Sig':
        tmp_prof = rho_param_INT_Rho(r0, tmp_nu, gp)
    elif prof == 'sig':
        tmp_sig, tmp_kap = phys.sig_kap_los(r0, pop, tmp_rho, tmp_nu, tmp_beta, gp)
        tmp_prof = tmp_sig
    elif prof == 'kap':
        tmp_sig, tmp_kap = phys.sig_kap_los(r0, pop, tmp_rho, tmp_nu, tmp_beta, gp)
        tmp_prof = tmp_kap
    return tmp_prof
## \fn physical(r0, prof, pop, tmp_rho, tmp_nu, tmp_beta, gp)
# calculate physical representations of profiles
# @param r0 radius in [pc]
# @param prof string of which profile to calculate
# @param pop  component (0,1,2) for (overall, 1st stellar tracer, 2nd stellar tracer)
# @param tmp_rho  internal representation of rho
# @param tmp_nu   internal representation of nu
# @param tmp_beta internal representation of beta
# @param gp

def plot_labels(prof, pop):
    import gl_plot as gpl
    gpl.xlabel('$R\\quad[\\rm{pc}]$')
    if prof == 'rho':
        gpl.ylabel('$\\rho\\quad[\\rm{M}_\\odot/\\rm{pc}^3]$')
    elif prof == 'nr':
        gpl.ylabel('$n(r)$')
    elif prof == 'beta':
        gpl.ylabel('$\\beta_'+str(pop)+'$')
        gpl.ylim([-1,1])
    elif prof == 'betastar':
        gpl.ylabel('$\\beta^*_'+str(pop)+'$')
        gpl.ylim([-1.05,1.05])
    elif prof == 'Sig':
        gpl.ylabel('$\\Sigma_'+str(pop)+'\\quad[\\rm{M}_\\odot/\\rm{pc}^2]$')
    elif prof == 'sig':
        analytic = 0.*gp.xipol
        gpl.ylabel('$\\sigma_{\\rm{LOS},'+str(pop)+'}\quad[\\rm{km}/\\rm{s}]$')
    return
## \fn plot_labels(prof)
# draw x and y labels
# @param prof string of profile to look at
# @param pop which population to analyze: 0 (all), 1, 2, ...


def ipol_rhalf_log(X, Y, rhalf):
    for i in range(len(X)):
        if X[i] < rhalf:
            lowr = X[i]
            lowi = i
        if X[len(X)-i-1] > rhalf:
            highr = X[len(X)-i-1]
            highi = len(X)-i-1
    dr = np.log(highr) - np.log(lowr)
    dy = np.log(Y[highi]) - np.log(Y[lowi])

    Yhalf = np.exp(np.log(Y[lowi]) + (np.log(rhalf)-np.log(lowr))*dy/dr)
    return Yhalf
# \fn ipol_rhalf_log(X, Y, rhalf)
# interpolate function in loglog space to value in between data points,
# used for determination of Sigma(rhalf)
# @param X radii, usually gp.xipol in our case, in [pc]
# @param Y values, usually Sigma(gp.xipol) here, in [arb. units]
# @param rhalf half-light radius in [pc]


def plot_data(r0, prof, pop, output):
    output.add('radius (data) [pc]', r0)
    if prof == 'Sig':
        # get 2D data here
        # Rbin [Rscale], Binmin [Rscale], Binmax [Rscale], Nu(R)/Nu(0) [1], error [1]
        dum,dum,dum,Nudat,Nuerr = np.transpose(np.loadtxt(gp.files.nufiles[pop], unpack=False, skiprows=1))
        Nudat *= gp.Nu0pc[pop]
        Nuerr *= gp.Nu0pc[pop]
        Nunorm = ipol_rhalf_log(gp.xipol, Nudat, gp.Rscale[pop])
        Nudat /= Nunorm; Nuerr /= Nunorm
        dat = Nudat; err = Nuerr
        output.add('data [Msun/pc^2]', dat)
        output.add('error [Msun/pc^2]', err)
        output.add('data - error [Msun/pc^2]', dat-err)
        output.add('data + error [Msun/pc^2]', dat+err)
    elif prof == 'sig':
        DATA = np.transpose(np.loadtxt(gp.files.sigfiles[pop], unpack=False, skiprows=1))
        sigdat = DATA[4-1];        sigerr = DATA[5-1]
        sigdat *= gp.maxvlos[pop]; sigerr *= gp.maxvlos[pop]
        dat = sigdat; err = sigerr
        output.add('data [km/s]', dat)
        output.add('error [km/s]', err)
        output.add('data - error [km/s]', dat-err)
        output.add('data + error [km/s]', dat+err)
    import gl_plot as gpl
    gpl.fill_between(r0, dat-err, dat+err, color='blue', alpha=0.2, lw=1)
    return output
## \fn plot_data(r0, prof, pop)
# plot data as blue shaded region
# @param r0 radii [pc]
# @param prof string of profile
# @param pop which population to analyze: 0 (all), 1, 2, ...




def unit(prof):
    if prof == 'rho':
        unit = '[Msun/pc^3]'
    elif prof == 'nr':
        unit = '[log(Msun/pc^3)/log(pc)]'
    elif prof == 'beta':
        unit = '[1]'
    elif prof == 'betastar':
        unit = '[1]'
    elif prof == 'Sig':
        unit = '[(Msun/pc^3)/Sigma(r_1/2)]'
    elif prof == 'sig':
        unit = '[km/s]'
    return unit
## \fn unit(prof)
# return string with units for profile
# @param prof string for profile


def plot_analytic(r0, prof, pop, gp, output):
    output.add('radius (analytic) [pc]', r0)
    from gl_project import rho_INT_Rho, rho_param_INT_Rho
    import gl_analytic as ga
    if gp.investigate == 'gaia':
        if prof == 'rho':
            analytic = ga.rhogaiatot_3D(r0)
        elif prof == 'nr':
            analytic = ga.nrgaiatot_3D(r0)
        elif prof == 'beta':
            analytic = ga.betagaia(r0)
        elif prof == 'betastar':
            b = ga.betagaia(r0)
            analytic = b/(2.-b)
        elif prof == 'Sig':
            analytic_3D = ga.nugaiatot_3D(r0)
            analytic = rho_INT_Rho(r0, analytic_3D)
            Nunorm = ipol_rhalf_log(r0, analytic_3D, gp.Rscale[pop])
            analytic /= Nunorm
        elif prof == 'sig':
            analytic = 0.*r0 # TODO: find analytical profile
    elif gp.investigate == 'walk':
        if prof == 'rho':
            analytic = ga.rhowalk_3D(r0, gp)[pop]
        elif prof == 'nr':
            analytic = ga.nrwalktot_3D_deriv(r0, gp)
            # TODO: enable for populations 1, 2
        elif prof == 'beta':
            analytic = ga.walker_delta(pop, r0, gp)
        elif prof == 'betastar':
            b = ga.walker_delta(pop, r0, gp)
            analytic = b/(2.-b)
        elif prof == 'Sig':
            analytic_3D = ga.rhowalk_3D(r0, gp)[pop]
            analytic = rho_INT_Rho(r0, analytic_3D)
            Nunorm = ipol_rhalf_log(r0, analytic, gp.Rscale[pop])
            analytic /= Nunorm
        elif prof == 'sig':
            analytic = 0.*r0            # TODO
    output.add('analytic ' + unit(prof), analytic)
    import gl_plot as gpl
    gpl.plot(r0, analytic, 'b')
    return output
## \fn plot_analytic(r0, prof, pop, gp)
# plot analytic curves in blue
# @param r0 radius in [pc]
# @param prof string of profile
# @param pop which population to analyze (attention, default is 0 here!)


def run(basename, prof, pop, output, gp):
    read_scale(basename, gp)
    # read in half light radi(i,us), Nu0[munit/Rscale^2], Nu0[munit/pc^2], munit
    print('r_{1/2, all} = ', gp.Rscale[0])
    # for all tracer particles, stored first ("0")
    # read in radial bins
    import gl_helper as gh
    Radii, Binmin, Binmax, Nudat1, Nuerr1 = gh.readcol5(gp.files.nufiles[pop])
    gp.xipol = Radii*gp.Rscale[0]       # [pc]
    maxR = max(Radii)
    Radii = np.hstack([Radii, 2*maxR, 4*maxR, 8*maxR])
    gp.xepol = Radii*gp.Rscale[1]
    
    models = read_models(basename)
    
    # use physical representation for profiles
    profs = []
    for i in range(len(models)):
        tmp_rho = []; tmp_nu = []; tmp_beta = []
        M = models[i]
        tmp_rho = M[0:gp.nepol]
        off = gp.nepol
        for comp in range(gp.pops):
            tmp_nu.append(M[off:off+gp.nepol])
            off += gp.nepol
            tmp_beta.append(M[off:off+gp.nbeta])
            off += gp.nbeta
        try:
            tmp_prof = physical(gp.xepol, prof, pop, tmp_rho, tmp_nu[0], tmp_beta[0], gp)
            profs.append(tmp_prof)
        except Exception as detail:
            print('handling error in profile: ', detail)
            continue
    sortedprofs = gh.sort_profiles_binwise(np.array(profs).transpose())

    # get 1,2 sigma ranges, and median
    Mmin, M95lo, M68lo, Mmedi, \
      M68hi, M95hi, Mmax = gh.get_median_1_2_sig(sortedprofs)
    if prof == 'Sig':
        Nunorm = ipol_rhalf_log(gp.xepol, Mmedi, gp.Rscale[0])
        Mmin  /= Nunorm;   M95lo /= Nunorm;   M68lo /= Nunorm
        Mmedi /= Nunorm
        M68hi /= Nunorm;   M95hi /= Nunorm;   Mmax  /= Nunorm
        
    # set up plotting area
    import gl_plot as gpl
    fig = gpl.start(); gpl.xscale('log')
    if prof=='rho' or prof=='Sig':
        gpl.yscale('log')

    plot_labels(prof, pop)
    if prof=='nr':
        gpl.xscale('log')
        # gpl.xlim([0.,2.*gp.xipol[-1]])
        gp.xepol = np.hstack([gp.xepol[0]/2., gp.xepol, gp.rinfty*gp.xepol[-1]])

    if prof == 'Sig' or prof=='sig':
        output = plot_data(gp.xipol, prof, pop, output)
    output = plot_analytic(gp.xepol, prof, pop, gp, output)

    radfill = gp.xepol
    if prof == 'sig' or prof=='Sig':
        radfill = gp.xipol
    uni = unit(prof)
    output.add('radius (models) [pc]', radfill)
    output.add('M 95% CL low ' + uni, M95lo)
    output.add('M 68% CL low ' + uni, M68lo)
    output.add('M median '     + uni, Mmedi)
    output.add('M 68% CL high '+ uni, M68hi)
    output.add('M 95% CL high '+ uni, M95hi)
    
    fill_nice(radfill, M95lo, M68lo, Mmedi, M68hi, M95hi, gp)
    output.write(basename+'/prof_'+prof+'_'+str(pop)+'.ascii')
    gpl.savefig(basename+'/prof_'+prof+'_'+str(pop)+'.png')

    # save as pdf
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(basename+'/prof_'+prof+'_'+str(pop)+'.pdf')
    pp.savefig(fig)
    pp.close()
    
    gpl.ion(); gpl.show(); gpl.ioff()
    return
## \fn run(basename, prof, pop, gp)
# main functionality
# @param basename
# @param prof string of profile to look at
# @param pop population number
# @param gp



def main():
    import select_run
    timestamp, basename, prof, pop = select_run.run()

    # include runtime gl_params, but other files all from general file
    import sys
    import import_path as ip
    ip.insert_sys_path(basename)
    import gl_params
    del sys.path[0] # delete first entry in sys.path (output dir, only used for gl_params)
    
    gp = gl_params.Params(timestamp) # set geometry here (update gl_params in old output dirs)

    import gl_file
    gp.dat = gl_file.get_data(gp)

    gp.pops = select_run.get_pops(basename)
    print('working with ', gp.pops, ' populations')

    import gl_output as go
    output = go.Output()

    run(basename, prof, pop, output, gp)

    import gl_plot as gpl
    gpl.show()

if __name__ == '__main__':
    main()
