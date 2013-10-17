#!/usr/bin/env ipython

##
# @file
# all functions to work with plots
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import pdb
from pylab import *
from matplotlib.ticker import MaxNLocator

import gl_params as gp
if gp.geom == 'sphere':
    import physics_sphere as phys
else:
    import physics_disc as phys

import gl_funs as gfun
from gl_analytic import *
from gl_int import *
from gl_project import rho_INT_Rho
import gl_units as units

global f, axs

## start new figure. to be used in debug mode
def start():
    clf()
    yscale('log')
    return

## prepare output plots, defining geometry
def prepare_plots():
    global f,axs
    if not gp.showplot: 
        return
    ion()
    f = figure(0)
    ax1 = f.add_subplot(421)                         # nu1
    ax2 = f.add_subplot(422, sharex=ax1, sharey=ax1) # nu2
    ax3 = f.add_subplot(423, sharex=ax1)             # sig1
    ax4 = f.add_subplot(424, sharex=ax1, sharey=ax3) # sig2
    ax5 = f.add_subplot(425, sharex=ax1) # kurtosis [sphere] or kappa_i [disc]
    ax6 = f.add_subplot(426, sharex=ax1, sharey=ax5) 
    # ^-- kurtosis [sphere] or empty [disc]
    ax7 = f.add_subplot(427, sharex=ax1) # rho or M
    ax8 = f.add_subplot(428, sharex=ax1) # delta = beta [sphere] or tilt [disc]
    axs=[[ax1, ax2],[ax3,ax4],[ax5,ax6],[ax7,ax8]]

    # define visibility of labels
    setp(ax1.get_xticklabels(), visible=False)
    setp(ax2.get_xticklabels(), visible=False)
    setp(ax3.get_xticklabels(), visible=False)
    setp(ax4.get_xticklabels(), visible=False)
    setp(ax5.get_xticklabels(), visible=False)
    setp(ax6.get_xticklabels(), visible=False)

    setp(ax2.get_yticklabels(), visible=False)
    setp(ax4.get_yticklabels(), visible=False)
    setp(ax6.get_yticklabels(), visible=False)
    setp(ax8.get_yticklabels(), visible=False)
    draw()
    return

## set limits on axis ax
def setlims(ax,xlim,ylim):
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    return

## plot first bit of data with green color, including 1sigma error bars
def plot_data():
    if not gp.showplot: return

    # plot nu
    axs[0][0].cla()
    x = gp.ipol.nux1_2D                 # [pc]
    if gp.geom == 'disc':
        x = gp.ipol.nux1
    if gp.analytic: 
        axs[0][0].plot(x, Sigma_anf(x), c='blue', lw=1)
    lbound = gfun.compare_nu(1,True,False) - gfun.compare_nu(1,True,True) 
    # [munit/pc**2]
    # gives gp.ipol.nudat1_2D
    ubound = gfun.compare_nu(1,True,False) + gfun.compare_nu(1,True,True) 
    # [munit/pc**2]
    axs[0][0].fill_between(x,lbound,ubound,alpha=0.5,color='g') 
    # [pc], 2*[munit/pc**2]
    if gp.lograd:
        axs[0][0].set_xscale('log')
    axs[0][0].set_yscale('log')
    blim1 = [min(lbound)/1.5, 1.5*max(ubound)]
    setlims(axs[0][0],[0,max(x)],blim1)
    setlims(axs[0][1],[0,max(x)],blim1)
    axs[0][0].set_ylabel('$\\nu_i\\quad[\\rm{M}_\odot/\\rm{pc}^2]$')
    axs[0][0].xaxis.set_major_locator(MaxNLocator(4))
    
    # plot sigma_LOS 1
    axs[1][0].cla()
    if gp.analytic: axs[1][0].plot(gp.ipol.sigx1, sig_los_anf(gp.ipol.sigx1),\
                                   c='blue', lw=1)
    x = gp.ipol.sigx1                                         # [pc]
    lbound = (gp.ipol.sigdat1 - gp.ipol.sigerr1)              # [km/s]
    ubound = (gp.ipol.sigdat1 + gp.ipol.sigerr1)              # [km/s]
    axs[1][0].fill_between(x,lbound,ubound,alpha=0.5,color='green') # [pc],2*[km/s]
    blimsig1 = [0,1.5*max(ubound)]
    setlims(axs[1][0],[0,max(x)],blimsig1)
    setlims(axs[1][1],[0,max(x)],blimsig1)
    axs[1][0].set_ylabel('$\\sigma_{i,\\rm{LOS}}\\quad[\\rm{km/s}]$')
    axs[1][0].yaxis.set_major_locator(MaxNLocator(4))
    axs[1][0].xaxis.set_major_locator(MaxNLocator(4))

    # plot kappa_LOS 1
    axs[2][0].cla()
    x = gp.ipol.kapx1                            # [pc]
    lbound = (gp.ipol.kapdat1 - gp.ipol.kaperr1) # [1]
    ubound = (gp.ipol.kapdat1 + gp.ipol.kaperr1) # [1]
    axs[2][0].fill_between(x,lbound,ubound,alpha=0.5,color='green') # [pc], 2*[1]
    setlims(axs[2][0],[0,max(x)],[0.,5.])
    if gp.analytic:
        axs[2][0].plot(x, kappa_anf(x), 'blue')

    axs[2][0].set_ylabel('$\\kappa_{i,\\rm{LOS}}\\quad[1]$')
    axs[2][0].yaxis.set_major_locator(MaxNLocator(4))
    axs[2][0].xaxis.set_major_locator(MaxNLocator(4))

    
    if gp.pops==2:
        # plot nu 2
        axs[0][1].cla()
        x = gp.ipol.nux2_2D #[pc]
        if gp.geom == 'disc':
            x = gp.ipol.nux2
        lbound = gfun.compare_nu(2,True,False) - gfun.compare_nu(2,True,True) # [munit/pc**2]
        ubound = gfun.compare_nu(2,True,False) + gfun.compare_nu(2,True,True) # [munit/pc**2]
        axs[0][1].fill_between(x, lbound, ubound, alpha=0.5, color='green')
        axs[0][1].set_yscale('log')
        blim2 = [min(blim1[0],min(lbound)/1.5), max(blim1[1], max(ubound)*1.5)]
        setlims(axs[0][1],[0,max(x)],blim2)
        setlims(axs[0][0],[0,max(x)],blim2)
        
        # plot sigma_los_2
        axs[1][1].cla()
        x = gp.ipol.sigx2                          # [pc]
        lbound = (gp.ipol.sigdat2 - gp.ipol.sigerr2) # [km/s]
        ubound = (gp.ipol.sigdat2 + gp.ipol.sigerr2) # [km/s]
        axs[1][1].fill_between(x, lbound, ubound, alpha=0.5, color='green')
        blimsig2 = [0, max(blimsig1[1],max(ubound))]
        setlims(axs[1][1],[0,max(x)],blimsig2)

        # plot kurtosis_los_2
        axs[2][1].cla()
        x = gp.ipol.kapx2                            # [pc]
        lbound = (gp.ipol.kapdat2 - gp.ipol.kaperr2) # [1]
        ubound = (gp.ipol.kapdat2 + gp.ipol.kaperr2) # [1]
        axs[2][1].fill_between(x, lbound, ubound, alpha=0.5, color='green')

        setlims(axs[2][1],[0,max(x)],[0.,5.])

    if gp.plotdens:
        # plot 3D density
        axs[3][0].cla()
        x = gp.ipol.densx                          # [pc]
        if not gp.investigate=='fornax':
            lbound = (gp.ipol.densdat - gp.ipol.denserr) # [munit/pc**3]
            ubound = (gp.ipol.densdat + gp.ipol.denserr) # [munit/pc**3]
            axs[3][0].fill_between(x, lbound, ubound, alpha=0.5, color='green')

        if gp.geom == 'sphere':
            axs[3][0].set_ylabel('$\\rho\\quad[\\rm{M}_\\odot/\\rm{pc}^3]$')

        # show blue theoretical density line:
        if gp.model:
            x = gp.ipol.Mx_2D[:]                        # [pc]
            if gp.investigate == 'walker':
                rhodm, rhostar1, rhostar2 = rhowalker_3D(x) # 3*[munit/pc^3]
                rhotot = rhodm + rhostar1 + rhostar2
                axs[3][0].plot(x, (rhostar1+rhostar2), c='blue', lw=2, ls='-.')
            elif gp.investigate == 'gaia':
                rhotot = rhogaiatot_3D(x)
            elif gp.investigate == 'triaxial':
                rhotot = rhotriax(x)
            axs[3][0].plot(x, rhotot, c='blue', lw=2, ls='-.')
            setlims(axs[3][0],[0,max(x)],[min(rhotot)/5.,max(rhotot)*5.])

        if gp.analytic:
            x = gp.xipol[:]                        # [pc]
            rhotot = rho_anf(x)
            axs[3][0].plot(x, rhotot, c='blue', lw=2)
            setlims(axs[3][0],[0,max(x)],[min(rhotot)/5.,max(rhotot)*5.])

    else:
        # plot mass
        axs[3][0].cla()
        # show blue theoretical mass line
        if gp.model:
            x = gp.ipol.Mx_2D[:]                        # [pc]
            if gp.investigate == 'walker':
                Mtot=Mwalkertot(x)      # [Msun]
                axs[3][0].plot(x, Mtot, '-.', color='blue', lw=2)
                setlims(axs[3][0],[0,max(gp.xipol)],[0.5*min(Mtot),max(Mtot)*1.5])
                # axs[3][0].set_yscale('log') 

        x = gp.dat.Mx                          # [pc] TODO: gp.ipol.Mx_3D same here?
        lbound = (gp.dat.Mdat - gp.dat.Merr) # [munit,3D]
        ubound = (gp.dat.Mdat + gp.dat.Merr) # [munit,3D]

        # x = gp.ipol.Mx_2D                            # [pc]
        # lbound = (gp.ipol.Mdat_2D - gp.ipol.Merr_2D) # [munit,2D]
        # ubound = (gp.ipol.Mdat_2D + gp.ipol.Merr_2D) # [munit,2D]
        axs[3][0].fill_between(x, lbound, ubound, alpha=0.5, color='green')
        #axs[3][0].set_ylim([0.,1.2*max(ubound)])
        # axs[3][0].plot(x, M_anf(gp.ipol.Mx)*gp.totmass[0], c='blue', lw=1) # only for Hernquist

        if gp.geom == 'sphere':
            axs[3][0].set_ylabel('$M\\quad[\\rm{M}_\\odot]$')
        elif gp.geom == 'disc':
            axs[3][0].set_ylabel('$\\Sigma\\quad[\\rm{M}_\\odot/\\rm{pc}^2]$')


    axs[3][0].set_xlabel('$r\\quad[\\rm{pc}]$')
    axs[3][1].set_xlabel('$r\\quad[\\rm{pc}]$')

    if gp.log:   axs[3][0].set_yscale('log')
    plt.draw()
    return
    
## get data for 1 population
def get_plot_data_1():
    x = gp.xipol[:]
    r_tot = gp.xipol[:]                   # TODO: extension to 2* rmax
    grav = gp.dens_x if gp.plotdens else gp.M_x

    x1 = gp.ipol.nux1[:]
    nu1 = rho_INT_Rho(x1, gp.nu1_x)
    if gp.geom == 'disc': nu1 = gp.nu1_x
    sig1 = gp.sig1_x
    kap1 = gp.kap1_x
    return x1, r_tot, nu1, gp.sig1_x, gp.kap1_x, grav, gp.d1_x

## get data for 2 populations
def get_plot_data_2():
    x = gp.xipol[:]
    r_tot = gp.xipol[:]                   # TODO: extension to 2* rmax

    grav = gp.dens_x if gp.plotdens else gp.M_x
    x1 = gp.ipol.nux1[:]; x2 = gp.ipol.nux2[:]
    nu1 = rho_INT_Rho(x1, gp.nu1_x)
    nu2 = rho_INT_Rho(x2, gp.nu2_x)
    if gp.geom == 'disc':
        nu1 = gp.nu1_x
        nu2 = gp.nu2_x
        
    return x1, r_tot, nu1, gp.sig1_x, gp.kap1_x, grav, gp.d1_x,\
           nu2, gp.sig2_x, gp.kap2_x, gp.d2_x


## plot the first red line (values from mcmc_init)
def plot_first_guess():
    global f1,f2,f3,f4,f5,f6,f7,f8
    if not gp.showplot: return
    if gp.pops==1:
        xpl, x_tot, nu1, sig1, kap1, grav, delta1 = get_plot_data_1()
    elif gp.pops==2:
        xpl, x_tot, nu1, sig1, kap1, grav, delta1, nu2, sig2, kap2, delta2 = get_plot_data_2()
    f1, = axs[0][0].plot(xpl, nu1,  c='red')
    f2, = axs[1][0].plot(xpl, sig1, c='red')
    if gp.geom == 'sphere':
        f3, = axs[2][0].plot(xpl, kap1, c='red')
        axs[2][0].set_ylim([0.,5.])
        
    if (gp.pops == 2) : 
        f4, = axs[0][1].plot(xpl, nu2,  c='red')
        f5, = axs[1][1].plot(xpl, sig2, c='red')
        if gp.geom == 'sphere':
            f6, = axs[2][1].plot(xpl, kap2, c='red')
            axs[2][1].set_ylim([0.,5.])

    f7, = axs[3][0].plot(x_tot, grav,  c='red')
    f8, = axs[3][1].plot(xpl,   delta1, c='red')
    # if gp.geom == 'sphere':
    #     axs[3][0].set_ylim([5e-3,3.])
    if gp.geom == 'disc':
        f7, = axs[3][0].plot(x_tot, gp.Mmodel, c='blue',ls='dashed')

    axs[3][1].set_ylim([-0.5,1.])
    # axs[3][1].xaxis.set_major_locator(MaxNLocator(4))
    if (gp.pops == 2) : axs[3][1].plot(xpl,delta2,color='orange')
    plt.draw()
    return


## update plot to new model values
def update_plot():
    if not gp.showplot:
        return
    if gp.pops==1:
        xpl,x_tot,nu1,sig1,kap1,grav,d1 = get_plot_data_1()
    elif gp.pops==2:
        xpl,x_tot,nu1,sig1,kap1,grav,d1,nu2,sig2,kap2,d2 = get_plot_data_2()
    f1.set_ydata(nu1)
    f2.set_ydata(sig1)
    if gp.geom == 'sphere': f3.set_ydata(kap1)
    if gp.pops==2:
        f4.set_ydata(nu2)
        f5.set_ydata(sig2)
        if gp.geom == 'sphere': f6.set_ydata(kap2)
    f7.set_ydata(grav)

    # plot delta
    axs[3][1].cla()
    axs[3][1].plot(xpl,d1,color='r')
    axs[3][1].yaxis.set_major_locator(MaxNLocator(4))
    axs[3][1].yaxis.tick_right()
    axs[3][1].yaxis.set_label_position("right")
    axs[3][1].set_xlabel('$r\\quad[\\rm{pc}]$')
    axs[3][1].set_ylabel('$\\delta_{1,2}$')

    # plot disc stuff
    if gp.geom == 'disc':
        # plot kappa
        axs[2][0].cla()
        axs[2][0].plot(xpl, gp.parst.dens,'r')    # overall kappa
        axs[2][0].plot(xpl, gp.parst.dens - phys.kappa(gp.xipol, -gp.blow*2.*np.pi*gp.G1), 'k')
        # scale to available range
        axs[2][0].set_ylim([-1.e-3,7.e-3])
        axs[2][0].yaxis.set_major_locator(MaxNLocator(4))
        axs[2][0].set_xlabel('$r\\quad[\\rm{pc}]$')
        axs[2][0].set_ylabel('$\\kappa_{tot,DM}$')

    if gp.pops==2: axs[3][1].plot(xpl,d2,color='orange')
    if gp.geom == 'sphere': axs[3][1].set_ylim([-0.5,1.0])
    axs[0][0].set_xlim([np.min(gp.xipol),np.max(gp.xipol)])
    axs[0][0].xaxis.set_major_locator(MaxNLocator(4))
    if gp.lograd: axs[0][0].set_xscale('log')
    plt.draw()
    # plt.savefig('first.png')
    if not gp.initphase: save_plot()
    return
    
## save plot to png file
def save_plot():
    plt.savefig(gp.files.get_outpng())
    return

## show all plots, stop execution of program. Have to close via (x) in window
def show_plots():
    if not gp.showplot: 
        return
    ioff();show()
    return


## plot profiles. TODO: seems not to be used anywhere.. delete
def plot_disc():
    if not showplot: 
        return

    # Calculate profiles:
    nu_z  = phys.nu(gp.xipol,gp.xipol,gp.parst.nu1,gp.quadratic,gp.monotonic)
    Rsun = 8.; hr = 3.0; hsig = 3.0
    sig_z = phys.sigma_z(gp.xipol,gp.xipol,gp.parst.dens,gp.blow,gp.parst.nu1,gp.parst.delta1,[Rsun,hr,hsig],gp.quadratic,gp.monotonic)
    kz_z  = phys.kz(gp.xipol,gp.xipol,gp.parst.dens,gp.blow,gp.quadratic)

    # Convert kz_z to Sigma in Msun / pc^2:
    Msigma_zpl = abs(kz_z) / (2.*np.pi*G1) / 1000.^2.
    
    # Do the plots: 
    plot(zpl,nu_zpl/max(nu_zpl),yrange=[numin,numax],xrange=[0.,zplmax],title='!6',xtitle='z(kpc)',ytitle='nu [units]')
    # TODO: check which ax, assign yscale('log')

    sel = (z_dat_bin > 0)
    zdp = z_dat_bin[sel];    ndp = nu_dat_bin[sel]
    # if gp.initphase:
    #     ndp = nu_dat_bin[sel]
    # else:
    #     ndp = nu_dat1
    ndpe = nu_dat_err_bin[sel]
    mnorm = 1.*max(ndp)
    ndp  = ndp/mnorm;     ndpe = ndpe/mnorm
        
    errorbar(zdp,ndp,ndpe)
    plot(zpl,nu_zpl/max(nu_zpl))
    if investigate == 'simple': plot(zth,nu_zth)
    
    plot(zpl,sig_zpl,xrange=[0,zplmax],yrange=[0,55],\
         title='!6',xtitle='z(kpc)',ytitle='sigma_z [km/s]')
    
    sdp = sig_dat_bin[sel]
    # if gp.initphase: sdp = sig_dat_bin[sel]
    #     else: sdp = sig_dat1
    sdpe = sig_dat_err_bin[sel]
    
    errorbar(zdp,sdp,sdpe)
    plot(zpl,sig_zpl)
    if investigate == 'simple': plot(zstar,sigzstar)
    
    plot(zpl,Msigma_zpl,title='!6',xtitle='z(kpc)',ytitle='Sigma_z [Msun/pc^2]',\
         yrange=[0,80],xrange=[0,zplmax])
    if investigate != 'simple':
        errorbar(zvis,sigusevis,siguseviserr)
        errorbar(zsurf,sigdm_surf,sigdm_surferr)
        errorbar(zsurf,sigdm_surf/2.,sigdm_surferr)
        plot(zpl,Msigma_zpl)
        plot(zsurf,sigs_surf+sigdm_surf)
        plot(zsurf,sigs_surf+sigdm_surf/2.)
    else:
        # Overplot theory curve: 
        plot(zth,abs(Kz_zth) / (2.*np.pi*G1) / 1000^2.)
        if bprior: plot(zvis,sigusevis)
        
    # plot tilt
    if tparsRin(0) > 0:  
        sig_Rz = sigma_rz(zpl,zp_kz,tpars)
        plot(zpl,sig_Rz,xtitle='z(kpc)',ytitle='Sig_Rz',yrange=[-200,200])
        errorbar(z_dat_bin,sigRz_dat_bin,sigRz_dat_err_bin)
    else:
        if gp.vcorrect :  
            vz_mean = vzmean(zpl,zp_kz,vmpars)
            plot(zpl,vz_mean,xtitle='z(kpc)',ytitle='<v_z>(km/s)',yrange=[-10,10])
        else:
            # Dark matter density here: 
            kzparsu = abs(kzpars) 
            denarr = np.zeros(len(kzparsu))
            denarr[0] = kzparsu[0]
            for i in range(1,len(kzparsu)):
                denarr[i] = denarr[i-1] + kzparsu[i]
            denarr = reverse(denarr)
            
            plot(zp_kz,denarr / (4.0*np.pi*G1) / 1000^3.,yrange = [0,0.03],\
                 xtitle='z(kpc)', ytitle='den [Msun/pc^3]')

            if investigate == 'simple':  
                # Overplot theory curve: 
                Kz_zthd = -2.0 * F * zth
                if adddarkdisc :
                    Kz_zthd = Kz_zthd - KD*zth/sqrt(zth^2.+DD^2.)
                Sigz_zth = abs(Kz_zthd) / (2.0*np.pi*G1) / 1000^2. 
                denth = deriv(zth,Sigz_zth) / 1000.
                plot(zth,denth)
    return
