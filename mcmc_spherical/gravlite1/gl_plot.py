#!/usr/bin/python
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
'''all functions to work with plots'''
import gl_params as gp
import gl_funs as gfun
from gl_analytic import *
from gl_int import *
if gp.geom == 'sphere':
    import physics_sphere as phys
elif gp.geom == 'disc':
    import physics_disc as phys
from pylab import *
import numpy as np
import gl_units as units
from matplotlib.ticker import MaxNLocator

import pdb

global f, axs





def prepare_plots():
    global f,axs
    if not gp.testplot: return
    ion()
    # f, axs = plt.subplots(3, 2, sharex=True, sharey=True)
    f = figure()
    ax1 = f.add_subplot(421)
    ax2 = f.add_subplot(422, sharex=ax1, sharey=ax1)
    ax3 = f.add_subplot(423, sharex=ax1)
    ax4 = f.add_subplot(424, sharex=ax1, sharey=ax3)
    ax5 = f.add_subplot(425, sharex=ax1)
    ax6 = f.add_subplot(426, sharex=ax1)
    ax7 = f.add_subplot(427, sharex=ax1)
    ax8 = f.add_subplot(428, sharex=ax1, sharey=ax7)
    axs=[[ax1, ax2],[ax3,ax4],[ax5,ax6],[ax7,ax8]]
    xticklabels = ax1.get_xticklabels()+ax2.get_xticklabels()
    xticklabels += ax3.get_xticklabels()+ax4.get_xticklabels()
    xticklabels += ax5.get_xticklabels()+ax6.get_xticklabels()
    setp(xticklabels, visible=False)
    setp(ax2.get_yticklabels(), visible=False)
    setp(ax4.get_yticklabels(), visible=False)
    setp(ax6.get_yticklabels(), visible=False)
    setp(ax8.get_yticklabels(), visible=False)
    draw()
    return





def setlims(ax,xlim,ylim):
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)






def plot_data():

    if not gp.testplot: return
    ##### plot nu 1
    axs[0][0].cla()
    x = gp.ipol.nux1_2D                 # [pc]
    if gp.geom == 'disc': x = gp.ipol.nux1
    if gp.analytic: axs[0][0].plot(x, surfden_anf(x), c='blue', lw=4)
    lbound = gfun.compare_nu(1,True,False) - gfun.compare_nu(1,True,True) # [munit/pc**2]
    # gives gp.ipol.nudat1_2D
    ubound = gfun.compare_nu(1,True,False) + gfun.compare_nu(1,True,True) # [munit/pc**2]
    axs[0][0].fill_between(x,lbound,ubound,alpha=0.5,color='g') #[pc], 2*[munit/pc**2]
    if gp.lograd: axs[0][0].set_xscale('log')
    axs[0][0].set_yscale('log')
    blim1 = [min(lbound)/1.5, 1.5*max(ubound)]
    setlims(axs[0][0],[0,max(x)],blim1)
    setlims(axs[0][1],[0,max(x)],blim1)
    axs[0][0].set_ylabel('$\\nu_i\\quad[\\rm{M}_\odot/\\rm{pc}^2]$')
    axs[0][0].xaxis.set_major_locator(MaxNLocator(4))

    ##### sigma_LOS 1
    axs[1][0].cla()
    if gp.analytic: axs[1][0].plot(gp.ipol.sigx1, sig_los_anf(gp.ipol.sigx1), c='blue', lw=5)
    x = gp.ipol.sigx1                                         # [pc]
    lbound = (gp.ipol.sigdat1 - gp.ipol.sigerr1)                # [km/s]
    ubound = (gp.ipol.sigdat1 + gp.ipol.sigerr1)                # [km/s]
    axs[1][0].fill_between(x,lbound,ubound,alpha=0.5,color='green') # [pc], 2*[km/s]
    blimsig1 = [0,1.5*max(ubound)]
    setlims(axs[1][0],[0,max(x)],blimsig1)
    setlims(axs[1][1],[0,max(x)],blimsig1)
    axs[1][0].set_ylabel('$\\sigma_{i,\\rm{LOS}}\\quad[\\rm{km/s}]$')
    axs[1][0].yaxis.set_major_locator(MaxNLocator(4))
    axs[1][0].xaxis.set_major_locator(MaxNLocator(4))
    
    if gp.pops==2:
        #####  nu 2
        axs[0][1].cla()
        x = gp.ipol.nux2_2D #[pc]
        if gp.geom == 'disc': x = gp.ipol.nux2
        lbound = gfun.compare_nu(2,True,False) - gfun.compare_nu(2,True,True) # [munit/pc**2]
        ubound = gfun.compare_nu(2,True,False) + gfun.compare_nu(2,True,True) # [munit/pc**2]
        axs[0][1].fill_between(x, lbound, ubound, alpha=0.5, color='green')
        axs[0][1].set_yscale('log')
        blim2 = [min(blim1[0],min(lbound)/1.5), max(blim1[1], max(ubound)*1.5)]
        setlims(axs[0][1],[0,max(x)],blim2)
        setlims(axs[0][0],[0,max(x)],blim2)
        
        ##### sigmalos 2
        axs[1][1].cla()
        x = gp.ipol.sigx2                          # [pc]
        lbound = (gp.ipol.sigdat2 - gp.ipol.sigerr2) # [km/s]
        ubound = (gp.ipol.sigdat2 + gp.ipol.sigerr2) # [km/s]
        axs[1][1].fill_between(x, lbound, ubound, alpha=0.5, color='green')
        blimsig2 = [0, max(blimsig1[1],max(ubound))]
        setlims(axs[1][1],[0,max(x)],blimsig2)

    if gp.plotdens:
        # #### density 2D
        # x = gp.ipol.densx_2D * gp.rcore[0]        # [pc]
        # lbound = (gp.ipol.densdat_2D - gp.ipol.denserr_2D)*gp.dens0pc_2D[0] #[munit/pc^2]
        # ubound = (gp.ipol.densdat_2D + gp.ipol.denserr_2D)*gp.dens0pc_2D[0] #[munit/pc^2]
        # axs[2][0].fill_between(x, lbound, ubound, alpha=0.5, color='g')


        # # if gp.lim:   setlims(axs[2][0],[0,max(x)/2],[1e-4,1.5*max(den)])
        # if gp.model:
        #     x = gp.ipol.Mx_2D * gp.rcore[0]                      # [pc]
        #     rhodm, rhostar1, rhostar2 = rhowalker_2D(x) # 3*[munit/pc^2]
        #     axs[2][0].plot(x, (rhodm+rhostar1+rhostar2), c='blue', lw=4)
        #     axs[2][0].plot(x, (rhostar1+rhostar2), c='blue', lw=3, ls='-.')


        #### density 3D
        axs[2][0].cla()
        x = gp.ipol.densx                          # [pc]
        lbound = (gp.ipol.densdat - gp.ipol.denserr) # [munit/pc**3]
        ubound = (gp.ipol.densdat + gp.ipol.denserr) # [munit/pc**3]
        axs[2][0].fill_between(x, lbound, ubound, alpha=0.5, color='green')

        if gp.geom == 'sphere':
            axs[2][0].set_ylabel('$\\rho\\quad[\\rm{M}_\\odot/\\rm{pc}^3]$')
            

        # if gp.lim: setlims(axs[2][0],[0,max(x)],[min(lbound),max(ubound)])

        if gp.model:
            x = gp.ipol.Mx_2D[:]                        # [pc]
            rhodm, rhostar1, rhostar2 = rhowalker_3D(x) # 3*[munit/pc^3]
            rhotot = rhodm + rhostar1 + rhostar2
            axs[2][0].plot(x, rhotot, c='blue', lw=4)
            setlims(axs[2][0],[0,max(x)],[min(rhotot),max(rhotot)])
            axs[2][0].plot(x, (rhostar1+rhostar2), c='blue', lw=3, ls='-.')

        if gp.analytic:
            x = gp.xipol[:]                        # [pc]
            rhotot = rho_anf(x)
            axs[2][0].plot(x, rhotot, c='blue', lw=4)
            setlims(axs[2][0],[0,max(x)],[min(rhotot),max(rhotot)])

    else:
        #### mass
        axs[2][0].cla()
        x = gp.dat.Mx                          # [pc] TODO: gp.ipol.Mx_2D same here?
        lbound = (gp.dat.Mdat - gp.dat.Merr) # [munit,2D]
        ubound = (gp.dat.Mdat + gp.dat.Merr) # [munit,2D]
        #### mass
        #x = gp.ipol.Mx_2D                            # [pc]
        #lbound = (gp.ipol.Mdat_2D - gp.ipol.Merr_2D) # [munit,2D]
        #ubound = (gp.ipol.Mdat_2D + gp.ipol.Merr_2D) # [munit,2D]
        axs[2][0].fill_between(x, lbound, ubound, alpha=0.5, color='green')
        #axs[2][0].plot(x, M_anf(gp.ipol.Mx)*gp.totmass[0], c='blue', lw=4) # only for Hernquist
        if gp.geom == 'sphere':
            axs[2][0].set_ylabel('$M\\quad[\\rm{M}_\\odot]$')
        elif gp.geom == 'disc':
            axs[2][0].set_ylabel('$\\Sigma\\quad[\\rm{M}_\\odot/\\rm{pc}^2]$')

        if gp.lim:
            axs[2][0].set_ylim([0.,1.2*max(ubound)])
            #setlims(axs[2][0],[0,max(x)],[1e-4,1.5*max(ubound)])

    axs[2][0].set_xlabel('$r\\quad[\\rm{pc}]$')
    axs[2][1].set_xlabel('$r\\quad[\\rm{pc}]$')

    if gp.log:   axs[2][0].set_yscale('log')
    plt.draw()





    
def get_plot_data():
    '''determine physical quantities from densities and sigma_LOS'''
    x = gp.xipol[:]
    r_tot = gp.xipol[:]                   # TODO: extension to 2* rmax

    if gp.plotdens:
        M_tot = gp.dens_x
    else:
        if gp.geom == 'sphere':
            M_tot = int_project(x, gp.dens_x)
        else:
            M_tot = gp.M_x              # [Msun/pc**2] <= /kpc**2? TODO

    if gp.pops==1:
        x1 = gp.ipol.nux1[:]
        nu1 = int_surfden(x1, gp.nu1_x)
        if gp.geom == 'disc': nu1 = gp.nu1_x
        sig1 = gp.sig1_x
        return x1, r_tot, nu1, sig1, M_tot, gp.d1_x
    elif gp.pops==2:
        x1 = gp.ipol.nux1[:]; x2 = gp.ipol.nux2[:]

        nu1 = int_surfden(x1, gp.nu1_x)
        nu2 = int_surfden(x2, gp.nu2_x)
        if gp.geom == 'disc':
            nu1 = gp.nu1_x
            nu2 = gp.nu2_x
        
        sig1 = gp.sig1_x; sig2 = gp.sig2_x
        return x1, r_tot, nu1, sig1, M_tot, gp.d1_x, nu2, sig2, gp.d2_x





def plot_first_guess():
    global f1,f2,f3,f4,f5,f6
    if not gp.testplot: return
    if gp.pops==1:
        xpl, x_tot, nu1, sig1, M_tot, delta1 = get_plot_data()
    elif gp.pops==2:
        xpl, x_tot, nu1, sig1, M_tot, delta1, nu2, sig2, delta2 = get_plot_data()
    f1, = axs[0][0].plot(xpl, nu1,  c='red')
    f2, = axs[1][0].plot(xpl, sig1, c='red')

    if (gp.pops == 2) : 
        f3, = axs[0][1].plot(xpl, nu2,  c='red')
        f4, = axs[1][1].plot(xpl, sig2, c='red')

    f5, = axs[2][0].plot(x_tot, M_tot,  c='red')
    f6, = axs[2][1].plot(xpl,   delta1, c='red')
    if gp.geom == 'sphere':
        axs[2][0].set_ylim([1e-2,1.])
    if gp.geom == 'disc':
        f7, = axs[2][0].plot(x_tot, gp.Mmodel, c='blue',ls='dashed')
    # axs[2][1].set_ylim([-0.5,1.])
    # axs[2][1].xaxis.set_major_locator(MaxNLocator(4))
    if (gp.pops == 2) : axs[2][1].plot(xpl,delta2,color='orange')
    plt.draw()
    save_plot()
    return




    
def update_plot():
    if not gp.testplot: return
    if gp.pops==1:
        xpl,x_tot,nu1,sig1,M_tot,d1 = get_plot_data()
    elif gp.pops==2:
        xpl,x_tot,nu1,sig1,M_tot,d1,nu2,sig2,d2 = get_plot_data()
    f1.set_ydata(nu1)
    f2.set_ydata(sig1)
    if gp.pops==2:
        f3.set_ydata(nu2)
        f4.set_ydata(sig2)
    f5.set_ydata(M_tot)

    # plot delta
    axs[2][1].cla()
    axs[2][1].plot(xpl,d1,color='r')
    axs[2][1].yaxis.set_major_locator(MaxNLocator(4))
    axs[2][1].yaxis.tick_right()
    axs[2][1].yaxis.set_label_position("right")
    axs[2][1].set_xlabel('$r\\quad[\\rm{pc}]$')
    axs[2][1].set_ylabel('$\\delta_{1,2}$')
    
    # plot kappa
    axs[3][0].cla()
    axs[3][0].plot(xpl, gp.parst.dens,'r')    # overall kappa
    axs[3][0].plot(xpl, gp.parst.dens - phys.kappa(gp.xipol, -gp.blow*2.*np.pi*gp.G1), 'k')
    # scale to available range
    axs[3][0].set_ylim([-1.e-3,7.e-3])
    axs[3][0].yaxis.set_major_locator(MaxNLocator(4))
    axs[3][0].set_xlabel('$r\\quad[\\rm{pc}]$')
    axs[3][0].set_ylabel('$\\kappa_{tot,DM}$')

    if gp.pops==2: axs[2][1].plot(xpl,d2,color='orange')
    if gp.geom == 'sphere': axs[2][1].set_ylim([-0.5,1.0])
    axs[0][0].set_xlim([np.min(gp.xipol),np.max(gp.xipol)])
    axs[0][0].xaxis.set_major_locator(MaxNLocator(4))
    if gp.lograd: axs[0][0].set_xscale('log')
    plt.draw()
    plt.savefig('first.png')
    save_plot()
    return





    
def save_plot():
    plt.savefig(gp.files.get_outpng())
    return




def show_plots():
    if not gp.testplot: return
    ioff();show()
    return



    

def plot_disc():
    # plot profiles
    if not testplot: return

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
        
    # Tilt:
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
