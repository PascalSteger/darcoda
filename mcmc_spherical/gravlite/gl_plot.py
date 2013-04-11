#!/usr/bin/python2.7
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
'''all functions to work with plots'''
import gl_params as gp
import gl_funs as dfun
from gl_funs import compare_nu
from gl_analytic import *
from gl_int import *
if gp.geom == 'sphere':
    import physics_sphere as phys
elif gp.geom == 'disc':
    import physics_disc as phys
from pylab import *
import numpy as np
import gl_units as units
import pdb

global f, ax1, ax2, ax3, ax4, ax5, ax6





def prepare_plots():
    global f,ax1,ax2,ax3,ax4,ax5,ax6
    if not gp.testplot: return
    ion()
    if(gp.pops==1):
        f, ((ax1,ax2),(ax5,ax6)) = plt.subplots(2, 2)
    if(gp.pops==2):
        f, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3, 2)
    draw()
    return





def setlabels(ax,xtext,ytext):
    ax.set_xlabel(r'$'+xtext+'$')
    ax.set_ylabel(r'$'+ytext+'$')

def setlims(ax,xlim,ylim):
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)






def plot_data():

    if not gp.testplot: return

    ##### plot nu 1
    if gp.analytic: ax1.plot(gp.dat.nux1, rho_anf(gp.dat.nux1), c='blue', lw=5)
    x = gp.dat.nux1 # [pc]
    lbound = compare_nu(1,True,False) - compare_nu(1,True,True) #[munit/pc**2]
    ubound = compare_nu(1,True,False) + compare_nu(1,True,True) #[munit/pc**2]
    ax1.fill_between(x,lbound,ubound,alpha=0.5,color='g') #[pc], 2*[munit/pc**2]
    ax1.set_yscale('log')
    setlims(ax1,[0,max(x)],[min(lbound),max(ubound)])
    setlabels(ax1,'r [pc]','\nu [Msun/pc^2]')

    ##### sigma_LOS 1
    if gp.analytic: ax2.plot(gp.dat.sigx1, sig_los_anf(gp.dat.sigx1), c='blue', lw=5)
    x = gp.dat.sigx1 #[pc]
    lbound = (gp.dat.sigdat1 - gp.dat.sigerr1) # [km/s]
    ubound = (gp.dat.sigdat1 + gp.dat.sigerr1) # [km/s]
    ax2.fill_between(x,lbound,ubound,alpha=0.5,color='g') #[pc], 2*[km/s]
    setlims(ax2,[0,max(x)],[min(lbound),max(ubound)])
    ## TODO: setlabels(ax2,'r [\text{pc}]','\sigma_{{1,LOS}} [\text{km}/\text{s}]')
    
    if gp.pops==2:
        #####  nu 2
        x = gp.dat.nux2_2D #[pc]
        lbound = compare_nu(2,True,False) - compare_nu(2,True,True)# [munit/pc**2]
        ubound = compare_nu(2,True,False) + compare_nu(2,True,True) # [munit/pc**2]
        ax3.fill_between(x, lbound, ubound, alpha=0.5, color='g')
        ax3.set_yscale('log')
        setlims(ax3,[0,max(x)],[min(lbound),max(ubound)])
        ## TODO: setlabels(ax3,'r [\text{pc}]','\nu_2 [\text{M}_\odot/\text{pc}^2]')
        
        ##### sigmalos 2
        x = gp.dat.sigx2                                      # [pc]
        lbound = (gp.dat.sigdat2 - gp.dat.sigerr2)            # [km/s]
        ubound = (gp.dat.sigdat2 + gp.dat.sigerr2)            # [km/s]
        ax4.fill_between(x, lbound, ubound, alpha=0.5, color='g')
        setlims(ax4,[0,max(x)],[min(lbound),max(ubound)])
        ## TODO: setlabels(ax4,'r [pc]','\sigma_{los,2} [km/s]')

    if gp.plotdens:
        # #### density 2D
        # x = gp.dat.densx_2D * gp.rcore[0]        # [pc]
        # lbound = (gp.dat.densdat_2D - gp.dat.denserr_2D)*gp.dens0pc_2D[0] #[munit/pc^2]
        # ubound = (gp.dat.densdat_2D + gp.dat.denserr_2D)*gp.dens0pc_2D[0] #[munit/pc^2]
        # ax5.fill_between(x, lbound, ubound, alpha=0.5, color='g')


        # ## TODO: setlabels(ax5,'r [pc]','dens [Msun/pc^2]')
        # # if gp.lim:   setlims(ax5,[0,max(x)/2],[1e-4,1.5*max(den)])
        # if gp.model:
        #     x = gp.dat.Mx_2D * gp.rcore[0]                      # [pc]
        #     rhodm, rhostar1, rhostar2 = rhowalker_2D(x) # 3*[munit/pc^2]
        #     ax5.plot(x, (rhodm+rhostar1+rhostar2), c='blue', lw=4)
        #     ax5.plot(x, (rhostar1+rhostar2), c='blue', lw=3, ls='-.')


        #### density 3D
        x = gp.dat.densx        # [pc]
        lbound = (gp.dat.densdat - gp.dat.denserr) #[munit/pc**3]
        ubound = (gp.dat.densdat + gp.dat.denserr) #[munit/pc**3]
        ax5.fill_between(x, lbound, ubound, alpha=0.5, color='g')

        ## TODO: setlabels(ax5,'r [pc]','dens [Msun/pc^2]')
        # if gp.lim: setlims(ax5,[0,max(x)],[min(lbound),max(ubound)])

        if gp.model:
            x = gp.dat.Mx_2D[:]                         # [pc]
            rhodm, rhostar1, rhostar2 = rhowalker_3D(x) # 3*[munit/pc^3]
            rhotot = rhodm + rhostar1 + rhostar2
            ax5.plot(x, rhotot, c='blue', lw=4)
            setlims(ax5,[0,max(x)],[min(rhotot),max(rhotot)])
            ax5.plot(x, (rhostar1+rhostar2), c='blue', lw=3, ls='-.')


    else:
        #### mass
        x = gp.dat.Mx_2D #[pc]
        lbound = (gp.dat.Mdat_2D - gp.dat.Merr_2D) # [munit,2D]
        ubound = (gp.dat.Mdat_2D + gp.dat.Merr_2D) # [munit,2D]
        ax5.fill_between(x, lbound, ubound, alpha=0.5, color='g')
        #ax5.plot(x, M_anf(gp.dat.Mx)*gp.totmass[0], c='blue', lw=4) # only for Hernquist
        ## TODO: setlabels(ax5,'r [pc]','M [Msun]')
        if gp.lim: ax5.set_ylim([0.,1.2*max(ubound)])
        if gp.lim:   setlims(ax5,[0,max(x)],[1e-4,1.5*max(ubound)])

    if gp.log:   ax5.set_yscale('log')
    plt.draw()





    
def get_plot_data():
    '''determine physical quantities from densities and sigma_LOS'''
    x = gp.xipol[:]
    r_tot = gp.xipol[:]                   # TODO: extension to 2* rmax

    if gp.plotdens:
        M_tot = gp.dens_x
    else:
        M_tot = int_project(x, gp.dens_x)

    if gp.pops==1:
        x1 = gp.ipol.nux1[:]
        nu1 = int_surfden(gp.nipol, x1, gp.nu1_x)
        sig1 = gp.sig1_x
        return x1, r_tot, nu1, sig1, M_tot, gp.parst.delta1
    elif gp.pops==2:
        x1 = gp.ipol.nux1[:]; x2 = gp.ipol.nux2[:]

        nu1 = int_surfden(gp.nipol, x1, gp.nu1_x)
        nu2 = int_surfden(gp.nipol, x2, gp.nu2_x)

        sig1 = gp.sig1_x; sig2 = gp.sig2_x
        return x1, r_tot, nu1, sig1, M_tot, gp.parst.delta1, nu2, sig2, gp.parst.delta2





def plot_first_guess():
    global f1,f2,f3,f4,f5,f6
    if not gp.testplot: return
    if gp.pops==1:
        xpl, x_tot, nu1, sig1, M_tot, delta1 = get_plot_data()
    elif gp.pops==2:
        xpl, x_tot, nu1, sig1, M_tot, delta1, nu2, sig2, delta2 = get_plot_data()
    f1, = ax1.plot(xpl, nu1,  c='red')
    f2, = ax2.plot(xpl, sig1, c='red')

    if (gp.pops == 2) : 
        f3, = ax3.plot(xpl, nu2,  c='red')
        f4, = ax4.plot(xpl, sig2, c='red')

    f5, = ax5.plot(x_tot, M_tot,  c='red')
    f6, = ax6.plot(xpl,   delta1, c='red')
    ax6.set_ylim([-0.5,1.])
    if (gp.pops == 2) : ax6.plot(xpl,delta2,color='orange')
    plt.draw()
    save_plot()
    return




    
def update_plot():
    if not gp.testplot: return
    if gp.pops==1:
        xpl,x_tot,nu1,sig1,M_tot,delta1 = get_plot_data()
    elif gp.pops==2:
        xpl,x_tot,nu1,sig1,M_tot,delta1,nu2,sig2,delta2 = get_plot_data()
    f1.set_ydata(nu1)
    f2.set_ydata(sig1)
    if gp.pops==2:
        f3.set_ydata(nu2)
        f4.set_ydata(sig2)
    f5.set_ydata(M_tot)
    ax6.cla()
    ax6.plot(xpl,delta1,color='r')
    setlabels(ax6,'r(R_s)','delta_{{1,2}}')    
    if gp.pops==2: ax6.plot(xpl,delta2,color='orange')
    ax6.set_ylim([-0.5,1.0])
    plt.draw()
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
    # if initphase == 'start':
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
    # if initphase == 'start': sdp = sig_dat_bin[sel]
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
