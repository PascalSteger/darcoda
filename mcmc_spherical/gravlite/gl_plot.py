#!/usr/bin/python
'''all functions to work with plots'''
import gl_params as gp
import gl_funs as dfun
from gl_analytic import *
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
    x = gp.dat.nux1*gp.rcore[0] #[pc]
    lbound = (gp.dat.nudat1 - gp.dat.nuerr1)*gp.dens0pc[0] #[munit/pc**2]
    ubound = (gp.dat.nudat1 + gp.dat.nuerr1)*gp.dens0pc[0] #[munit/pc**2]
    ax1.fill_between(x,lbound,ubound,alpha=0.5,color='g') #[rcore], [dens0]
    ax1.set_yscale('log')
    setlabels(ax1,'r [pc]','\nu [Msun/pc^2]')

    ##### sigma_LOS 1
    if gp.analytic: ax2.plot(gp.dat.sigx1, sig_los_anf(gp.dat.sigx1), c='blue', lw=5)
    # ax2.plot(gp.dat.sigx1, gp.dat.sigdat1, c='green',lw=3)
    # ax2.errorbar(gp.dat.sigx1, gp.dat.sigdat1, yerr=gp.dat.sigerr1, c='green')
    x = gp.dat.sigx1 * gp.rcore[0] #[pc]
    lbound = (gp.dat.sigdat1 - gp.dat.sigerr1)*gp.maxvlos[0] #[km/s]
    ubound = (gp.dat.sigdat1 + gp.dat.sigerr1)*gp.maxvlos[0] #[km/s]
    ax2.fill_between(x,lbound,ubound,alpha=0.5,color='g')
    #setlims(ax2,[0,max(x)],[0,1.2*max(lbound)])
    ## TODO: setlabels(ax2,'r [\text{pc}]','\sigma_{{1,LOS}} [\text{km}/\text{s}]')
    
    if gp.pops==2:
        #####  nu 2
        # ax3.plot(gp.dat.nux2, gp.dat.nudat2, c='green', lw=3)
        # ax3.errorbar(gp.dat.nux2,gp.dat.nudat2,yerr=gp.dat.nuerr2, c='green')
        x = gp.dat.nux2*gp.rcore[1] #[pc]
        lbound = (gp.dat.nudat2 - gp.dat.nuerr2)*gp.dens0pc[1]
        ubound = (gp.dat.nudat2 + gp.dat.nuerr2)*gp.dens0pc[1]
        ax3.fill_between(x, lbound, ubound, alpha=0.5, color='g')
        ax3.set_yscale('log')
        ## TODO: setlabels(ax3,'r [\text{pc}]','\nu_2 [\text{M}_\odot/\text{pc}^2]')
        
        ##### sigmalos 2
        # ax4.plot(gp.dat.sigx2,     gp.dat.sigdat2, c='green', lw=3)
        # ax4.errorbar(gp.dat.sigx2, gp.dat.sigdat2, yerr=gp.dat.sigerr2, c='green')
        x = gp.dat.sigx2 * gp.rcore[1]
        lbound = (gp.dat.sigdat2 - gp.dat.sigerr2)*gp.maxvlos[1]
        ubound = (gp.dat.sigdat2 + gp.dat.sigerr2)*gp.maxvlos[1]
        ax4.fill_between(x, lbound, ubound, alpha=0.5, color='g')

        setlims(ax4,[0,max(x)],[0,1.2*max(lbound)])
        ## TODO: setlabels(ax4,'r [pc]','\sigma_{los,2} [km/s]')

    if gp.plotdens:
        #### density
        x = gp.rcore[0]*gp.dat.densx #[pc]
        den = gp.dens0pc[0]*gp.dat.densdat  #[Msun/pc^2]
        ax5.plot(x, den, c='green', lw=3)
        # ax5.plot(gp.dat.densx, phys.calculate_dens(M_anf(gp.dat.densx),gp.dat.densx),c='blue',lw=2)
        ## TODO: setlabels(ax5,'r [pc]','dens [Msun/pc^2]')
        #if gp.lim:   setlims(ax5,[0,max(x)/2],[1e-4,1.5*max(den)])
        if gp.model:
            rhodm, rhostar1, rhostar2 = rhowalker(gp.dat.Mx)
            x = gp.dat.Mx*gp.rcore[0]#[pc]
            ax5.plot(x, (rhodm+rhostar1+rhostar2)*gp.totmass[0], c='blue', lw=4)
            ax5.plot(x, (rhostar1+rhostar2)*gp.totmass[0],     c='blue', lw = 3, ls = '-.')
    else:
        #### mass
        x = gp.dat.Mx*gp.rcore[0] #[pc]
        lbound = (gp.dat.Mdat - gp.dat.Merr)*gp.totmass[0] #[munit]
        ubound = (gp.dat.Mdat + gp.dat.Merr)*gp.totmass[0] #[munit]
        ax5.fill_between(x, lbound, ubound, alpha=0.5, color='g')
        ax5.set_ylim([0.,1.2*max(ubound)])
        #ax5.plot(x, M_anf(gp.dat.Mx)*gp.totmass[0], c='blue', lw=4) # only for Hernquist
        ## TODO: setlabels(ax5,'r [pc]','M [Msun]')
        if gp.lim:   setlims(ax5,[0,max(x)],[1e-4,1.5*max(ubound)])

    if gp.log:   ax5.set_yscale('log')

    plt.draw()
    
def get_plot_data():
    x = gp.xipol
    r0,r_tot,dummy,dummy,dr,r_outer = phys.get_zarrays(x)
    densp = gp.ipol.densdat if gp.checksigma else gp.parst.dens
    M = phys.Mzdefault(densp)
    mprioru = gp.parst.Msl
    M_tot = phys.get_M( M, mprioru, r_outer, dr ) * gp.totmass[0]
    x *= gp.rcore[0]
    r0 *= gp.rcore[0]
    r_tot *= gp.rcore[0]
    dr *= gp.rcore[0]
    r_outer *= gp.rcore[0]

    if gp.analytic:
        M_tot = M_anf(r_tot)
        # sig1  = gp.ipol.sigdat1 # check: plot is really done as it should
        #mprioru  = (M[gp.nipol-1]-M[3*gp.nipol/4])/(r0[gp.nipol-1]-r0[3*gp.nipol/4])
        mprioru  = 0. # TODO: check right

    if gp.plotdens:
        r_tot = x
        M_tot = phys.densdefault(gp.parst.dens)
        if gp.checksigma:
            M_tot = gp.ipol.densdat*gp.dens0pc[0]  # density from data over gp.xipol, attention: only for spherical case
        
    x1, nu1, sig1 = units.get_physical(gp.xipol, gp.nu1_x, gp.sig1_x, 0)
    if gp.pops==1:
        return x1, r_tot, nu1, sig1, M_tot, gp.parst.delta1
    elif gp.pops==2:
        x2, nu2, sig2 = units.get_physical(gp.xipol, gp.nu2_x, gp.sig2_x, 1)
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
    ax6.set_ylim([-1.,1.])
    if (gp.pops == 2) : ax6.plot(xpl,delta2,color='orange')
    plt.draw()
    return

def update_plot():
    if not gp.testplot: return
    if gp.pops==1:
        xpl,x_tot,nu1,sig1,M_tot,delta1=get_plot_data()
    elif gp.pops==2:
        xpl,x_tot,nu1,sig1,M_tot,delta1,nu2,sig2,delta2=get_plot_data()
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
    ax6.set_ylim([-1.5,1.5])
    save_plot()
    plt.draw()

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
                
