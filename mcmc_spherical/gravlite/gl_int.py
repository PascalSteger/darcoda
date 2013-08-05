#!/usr/bin/python
# (c) 2013 Pascal S.P. Steger
'''all integrals from physics_sphere'''

import numpy as np
import scipy
# import scipy.integrate
from scipy.integrate import simps,trapz
import pdb

import gl_params as gp
import gl_funs as gfun
import gl_helper as gh
import gl_plot as gpl




# def int_rho(rho,r):
#     'solve mass integral'
#     tmp = np.zeros(len(r))
#     ri  = np.hstack([0,r,r[-1]+(r[-1]-r[-2])])
#     rhoi= gh.ipollog(r,rho,ri)
#     # tmp[0] = 4.*np.pi/3.*rho[0]*r[0]**3
#     for i in range(2,len(ri)):
#         xint = ri[:i]
#         yint = rhoi[:i]*ri[:i]**2
#         tmp[i-2] = simps(yint,xint,even=gp.even)
#     return 4.*np.pi*tmp # *0.92, factor needed to account for overestimation of simpson vs rectangles





# TODO: correct first bin
def int_delta_r_t(pnts, rtot, deltatot): #[1], [lunit], [1]
    'integrate beta (b) over r (integrals in front of and after sigma_r^2 integral)'
    tmp = np.zeros(2*pnts-1)
    for i in range(1,2*pnts-1):
        xint = rtot[:i]                             # [lunit]
        yint = deltatot[:i]/rtot[:i]                # [1/lunit]
        tmp[i] = 2.*simps(yint, xint, even=gp.even) # [1]
    return tmp                                      # [1]



def ant_delta_r_t(r0, delta): #[1], [lunit], [1]
    'integrate beta (b) over r (integrals in front of and after sigma_r^2 integral)'
    tmp = np.zeros(len(delta))
    for i in range(1,len(delta)):
        xint = r0[:i]                               # [lunit]
        yint = delta[:i]/r0[:i]                     # [1/lunit]
        tmp[i] = 2.*simps(yint, xint, even=gp.even) # [1]
    return tmp                                      # [1]






def int_sigr2(pnts, r_tot, delta_tot, M_tot, nu_tot): #  [1], [pc], [1], [munit], [munit/pc^3]
    'sigma_r^2'
    # the integrand is assumed to be zero above 2*rmax
    tmp = np.zeros(2*pnts-1)
    for i in range(2*pnts-1):
        xint  = r_tot[i:]               # [pc]
        yint  = gp.G1 * M_tot[i:]/r_tot[i:]**2 # [1/pc  munit/msun km^2/s^2] = [1/pc (km/s)^2]
        yint *= nu_tot[i:]                     # [munit/pc^4 (km/s)^2]
        yint *= np.exp(delta_tot[i:])          # [munit/pc^4 (km/s)^2]
        tmp[i] = np.exp(-delta_tot[i]) * simps(yint, xint, even=gp.even)/nu[i] # [(km/s)^2]

    # gp.LOG.info( 'compute sigma outer')
    # could calculate that from last integral already, trying to extrapolate anyhow
    # if (sigmaprior < -1): sigmaprior=-1.
    # if gp.checksigma: sigmaprior=-1.
    # nusig2_outer=nusig2[-1]+(np.arange(1,pnts))*sigmaprior*nusig2[pnts-1]/(xmax-xmin)*dx
    # nusigma2_tot=np.hstack((nusig2,nusig2_outer))
    # sigr2_tot = gh.ipollog(r_tot[:-1],sigr2[:-1],r_tot) # beware of interpolations, dude!
    tmp[-1] = tmp[-2]           # [(km/s)^2]
    return tmp                  # [(km/s)^2]


def int_sigr2(pnts, r_tot, delta_tot, M_tot, nu_tot): #  [1], [pc], [1], [munit], [munit/pc^3]
    'sigma_r^2'
    # the integrand is assumed to be zero above 2*rmax
    tmp = np.zeros(2*pnts-1)
    for i in range(2*pnts-1):
        xint  = r_tot[i:]               # [pc]
        yint  = gp.G1 * M_tot[i:]/r_tot[i:]**2 # [1/pc  munit/msun km^2/s^2] = [1/pc (km/s)^2]
        yint *= nu_tot[i:]                     # [munit/pc^4 (km/s)^2]
        yint *= np.exp(delta_tot[i:])          # [munit/pc^4 (km/s)^2]
        tmp[i] = np.exp(-delta_tot[i]) * simps(yint, xint, even=gp.even)/nu[i] # [(km/s)^2]

    # gp.LOG.info( 'compute sigma outer')
    # could calculate that from last integral already, trying to extrapolate anyhow
    # if (sigmaprior < -1): sigmaprior=-1.
    # if gp.checksigma: sigmaprior=-1.
    # nusig2_outer=nusig2[-1]+(np.arange(1,pnts))*sigmaprior*nusig2[pnts-1]/(xmax-xmin)*dx
    # nusigma2_tot=np.hstack((nusig2,nusig2_outer))
    # sigr2_tot = gh.ipollog(r_tot[:-1],sigr2[:-1],r_tot) # beware of interpolations, dude!
    tmp[-1] = tmp[-2]           # [(km/s)^2]
    return tmp                  # [(km/s)^2]



def int_poly_inf(r0,poly):
    '''integrate polynomial from r0[i] to infinity'''
    # assume log(igra) = A + B*r;    B = poly[0] < 0, A = poly[1]
    f = -1/poly[0]*np.exp(poly[1]+poly[0]*r0)
    return f


def ant_sigr2(r0, intdelta, M, nu): #  [1], [pc], [1], [munit], [munit/pc^3]
    'sigma_r^2'
    if min(nu)<0.:
        print 'error: negative nu'
        pdb.set_trace()

    # get full integrand first
    igra = np.log(gp.G1 * M * nu/r0**2 * np.exp(intdelta))

    # then fit polynomial
    x = r0[-gp.nipol/2:]
    y = igra[-gp.nipol/2:]
    # TODO: check that the following is not wiggling too much for (not checksigma)
    # x = r0[-2:]
    # y = r0[-2:]
    polyhilo = np.polyfit(x,y, 1)#, rcond=None, full=False, w=None,cov=False) # 3 is degree, w is weights

    # integrate polynomial from r0 (up to max(r_data)) to infinity
    intpoly = int_poly_inf(r0[-1],polyhilo)

    # the integrand is assumed to be zero above 2*rmax
    tmp = np.zeros(gp.nipol)
    for i in range(gp.nipol):
        xint  = r0[i:]               # [pc]

        yint  = gp.G1 * M[i:]/r0[i:]**2 # [1/pc munit/msun km^2/s^2] = [1/pc (km/s)^2]
        yint *= nu[i:]                     # [munit/pc^4 (km/s)^2]
        yint *= np.exp(intdelta[i:])          # [munit/pc^4 (km/s)^2]

        # case a: add integrated polynomial above 3r_S only
        tmp[i] = (np.exp(-intdelta[i])/nu[i]) * \
                 (simps(yint, xint, even=gp.even)  +  intpoly) # [(km/s)^2]

        # case b: return integrated polynomial directly.
        # forget case b, only valid if polynomial fit over whole of gp.xipol fine
        # but: we have exponential decay in semilog igra, so only take 2nd half for approx.
        # and have to forget all the inner values from extrapolation

    return tmp                  # [(km/s)^2]




    





def ant_siglos2surf(r0, delta, nu, sigr2): # [1], [pc], [1], [munit/pc^3], [(km/s)^2], gives back [TODO]
    'take nu and sig_r^2, give back sigma_LOS, with analytical integration to infinity'
    pnts = len(r0)-1
    tmp = np.zeros(pnts)

    for i in range(pnts):
        xint = r0[i:]                          # [pc]
        yint = (1-delta[i:]*(r0[i]/r0[i:])**2) # [1]
        yint = yint * nu[i:] * sigr2[i:] * r0[i:] # [munit/pc^2 (km/s)^2]
        yint = yint / np.sqrt(r0[i:]**2-r0[i]**2) # [munit/pc^3 (km/s)^2]

        # interpolate diverging bin, need 4 (sub)points to get 3rd order
        xipol = xint[0]+np.arange(4)*(xint[1]-xint[0])/4.
        dipol = delta[i]+np.arange(4)*(delta[i+1]-delta[i])/4.
        nuipol= nu[i]+np.arange(4)*(nu[i+1]-nu[i])/4.
        sigr2ipol = sigr2[i]+np.arange(4)*(sigr2[i+1]-sigr2[i])/4.
        yipol = (1-dipol*(r0[i]/xipol)**2) * nuipol\
                * sigr2ipol * xipol / np.sqrt(xipol**2-r0[i]**2)
        yint[0] = gh.ipol(xipol[1:4],yipol[1:4],xint[0])

        # yint[0] = gh.ipol(xint[1:3],yint[1:3],[xint[0]]) # -(xint[1]-xint[0]))        
        # ^-- rather conservative, lower than above method by factor 0.6

        # then fit polynomial on second half of radii
        ml = max(len(xint)/2,2)
        x = xint[-ml:]
        y = np.log(yint[-ml:])
        polyhilo = np.polyfit(x,y, 1)

        # integrate polynomial from r0 (up to max(r_data)) to infinity
        intpoly = int_poly_inf(r0[-1],polyhilo)

        tmp[i] = 2. * (simps(yint, xint, even=gp.even) + intpoly) # [munit/pc^2 (km/s)^2]

    # extrapolate to last 4 bins
    x = r0[-4:-1]                 # take values from half+some offset to last point known for tmp
    y = np.log(tmp[-3:])
    polyhilo = np.polyfit(x,y,1)
    
    return np.hstack([tmp,np.exp(r0[-1:]*polyhilo[0]+polyhilo[1])])    # [munit/pc^2 (km/s)^2]









def int_surfden(r0, nu):           # [lunit], [munit/lunit^3]
    'compute surface density, assume len(r) = 2*pnts, or at least pnts+4'
    pnts = len(r0)-1 # start with one missing bin, s.t. interpolation on sub-bin possible
    tmp = np.zeros(pnts)
    for i in range(pnts):
        xint = r0[i:]                                      # [lunit]
        yint = r0[i:] * nu[i:]/np.sqrt(r0[i:]**2-r0[i]**2) # [munit/lunit^3]

        # interpolate diverging bin, need 4 (sub)points to get 3rd order
        xipol = xint[0]+np.arange(4)*(xint[1]-xint[0])/4.
        nuipol= nu[i]+np.arange(4)*(nu[i+1]-nu[i])/4.
        yipol = xipol * nuipol/np.sqrt(xipol**2-r0[i]**2)
        yint[0] = gh.ipol(xipol[1:3],yipol[1:3],xint[0])#*0.5 # TODO: tweaking

        # then fit polynomial on second half of radii
        ml = max(len(xint)/2,2)
        x = xint[-ml:]
        y = np.log(yint[-ml:])
        polyhilo = np.polyfit(x,y, 1)

        # integrate polynomial from r0 (up to max(r_data)) to infinity
        intpoly = int_poly_inf(r0[-1],polyhilo)

        # attention: can be that intpoly is negative, if nu[i+1] > nu[i]
        # and even that integration + intpoly < 0 in the wash for tmp[i] for the last bin
        # circumvent that by requesting intpoly in [0,infty[
        if intpoly<0.: intpoly = 0.

        tmp[i]  = 2.*(simps(yint, xint, even=gp.even) + intpoly) # [munit/lunit^2]

    # better: extrapolation using polynomes
    # surfden = gh.ipollog(r_tot[1:-1],surfden[1:-1],r_tot)
    # tmp[0]  = gh.ipol(r0[1:3],tmp[1:3],[r0[0]]) # TODO: check that it's fine to disable


    # extrapolate to last 4 bins, which were neglected for sake of extrapolating yint[0]
    x = r0[-4:-1]                 # take values from half+some offset to last point known for tmp
    y = np.log(tmp[-3:])          # [log(munit/lunit^2)]
    polyhilo = np.polyfit(x,y,1)

    return np.hstack([tmp,np.exp(r0[-1:]*polyhilo[0]+polyhilo[1])])    # [munit/lunit^2]









def int_project(r0, rho): #[lunit], [dens0, 3D] TODO
    'take 3D (mass) density, convert it to surface density, and give back enclosed mass in rings'
    surf_tot = int_surfden(r0, rho) #gives [dens0, 2D]
    if gp.geom == 'disc':
        print 'attention: using spherical part of code for disc!'
        pdb.set_trace()
    import physics_sphere as phys
    surfmass = phys.Mr2D(r0, surf_tot)
    return surfmass #[munit, 2D]





def exto0(r, var):
    'extrapolate binned data to radius 0'
    if r[0] == 0.:
        return r, var
    var0 = gh.ipol(r[1:3],var[1:3],[0])
    var = np.hstack([var0, var])
    r = np.hstack([0, r])
    return r, var





# TODO: delete, not needed anymore, use phys.Mr2D(rdat, dens)
def int_projmass(pnts, r, Sigma): #[rcore], [dens0]
    'take surface density, calculate mass in rings' #TODO
    # there is another one doing this already, in fixed bins.
    tmp = np.zeros(pnts)
    r, Sigma = exto0(r, Sigma)          # get first bin at 0
    
    for i in range(pnts):
        xint = r[:i] #[rcore]
        yint = r[:i]*Sigma[:i] #[rcore*dens0]
        tmp[i] = 2. * np.pi * simps(yint,xint,even=gp.even) #[munit]
    return tmp #[munit]



def ext4log(r,rho):
    'extend array rho(r) by 4 entries, from last dr, with decaying rho extrapolated'
    dr = r[-1]-r[-2];
    rext4 = [r[-1]+dr, r[-1]+2.*dr, r[-1]+3.*dr, r[-1]+4.*dr]

    # constant decrease after end
    rhoext4 = [rho[-1]*0.9, rho[-1]*0.8, rho[-1]*0.7, rho[-1]*0.6]

    # rhoext4 = gh.ipol(r,rho,rext4) #direct extrapolation

    # log continuation
    #gom = rho[-2]/rho[-1]
    #rhoext4 = [rho[-1]/gom, rho[-1]/gom**2, rho[-1]/gom**3, rho[-1]]/gom**4;

    r = np.hstack([r, rext4])
    rhoext = np.hstack([rho, rhoext4])
    return r,rhoext


def int_2D3D(r, nu2d):                   # [munit/lunit^2], [lunit]
    'take surface density, deproject, 2D => 3D with radius r'
    # we miss all contributions from outside rmax! => TODO: check from data, how much that is
    pnts = len(nu2d)                    # [1]

    # assume linear change in r, linear decay in log space for nu
    r, nu2d = ext4log(r,nu2d)           # [lunit], [munit/lunit^2]
    dnubydR = gh.derivcoarse(nu2d,r)      # [munit/lunit^3] # use deriv, derivipol, or derivcoarse
    # correct last bin (lest it goes positive)
    dnubydR[-1] = dnubydR[-2]

    nu3d = np.zeros(pnts)
    for i in range(pnts):
        xint = r[i:]                                 #[lunit]
        yint = dnubydR[i:]/np.sqrt(r[i:]**2-r[i]**2) #[munit/lunit^3/lunit]
        yint[0] = gh.ipol(xint[1:4],yint[1:4],xint[0]) # TODO: too low point, could be 2x higher
        nu3d[i] = -1./np.pi * simps(yint,xint,even = gp.even) #[munit/lunit^3]
    return nu3d #[munit/lunit^3]



    
def int_siglos2surf(pnts, r0, delta, nu, sigr2): # [1], [pc], [1], [munit/pc^3], [(km/s)^2], gives back [TODO]
    'take nu and sig_r^2, give back sigma_LOS'
    tmp = np.zeros(pnts)
    for i in range(pnts):
        xint = r0[i:]                   # [pc]
        yint = (1-delta[i:]*(r0[i]/r0[i:])**2) # [1]
        
        yint = yint * nu[i:] * sigr2[i:] * r0[i:]  # [munit/pc^2 (km/s)^2]
        yint = yint / np.sqrt(r0[i:]**2-r0[i]**2) # [munit/pc^3 (km/s)^2]

        yint[0] = gh.ipol(xint[1:3],yint[1:3],[xint[0]]) #-(xint[1]-xint[0]))        

        tmp[i] = 2. * simps(yint, xint, even=gp.even) # [munit/pc^2 (km/s)^2]


    # tmp[0] = tmp[1]    # crude approximation: first bin has same sigma_los2 as second bin
    # tmp[0] = gh.ipol(r_tot[1:],tmp[1:],r_tot[0])
    # tmp[0] = 2.0/surfden[0]*nu_tot[0:]*sigma_tot[0:]*r_tot[0:]/r_tot[:]*rmin

    return tmp                          # [munit/pc^2 (km/s)^2]
