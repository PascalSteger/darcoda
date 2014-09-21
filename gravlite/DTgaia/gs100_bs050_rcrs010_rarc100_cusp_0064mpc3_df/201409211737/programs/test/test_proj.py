#!/usr/bin/env ipython3

##
# @file
# check deprojection and projection of nu
# with analytic Hernquist profiles

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import sys
import numpy as np
from scipy.interpolate import splrep, splev, splint
from pylab import *

import import_path as ip
ip.insert_sys_path('/home/psteger/sci/darcoda/gravlite/programs/datareduction')
ip.insert_sys_path('/home/psteger/sci/darcoda/gravlite/programs/sphere')

#from gl_timing import *
import time
import gl_params as gp
gp.rinfty = 5
gp.nexp = 3
#import gr_params as grp
#import gl_file as gfile
import gl_helper as gh
import gl_int as gi
import gl_project as glp
import gl_analytic as ga

# unitsXS
G1  = 6.67398e-11                # [m^3 kg^-1 s^-2]
pc  = 3.08567758e16              # [m]
msun= 1.981e30                   # [kg]
km  = 1000.                      # [m]
G1  = G1*msun/km**2/pc

Nfine = 30
MtoL = 100000. # constant

def introduce_points_in_between(r0):
    rmin = np.log10(min(r0))
    rmax = np.log10(max(r0))
    return np.logspace(rmin, rmax, Nfine)
## \fn introduce_points_in_between(r0)
# get gp.fine points logarithmically spaced points
# @param r0 3D radius


r0 = np.logspace(-2, 2, 10)
nuanf  = ga.rho_hern(r0, 1, 1) # 1/MtoL
Siganf = ga.Sig_hern(r0, 1, 1) # 1/MtoL

minr = min(r0)
maxr = max(r0)
repol = xepol = np.hstack([minr/8.,minr/4.,minr/2.,r0,2*maxr, 4*maxr, 8*maxr])
rfine = introduce_points_in_between(r0)

Sigdatnu, Sigerrnu = gh.complete_nu(r0, Siganf, Siganf/20, rfine)
#loglog(r0, Siganf, 'b.-')
#loglog(rfine, Sigdatnu, 'r.-')


loglog(r0, nuanf, 'b.-')
nudatnu = glp.Sig_INT_rho(rfine, Sigdatnu, gp)
loglog(rfine, nudatnu, 'r.-', lw=0.5)
#nudatnuold = glp.Sig_INT_rho_buggy(rfine, Sigdatnu, gp) #negative found
#loglog(rfine, nudatnuold, 'g.-')
pdb.set_trace()

dummy, nudatnu, nuerrnu, Mrdatnu = glp.Sig_NORM_rho(rfine, Sigdatnu, Sigerrnu, gp)


#loglog(r0, nuanf, 'b.-')
#loglog(rfine, nudatnu, 'r.-', alpha=0.5)


Sigdatproj_coarse = glp.rho_INTIPOL_Sig(r0, nuanf, gp)
Sigdatproj_fine   = glp.rho_INTIPOL_Sig(rfine, nudatnu, gp)

# loglog(r0, Siganf, 'b.-')
# loglog(r0, Sigdatproj_coarse, 'r')
# loglog(rfine, Sigdatproj_fine, 'g')



betanu = np.zeros(len(rfine))
rhonu  = ga.rho_hern(rfine, 1, 1)
nunu   = ga.rho_hern(rfine, 1, 1/MtoL)
idnu   = np.zeros(len(rfine))

######################## calculation of sigma_LOS ##########################


start_time = time.time()

# integrate enclosed 3D mass from 3D density
r0tmp = np.hstack([0.,rfine])
rhoint = 4.*np.pi*rfine**2*rhonu
# add point to avoid 0.0 in Mrnu(rfine[0])
rhotmp = np.hstack([0.,rhoint])

splpar = splrep(r0tmp, rhotmp, k=1, s=0.) # not necessarily monotonic
Mrnu = np.zeros(len(rfine))              # work in refined model
for i in range(len(rfine)):              # get Mrnu
    Mrnu[i] = splint(0., rfine[i], splpar)
gh.checkpositive(Mrnu, 'Mrnu')

Mr_hern = ga.M_hern(rfine, 1, 1)

# (sigr2, 3D) * nu/exp(-idnu)
xint = rfine                           # [pc]
yint = G1 * Mrnu / rfine**2         # [1/pc (km/s)^2]
yint *= nunu                          # [Munit/pc^4 (km/s)^2]
yint *= np.exp(idnu)                  # [Munit/pc^4 (km/s)^2]
gh.checkpositive(yint, 'yint sigr2')

# use quadinflog or quadinfloglog here
sigr2nu = np.zeros(len(rfine))
for i in range(len(rfine)):
    # TODO: check quadinflog with walker profiles
    sigr2nu[i] = np.exp(-idnu[i])/nunu[i]*\
                 gh.quadinflog(xint, yint, rfine[i], gp.rinfty*rfine[-1], True)
    #if sigr2nu[i] == np.inf:
    #    sigr2nu[i] = 1e-100
    # last arg: warn if no convergence found
gh.checkpositive(sigr2nu, 'sigr2nu in sigl2s')

sigr2anf = ga.sigr2_hern(rfine, 1, 1, G1)


# project back to LOS values
# sigl2sold = np.zeros(len(rfine)-gp.nexp

### TODO: find error in here!

sigl2s = np.zeros(len(rfine)-gp.nexp)
for i in range(len(rfine)-gp.nexp): # get sig_los^2
    xnew = np.sqrt(rfine[i:]**2-rfine[i]**2)             # [pc]
    ynew = 2.*(1-betanu[i]*(rfine[i]**2)/(rfine[i:]**2)) # TODO check
    ynew *= nunu[i:] * sigr2nu[i:]
    gh.checkpositive(ynew, 'ynew in sigl2s') # is hit several times..
    # check sigr2nu: has too many entries of inf!

    # stop integration at xnew[-1] instead of at np.inf
    # to circumvent inf when nunu has increase at fudge radii
    sigl2s[i] = gh.quadinflog(xnew, ynew, 0, gp.rinfty*xnew[-1], False)
# for last 3 bins, we are up to a factor 2 off
gh.checkpositive(sigl2s, 'sigl2s')

# calculate surface density on the same rfine as the sigl2s
surfden = glp.rho_INT_Sig(rfine, nunu, gp)
siglos = np.sqrt(sigl2s/surfden[:-gp.nexp])

elapsed_time = time.time() - start_time
print("elapsed time: "+str(elapsed_time))

surfdenanf = ga.Sig_hern(rfine, 1, 1/MtoL)
rfine = rfine[:-gp.nexp]
siglanf = ga.sig_los_hern(rfine, 1, 1, G1)

plot(rfine, siglanf, 'b.-')
xscale('log')
plot(rfine, siglos, 'r.-')
pdb.set_trace()
