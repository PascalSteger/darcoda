#!/usr/bin/env python3

##
# @file
# all parameters for the gravlite MCMC, gaia investigation

# (c) 2013 ETHZ Pascal S.P. Steger

import numpy as np
import sys
import os
import pdb
import logging

loglevel = logging.WARNING
logging.basicConfig(stream=sys.stderr, level=loglevel)
LOG = logging.getLogger(__name__)

machine = 'darkside'
if not os.path.exists("/home/ast/read"):
    machine = 'local'        # machine: 'local' or 'darkside'

investigate  = 'gaia'  # determine which data set to work on
                         # 'simple': set up simple model for disc
                         # 'sim': read in disc simulation
                         # 'checkdwarf': check int for analytic dwarf, sig_LOS
                         # 'hernquist': check simple Hernquist prof. from simwiki
                         # 'walker': check with full obs. cont. data from Walker
                         # 'gaia': work on xyzvxyz data from gaia challenge
                         # 'fornax': real data from Fornax dwarf galaxy

case = 1           # choose different Walker models (0-2 so far)

# Set number of tracer stars to look at
# take all particles                       # case 0
# want to set ntracer = 3e3              # case 1
#             ntracer = 1e4              # case 2
cas = 1

getnewdata = False # get new data computed from observations before burn-in
metalpop   = False # split metallicities with a separate MCMC
lograd  = False # log steps for radial bin in readout, show x-axis in log scale
consttr = True # set radial bin by constant number of tracer particles

# units
G1  = 6.67398e-11                             # [m^3 kg^-1 s^-2]
pc  = 3.08567758e16                           # [m]
msun= 1.981e30                                # [kg]
km  = 1000.                                   # [m]
G1  = G1*msun/km**2/pc                        # [pc msun^-1 km^2 s^-2]

global geom
geom = 'sphere'
# TODO remove following scales for Gaia challenge. set to 1 for no scaling anywhere..
Mscale = 1.                      # [Msun], scale for dimensionless eqs
                                 # from Hernquist, Baes&Dejonghe
ascale = 1.                      # [pc]

########## plotting options
showplot   = True                       # show plots?
plotdens   = True # plot dens instead of M or Sigma_z in lower left plot
lim        = False         
log = True if plotdens and geom == 'sphere' else False
# ^--  adjust range of plots to predetermined values


########## checking options
checkint = False # check sigma_los integration
# (set nu, M to interpolated data values, only run for one iteration)

analytic   = False # calc sig_los from analytic Hernquist profiles for nu, M
if investigate != 'hernquist': analytic = False

model = True # for mock data: plot model. set beta to Osipkov-Merrit profile
if investigate != 'walker' and investigate != 'gaia': model = False


########## density options
poly       = True  # use polynomial representation of dens during init
if analytic: poly = False

densstart = -0.9 # -2.6 for Hernquist, -2.3 for Walker cusp, -1.8 for core
scaledens = 1. # fraction of max radius from data, for which the poly is scaled
scalepower = 1.8 # 0. 9 5 for Hernquist, 1.3 for Walker cusp, 2.2 for core


########## integration options
densint    = False # use integral of dens to find mass, not binned sum
even       = 'avg' # for simps integration (everywhere): 'avg', 'first', 'last'
usekappa   = False # switch to turn on (True) or off the calculation of kappa

pops      = 1
if analytic: pops = 1 # only work with 1 pop if analytic in hernquist case is set



# Set number of terms for enclosedmass+tracer+anisotropy models:   
nipol = 12

# TODO: better: set scale for central part (first and second bin),
# then determine number of particles inside: ninside = len(r[r<centralscale])
# then get nipol = min(ntotal,ntracer)/ninside

nconstraints = 4                        # number of constraints
dof = 1. # degrees of freedom: 1 for prop. chi2, 6*nipol - nconstraints # for reduced chi2

########## priors
gprior = -1                          # Gradient prior [-1 = no prior]:
cprior = -1 # Central pixel prior [prior for first point in M_r; def= 0.]:
if cprior<0:
    cprior = 1e30

# Convert gprior into appropriate units:
if gprior > 0: gpriorconv = gprior * (2.*np.pi*G1) * 1000.**2. 
else: gpriorconv = 1e30
if cprior >= 0: cpriorconv = cprior * (2.*np.pi*G1) * 1000.**2.
else: cpriorconv = 1e30

bprior   = False                       # Baryon minimum surfden prior
blow = np.zeros(nipol);   Mmodel = np.zeros(nipol)
baryonmodel = 'sim' # read in surface density from corresponding surfden file

mirror   = False # Mirror prior: TODO
nulog    = True  # sample nu (only) in logarithmic space.
denslog  = True  # after init: sample dens (only) in logarithmic space
mprior = -1      # Mass prior
deltaprior  = False # Deltaprior: - beta (velocity anisotropy) in spherical
                    #             - tilt in disc geometry
d1wild = False; d2wild = False # show message if jumpy delta, only once
dens2wild = False; sig2wild = False
# nu2wild gives Nr of possible nutol prior violations before fallback
b2wild = False; nu2wild = 1000    
lasterr = 'None'

delta0 = np.zeros(nipol)
# sigmaprior, default: -1
sigmaprior1 = 0.3; sigmaprior2 = 0.3

sprior  = False                       # rising sig_z needed in disc case
if geom=='sphere': sprior = False
constdens = False # constant DM density
rprior  = True    # regularize Nuz 
nutol   = 1.0     # (nu_(i+1) - nu_i) must be < nutol * nu_(i+1)
ktol    = 0.      # same as for nu, but for dens, 50% up is still fine
deltol  = 2./nipol                               # for delta
# norm1   = 17.**2 # offset of sigma[0]/nu[0], from int starting at zmin instead of 0
# norm2   = 10.**2 # and for the second component, if there is one
quadratic = False           # linear or quad interpol. 
monotonic = False           # mono-prior on nu(z)
uselike   = False           # use Likelihood function, or binned data?
adderrors = False

# last bin mass prior:
# exclude all models with a mass in last bin that exceeds
# lbtol times the integrated mass of all bins before
lbprior = False; lbtol = 0.33




########## file options
run_configs = []
init_configs = []

from gl_class_files import Files
files = Files()


### safe defaults: first guess for init, end of init for real run, to be called if 100 nutol error occurred
from gl_class_params import *
safepars = Params(0)
safeparstep = Params(0)
safechi2 = 1e6

nuerrcorr = 1.   # 2.   # scale error from grw_dens by this amount
sigerrcorr = 1.  # 1.5  # best done directly in the corresponding grw_dens, grw_siglos
                 # thus, both times 1. is best
kaperrcorr = 1.


########## global variables, not set to value here
global pars, parst, parstep
global dat, ipol, xipol
global nu1_x, nu2_x, d1_x, d2_x, sig1_x, sig2_x, M_x, dens_x, M_tot, dens_tot
global fnewoverf
global zmin, zmax    # Low/high-z range = min/max of data [-1 = default]: TODO: convert to xmin, xmax
xpmin = -1;  xpmax = -1                 # Default low/high-r range = min/max of data: 




########## MCMC parameters
niter = 1000000                  # Maximum number of iterations
# TODO: class for chi2
chi2   = 1e6                 # initial chi2 [large]
chi2t  = 1e6
chi2t1 = 1e6;    chi2t2 = 1e6
chi2t_nu  = 1e6; chi2t_sig  = 1e6; chi2t_kap = 1e6
chi2t_nu1 = 1e6; chi2t_nu2  = 1e6
chi2t_sig1= 1e6; chi2t_sig2 = 1e6
chi2t_kap1= 1e6; chi2t_kap2 = 1e6

prob = -1e6 if uselike else 1e6     # Initial log likeli. 
# [small in case we need it in disc]

ini     = 0                    # TODO: meaning
inimax  = 5000                 # 
counter = 0                    #

rollsize = 20 # number of stacked acceptance/rejection rates
              # as well as number of steps before next step adaption done
adaptstepwait = rollsize
from gl_class_rate import Rate
# new acceptance/rejection rate object, mean 0.25, acceptable +-5%
accrate = Rate(0.25, 0.10)
chi2tol = 1.*(6*nipol-nconstraints) 
# ^--- tolerance in (chi2) or (without paranteses) reduced chi2
endcount = 3*rollsize # 30*
# ~900 accepted models which chi2<chi2tol means initialization phase is over
# better measure: 1./(min stepsize),
# as this gives the time needed to get convergence on this parameter

stepcorr= 1.1  # factor to adapt stepsize if not 0.24 < acc/rec < 0.26
farinit = 8. # 5 times chi2 is too far off in init phase: start new from last point
stepafterrunaway = 0.9 # mult. stepsize by this amount if too low fnewoverf 2.5
farover = 5.         # 5 times chi2 is too high after init phase
scaleafterinit   = 1.0 # <= cheat: multiply stepsize by this amount if init is over

# Parameters to end initphase 
initphase = True # initialisation phase flag, first True, if over: False
endgame  = False # Ending flag


# Units: 
# rcore in [pc], surfdens_central (=dens0) in [munit/rcore**2], 
# and in [munit/pc**2], and totmass [munit], and max(v_LOS) in [km/s]

rcore=[]; dens0rcore=[]; dens0pc=[]; totmass=[]; maxvlos=[] # unit system
rcore_2D=[];dens0rcore_2D=[];dens0pc_2D=[]
if investigate != 'walker' and investigate != 'triaxial' and investigate != 'gaia':
    # TODO: adapt to physical units
    # each is set for all components and first component by default
    rcore.append(1.);      rcore_2D.append(1.)
    rcore.append(1.);      rcore_2D.append(1.)
    dens0rcore.append(1.); dens0rcore_2D.append(1.)
    dens0rcore.append(1.); dens0rcore_2D.append(1.)
    dens0pc.append(1.);    dens0pc_2D.append(1.)
    dens0pc.append(1.);    dens0pc_2D.append(1.)
    totmass.append(1.)
    totmass.append(1.)
    maxvlos.append(1.)
    maxvlos.append(1.)

