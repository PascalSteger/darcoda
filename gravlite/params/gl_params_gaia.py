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

investigate  = 'walk'  # determine which data set to work on
                       # 'discmock': set up simple model for disc
                       # 'discsim': read in disc simulation
                       # 'checkdwarf': check int for analytic dwarf, sig_LOS
                       # 'hern': check simple Hernquist prof. from simwiki
                       # 'walk': check with full obs. cont. data from Walker
                       # 'gaia': 6D data (x,y,z,vx,vy,vz) from gaia challenge, 1 pop only
                       # 'fornax': real data from Fornax dwarf galaxy

case = 2 # choose gaia models (1..8) or Walker (0..2,4,5) models

# Set number of tracer stars to look at
# take all particles                       # case 0
# want to set ntracer = 3e3              # case 1
#             ntracer = 1e4              # case 2
cas = 2

getnewdata = True # get new data computed from observations before burn-in
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
# TODO remove following scales for Gaia. set to 1 for no scaling anywhere..
Mscale = 1.                      # [Msun], scale for dimensionless eqs
                                   # from Hernquist, Baes&Dejonghe
ascale = 1.                      # [pc]

err = -1e30

model = False # for walker mock data: plot model. set beta to Osipkov-Merrit profile
if investigate != 'walk':
    model = False

########## integration options
densint    = False # use integral of dens to find mass, not binned sum
even       = 'avg' # for simps integration (everywhere): 'avg', 'first', 'last'
usekappa   = False # switch to turn on (True) or off the calculation of kappa

sim       = 2 # use hernquist model with 1 or 2 particle types. do not use second type (DM) as population
pops      = 1

global rstarhalf
    
# Set number of terms for enclosedmass+tracer+anisotropy bins = model parameters:
nipol = 7
nexp  = 3    # 3 more parameters for n(r) at 2, 4, and 8*rmax
nepol = nipol + nexp + 3 # +3 means 1 more parameter for the density at half light radius
                         #          1 more parameter for the asymptote to 0
                         #          1 more parameter for the asymptote to \infty
rinfty = 10. # interpolate from last slope to slope at 10*max(xipol), where asymptote to \infty is reached
nbeta = 2   # number of parameters for beta, in sum of polynomials
# next: # live points, > ndim, < 2^ndim, about number of ellipsoids in phase space to be found
nlive = 2*(nepol + pops*nepol + pops*nbeta)
# if more than nlive points trigger the same penalty value => Multinest finishes

########## priors
bprior   = True    # Baryon minimum surfden prior: rho_*=sum_i nu_i only for mock!
Mmodel = np.zeros(nipol)
baryonmodel = 'sim' # read in surface density from corresponding surfden file

delta0 = np.zeros(nipol)
# sigmaprior, default: -1
sigmaprior1 = 0.3; sigmaprior2 = 0.3

sprior  = False                       # rising sig_z needed in disc case
if geom=='sphere': sprior = False
constdens = False # constant DM density
rprior  = True    # regularize Nuz
nutol   = 1.0     # (nu_(i+1) - nu_i) must be < nutol * nu_(i+1)
rhotol  = 1.0
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

global xipol, xepol
xipol = np.array([]); xepol = np.array([])
from gl_class_files import Files
files = Files()

maxrhoslope = 5.
maxnuslope = 5.
maxbetaslope = 2.

nuerrcorr = 1.   # 2.   # scale error from grw_dens by this amount
sigerrcorr = 1.  # 1.5  # best done directly in the corresponding grw_dens, grw_siglos
                 # thus, both times 1. is best
kaperrcorr = 1.


########## global variables, not set to value here
global pars, parst, parstep
global dat
global nu1_x, nu2_x, d1_x, d2_x, sig1_x, sig2_x, M_x, rho_x, M_tot, rho_tot
global fnewoverf
global zmin, zmax    # Low/high-z range = min/max of data [-1 = default]: TODO: convert to xmin, xmax
xpmin = -1;  xpmax = -1                 # Default low/high-r range = min/max of data: 

prob = -1e6 if uselike else 1e6     # Initial log likeli. 
# [small in case we need it in disc]
dof = 1

ini     = 0                    # TODO: meaning
inimax  = 5000                 # 
counter = 0                    #

# Units: 
# rscale in [pc], surfdens_central (=Nu0) in [munit/rscale**2],
# and in [munit/pc**2], and totmass [munit], and max(v_LOS) in [km/s]

rscale=[]; nu0rscale=[]; nu0pc=[]; totmass=[]; maxvlos=[] # unit system
Rscale=[]; Nu0rscale=[]; Nu0pc=[]
if investigate != 'walk' and investigate != 'triax' and investigate != 'gaia':
    # TODO: adapt to physical units
    # each is set for all components and first component by default
    rscale.append(1.);      Rscale.append(1.)
    rscale.append(1.);      Rscale.append(1.)
    nu0rscale.append(1.);   Nu0rscale.append(1.)
    nu0rscale.append(1.);   Nu0rscale.append(1.)
    nu0pc.append(1.);       Nu0pc.append(1.)
    nu0pc.append(1.);       Nu0pc.append(1.)
    totmass.append(1.)
    totmass.append(1.)
    maxvlos.append(1.)
    maxvlos.append(1.)

