#!/usr/bin/python2.7
# (c) 2013 Pascal S.P. Steger
'''parameters for the gravlite MCMC'''

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




investigate  = 'walker'  # determine which data set to work on
                                        # 'simple': set up simple model for disc
                                        # 'sim': read in disc simulation
                                        # 'checkdwarf': checksigma for analytic dwarf, sig_LOS
                                        # 'hernquist': check simple Hernquist prof. from simwiki
                                        # 'walker': check with full obs. cont. data from Walker
                                        # 'fornax': real data from Fornax dwarf galaxy

walkercase = 1           # choose different Walker models (0-2 so far)

# Set number of tracer stars to look at in Hernquist profile
# take all particles                       # case 0
# want to set ntracers1 = 1e3              # case 1
#             ntracers1 = 1e4              # case 2
#             ntracers1 = ntracers2 = 5e3  # case 3
cas = 2
getnewdata = False       # get new data computed from observations


global geom
if investigate == 'sim' or investigate == 'simple':
    geom = 'disc'
    G1 = 6.67398e-11 * 1.989e30 / 3.086e19 / 1000.**2. #[kpc (km/s)^2 1/Msun]
    # means: measure mass in [Msun], distances in [kpc], velocities in [km/s]
    addpoptwo   = False  # add a second population if investigate = 'simple'
    adddarkdisc = False   # add a disk of DM particles in simple model
    slicedata   = False    # slice and dice data for sim case
    nusimpstart = True   # start from nu near data
    kzsimpfix   = False    # calculate kz from simple model directly
    nusimpfix   = False    # TODO: meaning?
    qtest = False # plot during physics_disc
    # Min/max Kz for MCMC search [only affects logprior].
    # If positive assume constant;
    # if negative take fraction of local baryonic value for that bin: 
    # kzmin = 0.0075 * (4.0 * !PI * G1) * 1000.^3. # kzmax = 100.*kzmin
    kzmin = 0.0;  kzmax = 1e5
    numin = 1e-3; numax = 1.
    patch = '0' # or 180 or ... for disc_sim case
else:
    geom = 'sphere'
    G1  = 6.67398e-11                       # [m^3 kg^-1 s^-2]
    pc  = 3.08567758e16                     # [m]
    msun= 1.981e30                          # [kg]
    km  = 1000.                             # [m]
    G1  = G1*msun/km**2/pc                  # [pc msun^-1 km^2 s^-2]
    Mscale = 1.e6                           # [Msun], scale for dimensionless eqs
                                            # from Hernquist, Baes&Dejonghe
    ascale = 1000.                          # [pc]

########## plotting options
testplot   = True          # show plots?
testplot_dwarfs = False    # show plots for dwarfs pos/com/density/siglos as well before init?
fplot      = 1./1.    # only plot every 1/Nth plot, in a random fashion
plotdens   = True          # plot dens instead of M or Sigma_z in lower left plot
lim        = False         
#log       = False          # yaxis in (M/dens) scaled in log?
log = True if plotdens and geom == 'sphere' else False
# ^--  adjust range of plots to predetermined values


########## checking options
checksigma = False # check sigma_los integration
# (set nu, M to interpolated data values, only run for one iteration)

analytic   = False         # calc sig_los from analytic Hernquist profiles for nu, M
if investigate != 'hernquist': analytic = False

model      = True # for Walker mock data: plot model
if investigate != 'walker': model = False


########## density options
poly       = True              # use polynomial representation of dens during init
if analytic: poly = False

densstart = -1.8              # -2.6 for Hernquist, -2.3979 for Walker
scaledens = 1. # percentage of maximum radius from data, for which the poly is scaled
scalepower = 2.4                  # 0.95 for Hernquist, 1.5 for Walker


########## integration options
densint    = False         # use integral of dens to find mass, not binned sum
even       = 'avg'        # for simps integration (everywhere): 'avg', 'first', 'last'




pops      = 2
if analytic: pops = 1 # only work with 1 pop if analytic in hernquist case is set





# Set number of terms for enclosedmass+tracer+anisotropy models:   
nipol = 12






########## priors
gprior = -1                             # Gradient prior [-1 = no prior]:
cprior = -1                           # Central pixel prior [prior for first point in M_r; def= 0.]:
if cprior<0:
    cprior = 1e30

# Convert gprior into appropriate units:
if gprior > 0: gpriorconv = gprior * (2.*np.pi*G1) * 1000.**2. 
else: gpriorconv = 1e30
if cprior >= 0: cpriorconv = cprior * (2.*np.pi*G1) * 1000.**2.
else: cpriorconv = 1e30

bprior   = False                       # Baryon minimum surfden prior
blow = np.zeros(nipol)
baryonmodel = 'sim'                    # read in surface density from corresponding surfden file

mirror   = False                       # Mirror prior: TODO
logprior = False                       # Logarithmic prior: sample in logarithmic space
nulog    = False                        # sample nu only in logarithmic space
mprior = -1                            # Mass prior
deltaprior  = False # Deltaprior: - beta (velocity anisotropy) in spherical
                    #             - tilt in disc geometry
delta0 = np.zeros(nipol)
# sigmaprior, default: -1
sigmaprior1 = 0.3; sigmaprior2 = 0.3

sprior  = False                   # rising sig_z, not recommended
if geom=='sphere': sprior = False
constdens = False                 # constant DM density
rprior  = True                   # regularize Nuz 
nutol   = 1.0     # (nu_(i+1) - nu_i) must be < nutol * nu_(i+1)
ktol    = 0.      # same as for nu, but for dens, 50% up is still fine
norm1   = 17.**2 # offset of sigma[0]/nu[0], from int starting at zmin instead of 0
quadratic = False                 # linear or quad interpol. 
monotonic = False                 # mono-prior on nu(z)
uselike   = False          # use Likelihood function, or binned data?
adderrors = False

# last bin mass prior:
# exclude all models with a mass in last bin that exceeds
# lbtol times the integrated mass of all bins before
lbprior = False; lbtol = 0.33




########## file options
run_configs = []
init_configs = []

import gl_class_files as gcf
files = gcf.Files()

nuerrcorr = 1.   # 2.   # scale error from grw_dens by this amount
sigerrcorr = 1.  # 1.5  # best done directly in the corresponding grw_dens, grw_siglos
                 # thus, both times 1. is best



########## global variables, not set to value here
global pars,  lpars
global parstep, parst, lparst, lparstep
global dat, ipol, xipol
global nu1_x, nu2_x, d1_x, d2_x, sig1_x, sig2_x, M_x, dens_x, M_tot, dens_tot
global fnewoverf
global zmin, zmax    # Low/high-z range = min/max of data [-1 = default]: TODO: convert to xmin, xmax
xpmin = -1;  xpmax = -1                 # Default low/high-r range = min/max of data: 




########## MCMC parameters
niter = 100000                  # Maximum number of iterations
# TODO: class for chi2
chi2   = 1e300                 # initial chi2 [large]
chi2t1 = 1e300;    chi2t2 = 1e300;    chi2t = 1e300
chi2t_nu  = 1e300; chi2t_sig  = 1e300
chi2t_nu1 = 1e300; chi2t_nu2  = 1e300
chi2t_sig1= 1e300; chi2t_sig2 = 1e300

prob = -1e300 if uselike else 1e300     # Initial log likeli. [small in case we need it in disc]

ini     = 0                    # TODO: meaning
inimax  = 5000                 # 
counter = 0                    # 
accrej  = np.zeros(1000)       # 
ratio   = 0.                   # 
account1= 0.                   # 

chi2tol = 50. if (pops == 1) else 90.  # more information in two tracer pops, but more errors as well
endcount = 60                  # 300 accepted models which chi2<chi2tol means initialization phase is over

rejcount = 1.                   # Rejection count
acccount = 0.                   # Acceptance count
accrejtollow  = 0.24            # Acceptance/rejection rate
accrejtolhigh = 0.26            #
farinit = 1./5. # 5 times chi2 is too high in init phase: start new from last point
stepafterrunaway = 1.5 # divide stepsize by this amount if too low fnewoverf 2.5
farover = 1./2.     # 2 times chi2 is too high after init phase 1./2.
scaleafterinit   = 1. # <= cheat: divide stepsize by this amount if init is over
stepcorr= 1.01   # adapt stepsize by this if not 0.24 < acc/rec < 0.26

# Parameters to end initphase 
initphase = True             # initialisation phase flag, first True, if over: False
endgame  = False                # Ending flag

# Units: 

# rcore in [pc], surfdens_central (=dens0) in [munit/rcore**2], and in [munit/pc**2], and totmass [munit], and max(v_LOS) in [km/s]

rcore=[]; dens0rcore=[]; dens0pc=[]; totmass=[]; maxvlos=[] # unit system
rcore_2D=[];dens0rcore_2D=[];dens0pc_2D=[]
if investigate != 'walker':
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
