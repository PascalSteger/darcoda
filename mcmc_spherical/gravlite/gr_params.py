#!/usr/bin/python2.7
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
'''global params for data analysis step 1 used to generate input for spherical MCMC'''

import numpy as np
import pdb
import gl_params as gp

from gl_class_files import *

def binparams():
  binlength=(rmax-rmin)/bins
  binmin = np.zeros(bins);  binmax = np.zeros(bins)
  rbin   = np.zeros(bins)
  for i in range(bins):
    binmin[i] = rmin+i*binlength;  binmax[i] = rmin+(i+1)*binlength
    rbin[i]   = binmin[i]+0.5*binlength
  return binmin, binmax, rbin

showplots = True

if gp.investigate == 'hernquist':
  sim = 1                                 # choose simulation
  
  # number of populations/tracers. both set to 0 means: take all particles
  pops = 2
  ntracers1 = 0
  ntracers2 = 5000
  
  if pops==1:
    totmass = 1.
    ntot =1.e6
  elif pops==2:
    if ntracers1>0: # look at dm
      totmass = 1.
      ntot = 1.e6
    else: # look at stars
      totmass = 1.e-5
      ntot = 1.e6
      
  munit = totmass/ntot
  
  rcut=1.e10                              # cutting radius
  
  bins = 12                               # set binning
  rmin = 0.; rmax = 3.
  
  
  rerror  = 0.1*(rmax-rmin)/bins         # distance error
  vrerror = 0.                            # 0.01 # velocity error. only raises sig_los!

  nit = 1000                                # set number of iterations
  procs = 16                              # number of processes for parallel execution


  dir = '/home/ast/read/dark/dwarf_data/data_hernquist/'
  #dir = '/home/psteger/sci/dwarf_data/data_hernquist/'
  
  if pops == 2:
    prename = 'dual_'
    if ntracers1 >0:
      comp = '_dm'
    else:
      comp = '_stars'
  else:
    prename = ''
    comp = ''
    
  simname = dir+'simulation/'+prename+'unit_hern_%i' %(sim)
  simpos = simname+comp+'_pos.txt'
  simvel = simname+comp+'_vel.txt' 

  fileposcartesian = dir+'simulation/'+prename+'unit_hern_%i_pos_%i_%i.txt'\
                     %(sim,ntracers1,ntracers2)
  fileposcenter = dir+'simulation/'+prename+'unit_hern_%i_centeredpos_%i_%i.txt'\
                  %(sim,ntracers1,ntracers2)
  filevelcartesian = dir+'simulation/'+prename+'unit_hern_%i_vel_%i_%i.txt' \
                     %(sim,ntracers1,ntracers2)
  
  fileposspherical = dir+'simulation/'+prename+'unit_hern_%i_sphericalpos_%i_%i.txt' \
                     %(sim,ntracers1,ntracers2)
  filevelspherical = dir+'simulation/'+prename+'unit_hern_%i_sphericalvel_%i_%i.txt' \
                     %(sim,ntracers1,ntracers2)
  filemass = dir+'enclosedmass/'+prename+'unit_hern_%i_enclosedmass_%i_%i.txt' \
             %(sim,ntracers1,ntracers2)
  
  filedenfalloff = dir+'densityfalloff/'+prename+'unit_hern_%i_falloffnotnorm_%i_%i.txt'\
                   %(sim,ntracers1,ntracers2)
  
  filesig = dir+'velocitydispersionlos/'+prename+'unit_hern_%i_veldisplos_%i_%i.txt' \
            %(sim,ntracers1,ntracers2)

  if ntracers1 == 0 and ntracers2 == 0:
    import commands
    out = commands.getstatusoutput('wc -l '+dir+'/simulation/unit_hern_1_pos.txt')
    import shlex
    my=shlex.split(out[1])
    ntracers1 = int(my[0])-1
    print 'ntracers1 = ',ntracers1
    fileposcartesian = dir+'simulation/unit_hern_%i_pos_my.txt' %(sim)
    filevelcartesian = dir+'simulation/unit_hern_%i_vel_my.txt' %(sim)
    fileposcenter    = dir+'simulation/unit_hern_%i_centeredpos.txt' %(sim)
    fileposspherical = dir+'simulation/unit_hern_%i_sphericalpos.txt' %(sim)
    filevelspherical = dir+'simulation/unit_hern_%i_sphericalvel.txt' %(sim)
    filemass         = dir+'enclosedmass/unit_hern_%i_enclosedmass.txt' %(sim)
    filedenfalloff   = dir+'densityfalloff/unit_hern_%i_falloffnotnorm.txt' %(sim)
    filesig          = dir+'velocitydispersionlos/unit_hern_%i_veldisplos.txt' %(sim)

#################### files for Walker's mock data ####################
elif gp.investigate == 'walker':
  # prior:  analysis is done out to rprior*r_core (half-light radius)
  # if < 0: include all particles, r=max(x**2+y**2)
  rprior = 3.
  r_DM  = 1000.

  rerror = 0.01
  vrerror= 0.001

  def rhodm(r):
    return rho0*(r/r_DM)**(-gamma_DM)*(1+(r/r_DM)**alpha_DM)**((gamma_DM-beta_DM)/alpha_DM)
  
  fi = Files()
  dir = fi.dir
  fil = dir+'mem2'
  
  nbins = 12
  ncomp = 4  # 3 possibilities: 0 (both), tracer pop 1, tracer pop 2, 3: foreground contamination
             # a value of 4 means: do analysis for all of them
  
  pmsplit = 0.1 # minimum probability of membership required for consideration for analysis
  fileposcartesian = dir+'simulation/pos.txt'
  filevelcartesian = dir+'simulation/vel_my.txt'

  import os; import os.path

  def newdir(bname):
    if not os.path.exists(bname):
      os.makedirs(bname)
    return

  def get_com_file(n):
    return dir+'centeredpos_%i.txt' %(n)

  def get_com_png(n):
    return dir+'centeredpos_%i.png' %(n)

  def get_params_file(n):
    return dir+'scale_%i.txt' %(n)

  def get_pos_sphere_file(n):
    return dir+'sphericalpos_%i.txt' %(n)

  def get_vel_sphere_file(n):
    return dir+'sphericalvel_%i.txt' %(n)

  def get_enc_mass_file(n):
    newdir(dir+'enclosedmass/')
    return dir+'enclosedmass/enclosedmass_%i.txt' %(n)

  def get_dens_file(n):
    newdir(dir+'nu/')
    return dir+'nu/nunotnorm_%i.txt' %(n)

  def get_dens_png(n):
    return dir+'nu/nunotnorm_%i.png' %(n)

  def get_siglos_file(n):
    newdir(dir+'siglos/')
    return dir+'siglos/siglos_%i.txt' %(n)

  def get_siglos_png(n):
    return dir+'siglos/siglos_%i.png' %(n)
