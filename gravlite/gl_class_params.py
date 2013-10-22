#!/usr/bin/env python

##
# @file
# parameters for the MCMC
# (c) 2013 ETHZ Pascal S.P. Steger

import pdb
import numpy as np
import numpy.random as npr
import random
import gl_params as gp
if gp.geom == 'sphere':
    import physics_sphere as phys
else:
    import physics_disc as phys
import gl_helper as gh




## Common base class for all parameter sets
class Params:

    ## constructor, with modes depending on locpop
    # @param locpop mode
    # @param nu1  preset tracer density falloff
    # @param dens preset overall density falloff
    # @param delta1 preset velocity anisotropy for the first component
    # @param norm1 preset norm factor for disc
    # @param nu2 preset tracer density falloff for the second component
    # @param delta2 preset velocity anisotropy for the second component
    # @param norm2 preset norm factor for disc
    def __init__ (self, locpop,nu1=None,dens=None,delta1=None,norm1=None,nu2=None,delta2=None,norm2=None):
        if locpop == -2: # set according to values
            ## tracer density of population 1
            self.nu1    = np.zeros(gp.nipol)+nu1
            ## tracer density of population 2
            self.nu2    = np.zeros(gp.nipol)+nu2
            ## overall mass density
            self.dens    = np.zeros(gp.nipol)+dens
            ## density of velocity anisotropy 1
            self.delta1  = np.zeros(gp.nipol)+delta1
            ## density of velocity anisotropy 2
            self.delta2  = np.zeros(gp.nipol)+delta2
            ## integration constant for disc case, pop 1
            self.norm1   = np.zeros(1)+norm1
            ## integration constant for disc case, pop 2
            self.norm2   = np.zeros(1)+norm2
        if locpop == -1: # set according to values
            self.nu1    = np.zeros(gp.nipol)+nu1
            self.dens    = np.zeros(gp.nipol)+dens
            self.delta1  = np.zeros(gp.nipol)+delta1
            self.norm1   = np.zeros(1)+norm1
        elif locpop==0: # set all to -1 or 0
            self.nu1    = np.zeros(gp.nipol)
            self.nu2    = np.zeros(gp.nipol)
            self.dens    = np.zeros(gp.nipol)
            self.delta1  = np.zeros(gp.nipol)
            self.delta2  = np.zeros(gp.nipol)
            self.norm1   = np.zeros(1)
            self.norm2   = np.zeros(1)
        elif locpop == 1: # fill with values for one population only
            self.nu1    = nu1
            self.dens   = dens
            self.delta1 = delta1
            self.norm1   = np.zeros(1)+norm1
            self.norm2   = np.zeros(1)
        elif locpop==2:   # fill with values for two populations
            self.nu1    = nu1
            self.nu2    = nu2
            self.dens    = dens
            self.delta1  = delta1
            self.delta2  = delta2
            self.norm1   = np.zeros(1)+norm1
            self.norm2   = np.zeros(1)+norm2

    ## set all profiles to a random value
    def setuniformrandom(self):
        self.nu1    = npr.uniform(-1,1,gp.nipol)
        self.dens    = npr.uniform(-1,1,gp.nipol)
        self.delta1  = npr.uniform(-1,1,gp.nipol)
        self.norm1   = npr.uniform(-1,1,1)
        if gp.pops==2:
            self.nu2    = npr.uniform(-1,1,gp.nipol)
            self.delta2  = npr.uniform(-1,1,gp.nipol)
            self.norm2   = npr.uniform(-1,1,1)

        return self

    ## only set tracer density of first component to random values
    def wiggle_nu1(self):
        self.nu1 = npr.uniform(-1, 1, gp.nipol)
        return self

    ## only set tracer density of second component to random values
    def wiggle_nu2(self):
        self.nu2 = npr.uniform(-1, 1, gp.nipol)
        return self

    ## change both tracer densities to random values
    def wiggle_nu(self):
        self.wiggle_nu1()
        if gp.pops==2:
            self.wiggle_nu2()
        return self

    ## set the overall density falloff to random values
    def wiggle_dens(self):
        if gp.constdens:
            self.dens = np.zeros(gp.nipol) + npr.uniform(-1,1)
        else:
            self.dens = npr.uniform(-1, 1, gp.nipol)
        return self

    ## set all velocity anisotropies to random values
    def wiggle_delta(self):
        self.delta1 = npr.uniform(-1,1,gp.nipol)
        if gp.pops ==2:
            self.delta2 = npr.uniform(-1,1,gp.nipol)
        return self

    ## Wiggle tilt parameters [if needed]:
    def wiggle_tilt(self):
        # rant = 2.0 * (0.5 - randomu(seed,/DOUBLE,ntparsR))
        # tparsRt = tparsR + rant * tparsRerr

        ## new variable tilt in disc case. TODO: not needed?
        self.tilt = npr.uniform(-1,1,gp.nipol)
        return self

    ## set all profiles at one given radial position to random values
    def wigglebin(self,p=0):
        self.nu1[p]   = npr.uniform(-1,1)
        self.dens[p]   = npr.uniform(-1,1)
        self.delta1[p] = npr.uniform(-1,1)
        if gp.pops==2:
            self.nu2[p]   = npr.uniform(-1,1)
            self.delta2[p] = npr.uniform(-1,1)
        return self

    ## set both norms to random values
    def wigglenorm(self):
        self.norm1   = npr.uniform(-1,1)
        if gp.pops==2:
            self.norm2   = npr.uniform(-1,1)
        return self

    ## dot product of all profiles with a constant factor nr
    # @param nr constant factor
    def dot(self,nr):
        self.nu1    *= nr;        self.dens    *= nr;        self.delta1  *= nr
        self.norm1    *= nr
        if gp.pops == 2:
            self.nu2    *= nr;    self.delta2  *= nr;        self.norm2   *= nr
        return self

    ## dot division of all profiles with a constant factor nr
    # @param nr constant factor
    def divdot(self,nr):
        self.nu1   /= nr;         self.dens   /= nr;         self.delta1 /= nr
        self.norm1   /= nr
        if gp.pops == 2:
            self.nu2   /= nr;     self.delta2 /= nr;         self.norm2 /= nr
        return self

    ## multiply all profiles with profiles from another class instance
    # @param pars class of parameters
    def mul(self, pars):
        self.nu1   *= pars.nu1;   self.dens   *= pars.dens;   self.delta1 *= pars.delta1
        self.norm1   *= pars.norm1
        if gp.pops == 2:
            self.nu2 *= pars.nu2; self.delta2 *= pars.delta2; self.norm2 *= pars.norm2
        return self

    ## add another class instance's profiles to our profiles
    # @param pars class of parameters
    def add(self, pars):
        self.nu1    += pars.nu1;  self.dens    += pars.dens;  self.delta1  += pars.delta1
        self.norm1    += pars.norm1
        if gp.pops==2:
            self.nu2 += pars.nu2; self.delta2 += pars.delta2; self.norm2 += pars.norm2
        return self

    ## take 10-base logarithm of all sensible profiles
    def getlog(self):
        if gp.pops==1:
            return Params(gp.pops,\
                          np.log10(self.nu1),\
                          np.log10(self.dens),\
                          self.delta1,\
                          self.norm1)

        elif gp.pops==2:
            return Params(gp.pops,\
                          np.log10(self.nu1),\
                          np.log10(self.dens),\
                          self.delta1,\
                          self.norm1,\
                          np.log10(self.nu2),\
                          self.delta2,\
                          self.norm2)

    ## take over values from another instance
    # @param other another class of parameters
    def assign(self,other):
        self.nu1    = other.nu1
        self.dens   = other.dens
        self.delta1 = other.delta1
        self.norm1   = other.norm1
        if gp.pops==2:
            self.nu2    = other.nu2
            self.delta2  = other.delta2
            self.norm2   = other.norm2
        return self

    ## short output in pretty format: rounded to 3 digits
    def output(self):
        print('  # nu1     = ',gh.pretty(self.nu1))
        print('  # dens     = ',gh.pretty(self.dens))
        print('  # delta1   = ',gh.pretty(self.delta1))
        print('  # norm1    = ',gp.pretty(self.norm1))
        if gp.pops==2:
            print('  # nu2     = ',gh.pretty(self.nu2))
            print('  # delta2   = ',gh.pretty(self.delta2))
            print('  # norm2    = ',gh.pretty(self.norm2))
        return
    

    ## multiply parameter set with worst contribution to chi2 with a constant factor,
    # to increase the stepsize of all profiles going into the calculation of said parameter set
    # @param mult constant factor
    def adaptworst(self,mult):
        mult *= mult # go up twice, decrease for all once
        if gp.pops == 1:
            if gp.chi2t_nu1 > gp.chi2t_sig1:
                self.nu1 *= mult
                self.norm1 *= mult
            else:
                # change one of the possible parameter steps
                if npr.rand() < 0.5:
                    self.nu1  *= mult
                    self.norm1 *= mult
                else:
                    self.dens *= mult
                self.delta1 *= mult

        if gp.pops==2:
            # here, we adapt stepsize first for worse nu, then for worse sigma
            # first adapt worse nu (s.t. not only one of (nu,delta) set are changed)
            if gp.chi2t_nu > gp.chi2t_sig:
                if gp.chi2t_nu1 > gp.chi2t_nu2:
                    self.nu1 *= mult
                    self.norm1 *= mult                
                else:
                    self.nu2 *= mult
                    self.norm2 *= mult
            else:
                # then adapt worse set
                if gp.chi2t1 > gp.chi2t2:
                    # TODO: do it only if sigma error > nu error
                    # if gp.chi2t_sig1 > gp.chi2t_nu1:
                    self.nu1    *= mult
                    self.dens   *= mult
                    self.delta1 *= mult
                else:
                    # if gp.chi2t_sig2 > gp.chi2t_nu2:
                    self.nu2    *= mult
                    self.dens   *= mult
                    self.delta2 *= mult

            # TODO: attention: has more emphasis on increasing nu now!

    ## multiply all profiles with a constant factor
    def adaptall(self,mult):
        self.nu1 *= mult
        self.dens *= mult
        self.norm1 *= mult
        self.delta1 *= mult
        if gp.pops == 2:
            self.nu2 *= mult
            self.norm2 *= mult
            self.delta2 *= mult
        return

    ## scale self depending on which chi2 is higher
    def scale_prop_chi2(self):
        if gp.pops == 1:
            # depending on nu parametrization: only
            self.nu1    *= gp.chi2t_nu1/gp.chi2t

            # depending on siglos, to be changed if mainly siglos contributing
            self.delta1 *= gp.chi2t_sig1/gp.chi2t
            self.dens   *= gp.chi2t_sig1/gp.chi2t
            
        elif gp.pops == 2:
            # depending on nu parametrization:
            self.nu1    *= np.sqrt(gp.chi2t_nu1/gp.chi2t)
            self.nu2    *= np.sqrt(gp.chi2t_nu2/gp.chi2t)


            # depending on siglos
            self.delta1 *= np.sqrt(gp.chi2t_sig1/gp.chi2t)
            self.delta2 *= np.sqrt(gp.chi2t_sig2/gp.chi2t)
                
            self.dens   *= np.sqrt(gp.chi2t_sig/gp.chi2t)
        return self

    ## determine whether tracer densities go negative somewhere
    def has_negative(self):
        if min(self.nu1) < 0. and not gp.nulog:
            return True
        if gp.geom=='sphere' and min(phys.dens(gp.xipol,self.dens)<0.):
            return True
        if gp.pops==2:
            if min(self.nu2) < 0. and not gp.nulog:
                return True
        return False

    ## return back all tracer densities
    def get_nu(self):
        if gp.pops==1:
            return [self.nu1]
        elif gp.pops==2:
            return [self.nu1, self.nu2]

    ## set all tracer densities
    def set_nu(self,nus):
        if gp.pops==1:
            self.nu1 = nus[0]
        elif gp.pops==2:
            self.nu1 = nus[0]
            self.nu2 = nus[1]
        return self

    ## get all velocity anisotropies
    def get_delta(self):
        if gp.pops==1:
            return [self.delta1]
        elif gp.pops==2:
            return [self.delta1, self.delta2]

    ## set all velocity anisotropies
    # @param deltas given velocity anisotropies
    def set_delta(self,deltas):
        if gp.pops==1:
            self.delta1 = deltas[0]
        elif gp.pops ==2:
            self.delta1 = deltas[0]
            self.delta2 = deltas[1]
        return self
