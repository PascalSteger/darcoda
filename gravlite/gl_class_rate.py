#!/usr/bin/env python3

##
# @file
# class to get a rolling acceptance/rejection rate

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import collections
import pdb
import gl_params as gp
if gp.geom == 'sphere':
    import physics_sphere as phys
elif gp.geom == 'disc':
    import physics_disc as phys
import gl_helper as gh


## Common base class for all filename sets
class Rate:

    ## constructor
    # @param aimrate rate of acceptance/total runs to aim for
    # @param eps allowed range around aimrate
    def __init__ (self, aimrate, eps):
        ## rate of acceptance/total runs to aim for
        self.aimrate = aimrate
        ## allowed range around aimrate
        self.eps     = eps
        ## low and high boundaries for rate
        self.lowaim, self.highaim = aimrate-eps, aimrate+eps
        ## buffer for boolean values for accepted/rejected models
        self.values = collections.deque(maxlen=gp.rollsize)
        ## old aimrate, for historical comparison
        self.oldrate = aimrate
        
        from random import shuffle
        x = [i<gp.rollsize*aimrate for i in range(gp.rollsize)]
        shuffle(x)
        
        for i in range(gp.rollsize):
            self.values.append(x[i])
            
        return

    ## calculate rate from values stored in class
    def rate(self):
        return sum(self.values)/(1.*gp.rollsize)

    ## update rate, store old ones
    def update(self, boo):
        self.oldrate = self.rate()
        self.values.append(boo)

    ## determine whether rate of acceptance lies within aimrate +/- eps
    def rightrate(self):
        ar = self.rate()
        if ar >= self.aimrate-self.eps and ar <= self.aimrate+self.eps:
            return True
        else:
            return False

    ## return True if diff to aimrate is smaller than before
    def getsbetter(self):
        if abs(self.rate()-self.aimrate) < abs(self.oldrate - self.aimrate):
            return True
        else:
            return False
