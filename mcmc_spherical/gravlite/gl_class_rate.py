#!/usr/bin/python2.7
# class to get a rolling acceptance/rejection rate


import collections
import pdb
import gl_params as gp
if gp.geom == 'sphere':
    import physics_sphere as phys
elif gp.geom == 'disc':
    import physics_disc as phys
import gl_helper as gh


    
class Rate:
    'Common base class for all filename sets'
    
    def __init__ (self, aimrate, eps):
        self.aimrate = aimrate
        self.eps     = eps
        self.lowaim, self.highaim = aimrate-eps, aimrate+eps
        self.values = collections.deque(maxlen=gp.rollsize)
        self.oldrate = aimrate
        
        from random import shuffle
        x = [i<gp.rollsize*aimrate for i in range(gp.rollsize)]
        shuffle(x)
        
        for i in range(gp.rollsize):
            self.values.append(x[i])
            
        return

    def rate(self):
        return sum(self.values)/(1.*gp.rollsize)

    def update(self, boo):
        self.oldrate = self.rate()
        self.values.append(boo)

    def rightrate(self):
        ar = self.rate()
        if ar >= self.aimrate-self.eps and ar <= self.aimrate+self.eps:
            return True
        else:
            return False
        
    def getsbetter(self):
        # True if diff to aimrate is smaller than before
        if abs(self.rate()-self.aimrate) < abs(self.oldrate - self.aimrate):
            return True
        else:
            return False
