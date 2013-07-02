#!/usr/bin/python
import numpy as np
import sys
import pdb
from pymc import *

from pylab import *
import gl_params as gp
import gr_params as gpr
from gl_helper import expDtofloat
from gl_class_files import *
from scipy.stats import norm


def bimodal_gauss(data,pm,dmin=0.3):
    '''run MCMC to get regression on bimodal normal distribution'''
    size = len(data[pm])

    

    ### set up model
    p = Uniform( "p", 0.2 , 0.8) #this is the fraction that come from mean1 vs mean2
    # p = distributions.truncated_normal_like('p', mu=0.5, tau=0.001, a=0., b=1.)
    # p = Normal( 'p', mu=(1.*sum(comp0==1))/size, tau=1./0.1**2 ) # attention: wings!, tau = 1/sig^2
    # p = Normal( 'p', mu=0.5, tau=1./0.1**2 ) # attention: wings!, tau = 1/sig^2
    
    ber = Bernoulli( "ber", p = p, size = size) # produces 1 with proportion p
    precision = Gamma('precision', alpha=0.01, beta=0.01)
    
    mean1 = Uniform( "mean1", -0.5, 1.0) # if not truncated
    sig1  = Uniform( 'sig1',  0.01, 1.)
    mean2 = Uniform( "mean2", mean1 + dmin, 1.5)
    sig2  = Uniform( 'sig2',  0.01, 1.)

    pop1  = Normal( 'pop1', mean1, 1./sig1**2) # tau is 1/sig^2
    pop2  = Normal( 'pop2', mean2, 1./sig2**2)


    @deterministic
    def bimod(ber = ber, pop1 = pop1, pop2 = pop2): # value determined from parents completely
        return ber*pop1 + (1-ber)*pop2

    obs = Normal( "obs", bimod, precision, value = data[pm], observed = True)
    model = Model( {"p":p, "precision": precision, "mean1": mean1, 'sig1': sig1, "mean2":mean2, 'sig2':sig2, "obs":obs} )
    
    from pymc import MCMC, Matplot


    M = MCMC(locals(), db='pickle', dbname='metals.pickle')
    iter = 10000; burn = 9000; thin = 10
    M.sample(iter=iter, burn=burn, thin=thin)
    M.db.commit()

    mu1 = np.mean(M.trace('mean1')[:])
    sig1= np.mean(M.trace('sig1')[:])
    mu2 = np.mean(M.trace('mean2')[:])
    sig2= np.mean(M.trace('sig2')[:])
    p   = np.mean(M.trace('p')[:])
    return p, mu1, sig1, mu2, sig2, M







def assign_pop(data, pm, p, mu1, sig1, mu2, sig2):
    size = len(data)
    pm1 = []; pm2 = []
    for i in range(size):
        x = 1.*data[i]
        if np.random.rand() < p*norm.pdf(x, loc=mu1, scale=sig1)/\
           (p*norm.pdf(x,loc=mu1,scale=sig1)+(1.-p)*norm.pdf(x,loc=mu2,scale=sig2)):
            pm1.append(True); pm2.append(False)
        else:
            pm1.append(False);pm2.append(True)
    return np.array(pm1)*pm, np.array(pm2)*pm # array needed to use it as subset indicator


def center(x,y):
    x = x-np.mean(x)
    y = y-np.mean(y)
    return x,y

def rhalf(x,y):
    r0 = np.sqrt(x**2+y**2)
    r0 = np.sort(r0)
    return r0[len(r0)/2]



if __name__=="__main__":
    x0,y0,z0,vz0,vb0,Mg0,PM0,comp0=np.genfromtxt(gpr.fil,skiprows=0,unpack=True,\
                                                 usecols=(0,1,2,11,12,13,19,20),\
                                                 dtype="d17",\
                                                 converters={0:expDtofloat,  # x0  in pc \
                                                             1:expDtofloat,  # y0  in pc \
                                                             2:expDtofloat,  # z0  in pc \
                                                             11:expDtofloat, # vz0 in km/s\
                                                             12:expDtofloat, # vb0(LOS due binary), km/s\
                                                             13:expDtofloat, # Mg0 in Angstrom\
                                                             19:expDtofloat, # PM0 [1]\
                                                             20:expDtofloat}) # comp0 1,2,3(background)

    pm = (PM0 >= gpr.pmsplit) * (comp0<3) # exclude outliers .  TODO: do that without cheating
    pm1 = pm*(comp0==1)                   # are updated below if gp.metalpop is set to True
    pm2 = pm*(comp0==2)                   # same same
    pm3 = pm*(comp0==3)

    data = Mg0
    p,mu1,sig1,mu2,sig2,M = bimodal_gauss(data,pm,0.3) # only take Mg0[pm] for the task
    pm1new, pm2new = assign_pop(Mg0,pm,p,mu1,sig1,mu2,sig2)

    ion()
    bins = np.linspace(min(data[pm]),max(data[pm]),60)
    h1 = hist(data[pm],bins=bins,color='k',alpha=0.2,normed=False)
    print 'overall mean = ', np.mean(data)
    
    hist(data[comp0==1],bins=bins,color='r',alpha=0.5,normed=False)
    axvline(x=np.mean(data[comp0==1]),color='r',ls='dashed')
    print 'mean 1 = ', np.mean(data[comp0==1])

    hist(data[comp0==2],bins=bins,color='b',alpha=0.5,normed=False)
    axvline(x=np.mean(data[comp0==2]),color='b',ls='dashed')
    print 'mean 2 = ', np.mean(data[comp0==2])


    x = np.linspace(-1., 1.5, 100)
    y1 = p*norm.pdf(x, loc=mu1, scale=0.1)
    y2 = (1-p)*norm.pdf(x,loc=mu2, scale=0.1)
    y = y1 + y2
    
    # plot(x, y1/max(y1)*max(h1[0]),color='k')
    plot(x, y1/max(y)*max(h1[0]),color='b')
    plot(x, y2/max(y)*max(h1[0]),color='r')
    plot(x, y/max(y)*max(h1[0]),color='k')

    ioff()
    savefig('data.png')

    ion()
    Matplot.plot(M)
    ioff()
    show()
