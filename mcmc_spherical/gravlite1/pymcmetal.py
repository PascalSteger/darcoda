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



def bimodal_gauss(data,pm):
    '''run MCMC to get regression on bimodal normal distribution'''
    m1 = np.mean(data[pm])/2.
    m2 = np.mean(data[pm])*2.
    dm = m2 - m1
    size = len(data[pm])

    ### set up model
    p = Uniform( "p", 0.2 , 0.8) #this is the fraction that come from mean1 vs mean2
    # p = distributions.truncated_normal_like('p', mu=0.5, tau=0.001, a=0., b=1.)
    # p = Normal( 'p', mu=(1.*sum(comp0==1))/size, tau=1./0.1**2 ) # attention: wings!, tau = 1/sig^2
    # p = Normal( 'p', mu=0.5, tau=1./0.1**2 ) # attention: wings!, tau = 1/sig^2
    
    ber = Bernoulli( "ber", p = p, size = size) # produces 1 with proportion p
    precision = Gamma('precision', alpha=0.01, beta=0.01)
    
    dmu = Normal( 'dmu', dm, tau=1./0.05**2 ) # [PS] give difference between means, finite
    # dmu = Lognormal( 'dmu', 0.3, tau=1./0.1**2)
    
    mean1 = Normal( "mean1", mu = m1,          tau = 1./0.1**2 ) # better to use Normals versus Uniforms,
                                                                 # if not truncated
    mean2 = Normal( "mean2", mu = mean1 + dmu, tau = 1./0.1**2 ) # tau is 1/sig^2
    
    @deterministic
    def mean( ber = ber, mean1 = mean1, mean2 = mean2):
        return ber*mean1 + (1-ber)*mean2

    
    obs = Normal( "obs", mean, precision, value = data[pm], observed = True)
    model = Model( {"p":p, "precision": precision, "mean1": mean1, "mean2":mean2, "obs":obs} )
    
    from pymc import MCMC, Matplot



    M = MCMC(locals(), db='pickle', dbname='metals.pickle')
    iter = 3000; burn = 2000; thin = 10
    M.sample(iter=iter, burn=burn, thin=thin)
    M.db.commit()

    mu1 = np.mean(M.trace('mean1')[:])
    mu2 = np.mean(M.trace('mean2')[:])
    p   = np.mean(M.trace('p')[:])
    return p, mu1, 0.1, mu2, 0.1







def assign_pop(data, pm, p, mu1, sig1, mu2, sig2):
    from scipy.stats import norm
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









    ### plotting
    # ion()
    # h1 = hist(data,50,color='k',alpha=0.0,normed=True)
    # print 'overall mean = ', np.mean(data)
    
    # hist(data[comp0==1],30,color='r',alpha=0.5,normed=True)
    # axvline(x=np.mean(data[comp0==1]),color='r',ls='dashed')
    # print 'mean 1 = ', np.mean(data[comp0==1])

    # hist(data[comp0==2],30,color='b',alpha=0.5,normed=True)
    # axvline(x=np.mean(data[comp0==2]),color='b',ls='dashed')
    # print 'mean 2 = ', np.mean(data[comp0==2])


    # from scipy.stats import norm
    # x = np.linspace(-1., 1.5, 100)
    # y = (p*norm.pdf(x, loc=mu1, scale=0.1) + (1-p)*norm.pdf(x,loc=mu2, scale=0.1))
    # plot(x, y/max(y)*max(h1[0]),color='k')

    # ioff()
    # savefig('data.png')

    # ion()
    # Matplot.plot(M)

    # ioff()
    # show()
