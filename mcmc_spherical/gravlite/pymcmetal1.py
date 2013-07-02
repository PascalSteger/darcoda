#!/usr/bin/python
import numpy as np
import sys
import pdb
import pymc as mc

from pylab import *
import gl_params as gp
import gr_params as gpr
from gl_helper import expDtofloat
from gl_class_files import *
from scipy.stats import norm


def bimodal_gauss(data):
    '''run MCMC to get regression on bimodal normal distribution'''
    # http://nbviewer.ipython.org/urls/raw.github.com/CamDavidsonPilon/Probabilistic-Programming-and-Bayesian-Methods-for-Hackers/master/Chapter3_MCMC/IntroMCMC.ipynb

    p = mc.Uniform( "p", 0.3, 0.7)
    assignment = mc.Categorical("assignment", [p, 1-p], size = data.shape[0] ) 
    print "prior assignment, with p = %.2f:"%p.value
    print assignment.value[:100], "..."

    taus = 1.0/mc.Uniform( "stds", 0., 0.2, size=2)**2
    centers = mc.Normal( "centers", [0.0, 0.3], [1/0.1**2, 1/0.1**2], size=2)
    pdb.set_trace()

    @mc.deterministic 
    def center_i( assignment = assignment, centers = centers ):
        return centers[ assignment] 
    
    @mc.deterministic
    def tau_i( assignment = assignment, taus = taus ):
        return taus[ assignment] 

    print "Random assignments: ", assignment.value[ :10 ], "..."
    print "Assigned center: ", center_i.value[ :10], "..."
    print "Assigned precision: ", tau_i.value[ :10], "..."


    # and to combine it with the observations:
    observations = mc.Normal( "obs", center_i, tau_i, value = data, observed = True )

    # below we create a model class
    model = mc.Model( [p, assignment, taus, centers, observations ] )

    M = mc.MCMC( model )#, db='pickle', dbname='metals.pickle')
    iter = 50000; burn = 40000; thin = 10
    M.sample( iter=iter, burn = burn, thin = thin)

    # M.db.commit()
    pdb.set_trace()
    mu1, mu2 = M.trace('centers')[:].mean(axis=0) 
    sig1, sig2 = M.trace('stds')[:].mean(axis=0)
    p = M.trace("p")[:].mean()
    

    return p, mu1, sig1, mu2, sig2, M


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

    data = np.array(Mg0[pm])
    # data = np.hstack([data, data+0.4])
    # import random
    # random.shuffle(data)
    p, mu1, sig1, mu2, sig2, M = bimodal_gauss(data)


    # for pretty colors later in the book.
    colors = ["#348ABD", "#A60628"]# if center_trace[-1,0] > center_trace[-1,1] else [ "#A60628", "#348ABD"]

    ion()
    bins = np.linspace(min(data),max(data),100)
    h1 = hist(data,bins=bins,color='k', alpha=0.2, normed=True)
    print 'overall mean = ', np.mean(data)
    
#    hist(data[comp0==1],bins=bins,color='r',alpha=0.5,normed=True)
#    axvline(x=np.mean(data[comp0==1]),color='r',ls='dashed')
#    print 'mean 1 = ', np.mean(data[comp0==1]), ', found ',mu1,', ',sig1

#    hist(data[comp0==2],bins=bins,color='b',alpha=0.5,normed=True)
#    axvline(x=np.mean(data[comp0==2]),color='b',ls='dashed')
#    print 'mean 2 = ', np.mean(data[comp0==2]), ', found ',mu2,', ',sig2


    x = bins
    y = p * norm.pdf(x, loc = mu1, scale = sig1)
    plt.plot(x, y, label = "posterior pop 1", lw = 2 )
    plt.fill_between( x, y, color = colors[0], alpha = 0.3 )
    
    y = (1-p) * norm.pdf(x, loc = mu2, scale = sig2)
    plt.plot(x, y, label = "posterior pop 2", lw = 2 )
    plt.fill_between( x, y, color = colors[1], alpha = 0.3 )

    y = p * norm.pdf(x, loc = mu1, scale = sig1) + (1-p) * norm.pdf(x, loc = mu2, scale = sig2)
    plt.plot(x, y, label = "overall fitted distribution", lw = 2 )
    
    plt.legend(loc = "upper left")
    plt.title( "Visualizing Clusters using posterior-mean parameters" );

    ioff()
    savefig('data.png')

    # ion()
    # mc.Matplot.plot(M)
    # ioff()

    # ion()
    # cmap = mpl.colors.LinearSegmentedColormap.from_list("BMH", colors )
    # assign_trace = M.trace("assignment")[:]
    # scatter( data,1-assign_trace.mean(axis=0), cmap=cmap, c=assign_trace.mean(axis=0), s = 50)
    # # plt.ylim( -0.05, 1.05 )
    # # plt.xlim( 35, 300)
    # plt.title( "Probability of data point belonging to cluster 0" )
    # plt.ylabel("probability")
    # plt.xlabel("value of data point" );
    # ioff()
    show()
