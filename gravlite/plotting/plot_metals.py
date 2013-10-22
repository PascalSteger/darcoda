#!/usr/bin/env python3

##
# @file
# plot metallicity plane for a Walker

# (c) 2013 ETHZ Pascal Steger, psteger@phys.ethz.ch

import gl_params as gp
import gr_params as gpr
from pylab import *
import numpy as np
import pdb
from gl_helper import expDtofloat

Mg,dMg,PM0,comp0=np.genfromtxt(gpr.fil,skiprows=0,unpack=True,\
                               usecols=(13,14,19,20),\
                                             dtype="d17",\
                                             converters={13:expDtofloat, # Mg0 in Angstrom\
                                                         14:expDtofloat, # delta Mg0 in Angstrom\
                                                         19:expDtofloat, # PM0 [1]\
                                                         20:expDtofloat}) # comp0 1,2,3(background)


global f, ax1


ion()
f, ax1 = plt.subplots(1,1)
draw()

pm = (PM0>0.8)


# plot components from data
pm1 = (comp0==1)*pm
pm2 = (comp0==2)*pm
pm3 = (comp0==3)*pm

plh = Mg
bins = 60
norm = False

hist(plh[pm2], bins,color='red',   alpha=0.2,  normed=norm)
hist(plh[pm1], bins,color='blue',  alpha=0.2,  normed=norm)
#hist(plh[pm3],bins,color='black', alpha=0.1)
hist(plh[pm],  bins,color='black', histtype='step', normed=norm)
xlabel('Mg')
ylabel('f')
draw()



# write metallicity file
fm = open('metals.csv', 'w')
Mgpm = Mg[pm]
for item in Mgpm:
      print(item, file=fm)
fm.close()

# call Gaussian Mixture Model, gmm, (Muratov,Gnedin 2010)
import os
os.system('GMM/gmm metals.csv 0 0.0 0.4')

# TODO: read gmm output data
fo = open('peakprob.out')


# TODO: draw population assignment

# TODO: plot histogram of GMM assigned population

# TODO: fraction of right determination?


plt.savefig('metalplane.png')
ioff();show()
