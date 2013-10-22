#!/usr/bin/python

from scipy.stats.kde import gaussian_kde
from scipy.stats import norm
from numpy import linspace,hstack
from pylab import plot,draw,show,ion,ioff,hist

# creating data with two peaks
sampD1 = norm.rvs(loc=-1.0,scale=1,size=300)
sampD2 = norm.rvs(loc=2.0,scale=0.5,size=300)
samp = hstack([sampD1,sampD2])

# obtaining the pdf (my_pdf is a function!)
my_pdf = gaussian_kde(samp)

# plotting the result
x = linspace(-5,5,100)
ion()
plot(x,my_pdf(x),'r') # distribution function
hist(samp,normed=1,alpha=.3) # histogram
ioff()
show()
