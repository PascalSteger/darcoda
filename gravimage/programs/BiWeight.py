######################################################################
# original IDL code by
# Glenn van de Ven
# Leiden Observatory, 2005, vdven@strw.leidenuniv.nl
# Princeton University, 2006, glenn@astro.princeton.edu
# Institute for Advanced Study, 2007, glenn@ias.edu
#
# python version by Silvia Garbari,
# ITP UZH Zurich silvia@physik.uzh.ch
#   May 2011
#
# extended for 4th order moment by Pascal Steger,
# ETH Zurich, pascal@steger.aero
#------------------------------------------------------------------
#
# PURPOSE:
#   Estimation of the central location ('mean') of a distribution,
#   resistant for outliers and robust for a broad range of
#   non-Gaussian underlying populations.
#
#   Estimation of the scale ('dispersion') of a distribution,
#   resistant for outliers and robust for a broad range of
#   non-Gaussian underlying populations.
#
# REFERENCES:
#   Beers, T. C., Flynn, K., & Gebhardt, K. 1990, AJ, 100, 32
#   Hoaglin, D. C., Mostelle, F., & Tukey, J. W. 1983, Understanding
#      Robust and Exploratory Data Analysis (Wiley, New York)
#
# MODIFICATION HISTORY:
#   V1.0 -- Created by Glenn van de Ven, Leiden, 12 November 2003.
#   python version by Silvia Garbari, ITP UZH Zurich silvia@physik.uzh.ch
#           May 2011
#   fourth order moment by Pascal Steger, ETH Zurich, pascal@steger.aero
#------------------------------------------------------------------


import pdb
import numpy as np
from scipy.stats import t as Tfunc
from scipy.stats import chi2 as CHI2



def stddevbiweight(vector, zero=None, eps=1e-20):
    c=9. #tuning constant
    if zero!=None:
        center=0.
    else:
        center = np.median(vector) #sample median, [vecunit]

    #median absolute deviation = deviation from the sample median
    mad = np.median(np.abs(vector-center)) #[vecunit]
    if mad<eps:
        scale = 0.
    else:
        u2 = ((vector-center)/(c*mad))**2 #[vecunit/vecunit^2]
        ind = np.argwhere(u2<=1.).flatten()
        count = len(ind) #[1]
        if count<3:
            scale=-1.
            raise Exception('stddev biweight: distribution is too strange (returning -1)')
        else:
            term1 = (1.-u2[ind])    # [1]
            term2 = (1.-5.*u2[ind]) # [1]
            n = len(vector)         # [1]
            scale = n*np.sum((vector[ind]-center)**2*term1**4) # [vecunit^2]
            scale = np.sqrt(scale)/np.sum(term1*term2) # [vecunit]

    return scale
## \fn stddevbiweight(vector, zero=None, eps=1e-20)
## PURPOSE: Estimation of the scale ('dispersion') of a distribution,
#        resistant for outliers and robust for a broad range of
# non-Gaussian underlying populations.
# -------------------------------------------------------------
# REFERENCES:
# Beers, T. C., Flynn, K., & Gebhardt, K. 1990, AJ, 100, 32
# Hoaglin, D. C., Mostelle, F., & Tukey, J. W. 1983, Understanding
# Robust and Exploratory Data Analysis (Wiley, New York)
# @param vector = distribution in vector form
# @param zero   = set this keyword to force the scale to be calculated
#                 w.r.t. zero central location, e.g. in the case the
#                 vector elements are residuals of some fit
# @param eps    = small number: zero scale if median absolute deviation MAD < EPS
# @return stddev = scale of the distribution, -1 in case of failure



def meanbiweight(vector,itmax=10,fracmin=[],eps=1e-24,\
                 ci_perc=95,ci_mean=True,ci_std=True):
    c = 6. #tuning constant
    n = len(vector) #[1]

    # fractional change in MAD
    if fracmin==[]:
        fracmin = 0.03*np.sqrt(0.5/(n-1)) #[1]

    #calculate scale
    STDDEV = stddevbiweight(vector) #[vecunit]

    #take median as first guess of the central location
    center = np.median(vector) #[vecunit]
    #calculate MAD
    mad = np.median(np.abs(vector-center)) #[vecunit]

    #if mad is already 0, no iteration is needed
    if mad < eps: #[vecunit], [1] !! attention, eps is though of in vecunits too
        #zero weight for elements >3 stddev away from median
        limit = 3.*STDDEV #[vecunit]
        wei = np.zeros(len(vector)) #[1]
        i3std = np.argwhere(np.abs(vector-center)<limit).flatten() #[1]
        wei[i3std] = 1./n #[1]
    else:
        frac = 1e30 #[1]
        it = 0
        while frac > fracmin and it < itmax:
            it += 1
            u2 = ((vector-center)/(c*mad))**2 #[(vecunit/vecunit)**2 = 1]
            u2[np.argwhere(u2>1.).flatten()]=1 # replace all values >1 with 1 #[1]
            wei = (1.-u2)**2 # zero weight if |u|=1 (u2=1) #[1]
            wei /= np.sum(wei) #normalize #[1]
            center = np.sum(wei*vector) #new estimate of the center location #[vecunit]
            mad_previous = mad #save previous value of mad #[vecunit]
            mad=np.median(np.abs(vector-center)) #new value of mad
            if mad < eps: #[vecunit, 1] ! attention, eps is thought of in vecunits!
                frac = 0. #no further iteration is necessary
            else:
                frac = np.abs(mad_previous-mad)/mad_previous #[1]

    ci_alpha = (1-ci_perc*1e-2) #[1]

    if ci_mean:
        ci_len = Tfunc.ppf(ci_alpha*0.5,0.7*(n-1))*STDDEV/np.sqrt(n) #[vecunit]
        ci_mean = np.zeros(2)
        ci_mean[0] = center - ci_len #[vecunit]
        ci_mean[1] = center + ci_len #[vecunit]
    else:
        ci_mean = np.zeros(2) #[vecunit]

    if ci_std:
        ci_std = np.empty(2)
        ci_std[0] = -STDDEV*(np.sqrt((n-1)/CHI2.ppf(1-ci_alpha*0.5,n-1))-1) #[vecunit]
        ci_std[1] =  STDDEV*(np.sqrt((n-1)/CHI2.ppf(ci_alpha*0.5,n-1))-1)   #[vecunit]
    else:
        ci_std = np.zeros(2)

    if type(ci_mean)!=np.ndarray or type(ci_std)!=np.ndarray:
        print('out=',center,STDDEV,ci_mean,ci_std)
        ci_mean = np.zeros(2)
        ci_std  = np.zeros(2)
    #[vecunit],    [vecunit] [vecunit]   [vecunit]   [vecunit]  [vecunit]
    return center, STDDEV,   ci_mean[0], ci_mean[1], ci_std[0], ci_std[1]
## \fn meanbiweight(vector,itmax=10,fracmin=[],eps=1e-24, ci_perc=95,ci_mean=True,ci_std=True):
# Estimation of the central location ('mean') of a distribution,
# resistant for outliers and robust for a broad range of
# non-Gaussian underlying populations.
# REFERENCES:
# Beers, T. C., Flynn, K., & Gebhardt, K. 1990, AJ, 100, 32
# Hoaglin, D. C., Mostelle, F., & Tukey, J. W. 1983, Understanding
# Robust and Exploratory Data Analysis (Wiley, New York)
# @param vector      = distribution in vector form
# @param itmax =10    = maximum number of iterations used in determination of
#                      the central location; default 10
# @param fracmin =[]   = minimum fractional change in median absolute deviation
#                      (MAD), used as convergence criterion in iteration;
#                      default 0.03*sqrt(0.5/(n-1)), with n number of vector
#                      elements
# @param eps =1e-24   = small number: iteration is stopped if MAD < EPS;
#                      default 1.0e-24
# @param ci_perc =95  = default confidence interval (95), 1sigma=68.4
# @param ci_mean =True -> if True (default) the code output the confidence
#                        interval for the mean as output[2] and output[3],
#                        else ouput[2]=output[3]=0
# @param ci_std =True -> if True (default) the code output the confidence
#                       interval for the standard deviation as output[4]
#                       and output[5], else output[4]=output[5]=0
# @return mean      = central location of the distribution
#         STDDEV    = scale of the distribution
#         ci_mean[0],ci_mean[1]= confidence interval for the mean (default 95%)
#         ci_std[0],ci_std[1]= confidence interval for the mean (default 95%)
