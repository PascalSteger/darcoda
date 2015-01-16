#!/usr/bin/env ipython3

## \file

# test inverse of complimentary error function in double precision

import numpy as np

def dierfc_alt(y):
    z = y
    if (y > 1):
        z = 2 - y
    w = 0.916461398268964 - np.log(z)
    u = np.sqrt(w)
    s = (np.log(u) + 0.488826640273108) / w
    t = 1 / (u + 0.231729200323405)
    x = u * (1 - s * (s * 0.124610454613712 + 0.5)) -\
        ((((-0.0728846765585675 * t + 0.269999308670029) * t +\
        0.150689047360223) * t + 0.116065025341614) * t +\
        0.499999303439796) * t
    t = 3.97886080735226 / (x + 3.97886080735226)
    u = t - 0.5
    s = (((((((((0.00112648096188977922 * u +\
        1.05739299623423047e-4) * u - 0.00351287146129100025) * u -\
        7.71708358954120939e-4) * u + 0.00685649426074558612) * u +\
        0.00339721910367775861) * u - 0.011274916933250487) * u -\
        0.0118598117047771104) * u + 0.0142961988697898018) * u +\
        0.0346494207789099922) * u + 0.00220995927012179067
    s = ((((((((((((s * u - 0.0743424357241784861) * u -\
        0.105872177941595488) * u + 0.0147297938331485121) * u +\
        0.316847638520135944) * u + 0.713657635868730364) * u +\
        1.05375024970847138) * u + 1.21448730779995237) * u +\
        1.16374581931560831) * u + 0.956464974744799006) * u +\
        0.686265948274097816) * u + 0.434397492331430115) * u +
        0.244044510593190935) * t -\
        z * np.exp(x * x - 0.120782237635245222)
    x += s * (x * s + 1)
    if (y > 1):
        x = -x
    return x
## \fn dierfc_alt(y)
# alternative inverse of complimentary error function
# @param y

def dierfc(y):
    y = np.array(y)
    qa=9.16461398268964e-01
    qb=2.31729200323405e-01
    qc=4.88826640273108e-01
    qd=1.24610454613712e-01
    q0=4.99999303439796e-01
    q1=1.16065025341614e-01
    q2=1.50689047360223e-01
    q3=2.69999308670029e-01
    q4=-7.28846765585675e-02
    pa=3.97886080735226000e+00
    pb=1.20782237635245222e-01
    p0=2.44044510593190935e-01
    p1=4.34397492331430115e-01
    p2=6.86265948274097816e-01
    p3=9.56464974744799006e-01
    p4=1.16374581931560831e+00
    p5=1.21448730779995237e+00
    p6=1.05375024970847138e+00
    p7=7.13657635868730364e-01
    p8=3.16847638520135944e-01
    p9=1.47297938331485121e-02
    p10=-1.05872177941595488e-01
    p11=-7.43424357241784861e-02
    p12=2.20995927012179067e-03
    p13=3.46494207789099922e-02
    p14=1.42961988697898018e-02
    p15=-1.18598117047771104e-02
    p16=-1.12749169332504870e-02
    p17=3.39721910367775861e-03
    p18=6.85649426074558612e-03
    p19=-7.71708358954120939e-04
    p20=-3.51287146129100025e-03
    p21=1.05739299623423047e-04
    p22=1.12648096188977922e-03

    # remove insensible ranges
    y[y==0] = 1e-15
    y[y==1] = 1.-1e-15

    z=1.*y
    w=qa-np.log(z)
    u=np.sqrt(w)
    s=(qc+np.log(u))/w
    t=1/(u+qb)
    x=u*(1-s*(0.5+s*qd))-((((q4*t+q3)*t+q2)*t+q1)*t+q0)*t
    t=pa/(pa+x)
    u=t-0.5
    s=(((((((((p22*u+p21)*u+p20)*u+p19)*u+p18)*u+p17)*u+p16)*u+p15)*u+p14)*u+p13)*u+p12
    s=((((((((((((s*u+p11)*u+p10)*u+p9)*u+p8)*u+p7)*u+p6)*u+p5)*u+p4)*u+p3)*u+p2)*u+p1)*u+p0)*t-z*np.exp(x*x-pb)
    x=x+s*(1+x*s)
    return x
## \fn dierfc(y)
# inverse of complimentary error function from MultiNest implementation
# @param y

def ginv(x, mu, sigma):
    # parameter(SqrtTwo=1.414213562d0)
    return mu+sigma*np.sqrt(2)*dierfc(2.*(1.-x))
## \fn ginv(x, mu, sigma)
# take uniform range [0,1] and return normal pdf centered around mu, with width sigma
# sampled with a normal distribution
# @param x random variable sampled by MultiNest, uniformly in range [0,1]
# @param mu central value
# @param sigma width of normal distribution

if __name__=="__main__":
    import pdb
    from pylab import *
    Nsample = 1000
    x = np.linspace(0, 1, Nsample)
    out = ginv(x, 0., 0.1)
    hist(out, np.sqrt(Nsample))
    show()
    pdb.set_trace()
