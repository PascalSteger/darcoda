#!/usr/bin/python
import sys
import math

f = open(sys.argv[1],"r")
m=[];x=[];y=[];z=[]; N=0
for line in f:
    val = line.split()
    m.append(float(val[0]))
    x.append(float(val[1]))
    y.append(float(val[2]))
    z.append(float(val[3]))
    N=N+1

totm=0.0;
xcm=0.0;
ycm=0.0;
zcm=0.0;
for i in range(N):
    xcm += m[i]*x[i]
    ycm += m[i]*y[i]
    zcm += m[i]*z[i]
    totm +=m[i]


xcm = xcm/totm
ycm = ycm/totm
zcm = zcm/totm

print xcm, ycm, zcm
