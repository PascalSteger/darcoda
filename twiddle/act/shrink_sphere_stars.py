#!/usr/bin/env python2

## \file
# determine center of mass by iteratively shrinking a sphere

# (c) 2014 ETHZ, Pascal Steger, pascal@steger.aero

import sys
from numpy import array,inner,sqrt
import lib.mysql as mys
import lib.initialize as my

i=len(sys.argv)
if(i!=3):
    print("shrink_sphere.py snap hid")
    exit(0)

snap=int(sys.argv[1])
hid=int(sys.argv[2])

eps=1e-5
frad = 0.9 # shrinked sphere has radius frad*maxr, i.e. frad = 0.90 means 10% smaller

#print("missing file")
if(not mys.exists_snap(snap)):
    print("snapshot "+str(snap)+" missing")
    exit(0)

xc,yc,zc,mvir,rvir=mys.getxyzmr(snap,1)
halo = my.open_file(mys.d(snap)+"stars/stars_"+str(hid)+".dat","r")
x=[];y=[];z=[];m=[]
for line in halo:
    val=line.split()
    m.append(float(val[0]))
    x.append(float(val[1]))
    y.append(float(val[2]))
    z.append(float(val[3]))

x = array(x); y = array(y); z = array(z); m = array(m)
#print(x)

def converged(xc,yc,zc,xc2,yc2,zc2):
    return  abs(xc-xc2)/abs(xc)   < eps \
        and abs(yc-yc2)/abs(yc)   < eps \
        and abs(zc-zc2)/abs(zc)   < eps #\
#        or rc < eps*max(r)
def converged(xc,yc,zc,xc2,yc2,zc2):
    return  abs(xc-xc2)/abs(xc)   < eps \
        and abs(yc-yc2)/abs(yc)   < eps \
        and abs(zc-zc2)/abs(zc)   < eps #\
#        or rc < eps*max(r)
## \fn converged(xc,yc,zc,xc2,yc2,zc2)
# determine whether two consecutive centers of mass are closer than epsilon
# @param xc com1
# @param yc com1
# @param zc com1
# @param xc2 com2
# @param yc2 com2
# @param zc2 com2


M  = sum(m);
xc = inner(m,x)/M; yc = inner(m,y)/M; zc = inner(m,z)/M
r  = sqrt((x-xc)**2+(y-yc)**2+(z-zc)**2)
rc = 0.9*max(r)

xc2 = yc2 = zc2 = -1e99
#while(True):
for i in range(100):
    rc = frad * rc
    r  = sqrt((x-xc)**2+(y-yc)**2+(z-zc)**2)
    order = r.argsort()
    x = x[order]; y = y[order]; z = z[order]; m = m[order]; r = r[order]
    for cut in range(len(r)):
        if(r[cut]>rc):
            break
    
    xc2 = inner(m[:cut],x[:cut])/sum(m[:cut])
    yc2 = inner(m[:cut],y[:cut])/sum(m[:cut])
    zc2 = inner(m[:cut],z[:cut])/sum(m[:cut])
    
    #print(cut,rc,xc2,yc2,zc2)
    
    if(cut<len(m)/100 or converged(xc,yc,zc,xc2,yc2,zc2)):
        print(xc,yc,zc,rc)
        mys.set_ss_stars(snap,hid,xc,yc,zc,rc)
        exit()
    xc,yc,zc = xc2,yc2,zc2
