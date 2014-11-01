#!/usr/bin/python
#get clusters with > 10**7 Msun, and not subhalos

import os

for i in range(270):
    cdc = "output_"+str(i+1).zfill(5)+"/"
    if(not os.path.exists(cdc+"halos")):
        continue;
    f = open(cdc+'halos',"r")
    out = open(cdc+'clusters',"w")
    xc = []; yc = []; zc = []; mvir = []; rvir = []
    N=0
    h=0.719
    for line in f:
        if(line[0]=="#"):
            continue
        N = N+1
        val = line.split()
        xc.append(float(val[3-1])/h)
        yc.append(float(val[4-1])/h)
        zc.append(float(val[5-1])/h)
        mvir.append(float(val[9-1]))
        rvir.append(float(val[10-1])/h/1000)
        
    for i in xrange(N):
        if(mvir[i]<10**7):
            continue
        # check whether it is a possible subhalo, then exclude it
        inside = 0
        for j in xrange(N):
            # don't look at same halo
            if(i==j):
                continue
            # don't look at smaller halo
            if(rvir[j]<rvir[i]):
                continue
            dx = xc[i]-xc[j]
            dy = yc[i]-yc[j]
            dz = zc[i]-zc[j]
            d2 = dx**2+dy**2+dz**2
            if(d2 < rvir[j]**2):
                inside = 1
                continue
        if(inside==0):
            out.write(str(xc[i])+" "+str(yc[i])+" "+str(zc[i])+" "+str(rvir[i])+" "+str(mvir[i])+"\n")

    f.close()
    out.close()
