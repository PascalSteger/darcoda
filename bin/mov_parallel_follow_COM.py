#!/usr/bin/python

import os
import linecache

import threading
semaphore = threading.Semaphore(48)

def run_command(cmd):
    with semaphore:
        os.system(cmd)

N=0
h=0.719
xc = []; yc = []; zc = []; rvir = []
xo =0.5; yo =0.5; zo =0.5; rviro=-1.0
# find center, and zoom region
for i in range(270):
    cdc = "output_"+str(i+1).zfill(5)+"/"
    fn  = cdc + "halom"
    if(not os.path.exists(fn)):
        xc.append(xo)
        yc.append(yo)
        zc.append(zo)
        rvir.append(rviro)
        print "file "+fn+" does not exist"
        continue;
    if(not os.path.getsize(fn) > 0):
        xc.append(xo)
        yc.append(yo)
        zc.append(zo)
        rvir.append(rviro)
        print "file empty"
        continue

    f = open(fn,"r")
    line = linecache.getline(fn,1)
    print line
    val = line.split()
    xc.append(float(val[2-1])); xo = xc[N]
    yc.append(float(val[3-1])); yo = yc[N]
    zc.append(float(val[4-1])); zo = zc[N]
    rvir.append(float(val[5-1])); rviro=rvir[N]
    N = N+1
    f.close()

res    = "512"
dirdmd = "mov/dmdens"
for i in range(270):
    stri = str(i+1).zfill(5)
    cdc = "output_"+stri+"/"
    fout= dirdmd+"/map_"+stri
    # get projections of dm density
    os.system("mkdir "+dirdmd)
    cmd = "part2map -inp "+cdc
    cmd += " -out "+fout
    cmd += " -dir z "
    cmd += " -xmi "+str(xc[i]-rviro)+" -xma "+str(xc[i]+rviro)
    cmd += " -ymi "+str(yc[i]-rviro)+" -yma "+str(yc[i]+rviro)
    cmd += " -zmi "+str(zc[i]-rviro)+" -zma "+str(zc[i]+rviro)
    cmd += " -nx "+res+" -ny "+res
    print "#### "+cmd
    threading.Thread(target=run_command, args=(cmd, )).start()


for i in range(270):
    stri = str(i+1).zfill(5)
    cdc = "output_"+stri+"/"
    fout= dirdmd+"/map_"+stri
    cmd = "map2img.py "+fout+" -o "+fout+".png"
    print "#### "+cmd
    threading.Thread(target=run_command, args=(cmd, )).start()

# get projections of gas density

# get temperature

# get projections of stars density

# get metallicity

