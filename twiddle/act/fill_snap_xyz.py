#!/usr/bin/env python2

## \file
# update MySQL table snapshot after AHF run

# (c) 2014 ETHZ, Pascal Steger, pascal@steger.aero

import sys
import lib.initialize as my

startsnap = int(sys.argv[1])
stopsnap  = int(sys.argv[2])

curhid=1
for nc in range(stopsnap-startsnap+1):
    nsnap=stopsnap-nc
    print('number of snapshot = ',nsnap)

    if(len(my.sql("select snap from halo where snap="+str(nsnap)))==0):
        print("missing snapshot")
        my.sql("update snapshot set "+
               "xm=-1,ym=-1,zm=-1,rm=-1,"+
               "xms=-1,yms=-1,zms=-1,rms=-1,"+
               "atime=now() "+
               "where snap="+str(nsnap))
        curhid=1
        continue
    
    o=my.sql("select xs,ys,zs,rvir,xs_star,ys_star,zs_star,rs_star,proghid"
             +" from halo "
             +"where snap="+str(nsnap)
             +" and hid="+str(curhid))
    x=o[0][0]; y=o[0][1]; z=o[0][2]; r=o[0][3];
    xs=o[0][4];ys=o[0][5];zs=o[0][6];rs=o[0][7];proghid=o[0][8]
    print(x,y,z,r,"progenitor:",proghid)

    my.sql("update snapshot set "+
           "xm="+str(x)+","+
           "ym="+str(y)+","+
           "zm="+str(z)+","+
           "rm="+str(r)+","+
           "xms="+str(xs)+","+
           "yms="+str(ys)+","+
           "zms="+str(zs)+","+
           "rms="+str(rs)+","+
           "atime=now() "+
           "where snap="+str(nsnap))
    curhid=proghid
    
    #my.sql("update snapshot set xmega=")
    #select snapshot.snap,halo.proghid,halo.ss_x from snapshot left join halo on snapshot.snap=halo.snap and halo.hid=1 where halo.snap=270;
