#!/usr/bin/python
# -*- coding: utf-8 -*-
import initialize as my
import numpy as np
import pdb

def clear():
    clear_sim();
    clear_snapshot();
    clear_halo();
## \fn clear()
# clear simulations, snapshots, and halos

def clear_sim():
    my.sql('DROP TABLE IF EXISTS sim;')
## \fn clear_sim()
# remove table sim

def clear_snapshot():
    my.sql('DROP TABLE IF EXISTS snapshot;')
## \fn clear_snapshot()
# remove table snapshot

def clear_halo():
    my.sql('DROP TABLE IF EXISTS halo;')
## \fn clear_halo()
# remove table halo

def setup():
    # sample statement
    # print('Server version:')
    # cmd = "select version();"
    # my.sqlout(cmd)
    # set up tables
    setup_sim();
    setup_snapshot();
    setup_halo();
    print('sim, snapshot, halo set up')
## \fn setup()
# setup sim, snapshot, halo

def setup_sim():
    my.sql('CREATE TABLE if not exists sim (\
                  name VARCHAR(255) NOT NULL PRIMARY KEY,\
                  active BOOL NOT NULL DEFAULT 0,\
                  dir VARCHAR(255) NOT NULL,\
                  lma INT NOT NULL DEFAULT 14,\
                  nsnap INT NOT NULL,\
                  nstart INT NOT NULL,\
                  nstop INT NOT NULL,\
                  dmonly BOOL DEFAULT FALSE,\
                  atime TIMESTAMP DEFAULT NOW()\
                  );')
    #                  id INT NOT NULL AUTO_INCREMENT,\
    return
## \fn setup_sim()
# create table sim

def setup_snapshot():
    my.sql('CREATE TABLE if not exists snapshot (\
                  snap INT NOT NULL PRIMARY KEY,\
                  sim VARCHAR(255),\
                  a float, z float,\
                  nhalo int,\
                  xm float, ym float, zm float, rm float,\
                  xms float, yms float, zms float, rms float,\
                  atime TIMESTAMP DEFAULT NOW()\
                  );')
## \fn setup_snapshot()
# create table snapshot

def setup_halo():
    my.sql('CREATE TABLE if not exists halo (\
                  id varchar(33) NOT NULL PRIMARY KEY,\
                  snap INT NOT NULL,\
                  sim VARCHAR(255),\
                  buhid INT NOT NULL,\
                  atime TIMESTAMP DEFAULT NOW(),\
                  hid INT NOT NULL DEFAULT 0,\
                  proghid INT NOT NULL DEFAULT 0,\
                  hosthid INT NOT NULL DEFAULT 0,\
                  npart INT,\
                  fmhires float,\
                  xc float, yc float, zc float,\
                  xs float, ys float, zs float, rs float,\
                  vxc float, vyc float, vzc float,\
                  mvir FLOAT, rvir FLOAT, rcore FLOAT,\
                  vmax float, rmax float, sigv float, lambda float,\
                  lx float, ly float, lz float,\
                  a float, eax float, eay float, eaz float,\
                  b float, ebx float, eby float, ebz float,\
                  c float, ecx float, ecy float, ecz float,\
                  ovdens float, Redge float,\
                  nbins int,\
                  Ekin float, Epot float,\
                  mbp_offset float, com_offset float, r2 float,\
                  lambdaE float, v_esc float, Phi0 float,\
                  n_gas int default 0,\
                  xs_gas float, ys_gas float, zs_gas float,\
                  rs_gas float default 0.0,\
                  M_gas float default 0.0, lambda_gas float,\
                  Lx_gas float, Ly_gas float, Lz_gas float,\
                  a_gas float, Eax_gas float, Eay_gas float, Eaz_gas float,\
                  b_gas float, Ebx_gas float, Eby_gas float, Ebz_gas float,\
                  c_gas float, Ecx_gas float, Ecy_gas float, Ecz_gas float,\
                  Ekin_gas float, Epot_gas float, lambdaE_gas float,\
                  n_star int default 0,\
                  xs_star float, ys_star float, zs_star float,\
                  rs_star float default 0.0,\
                  rhalf_star float default 0.0,\
                  M_star float, lambda_star float,\
                  Lx_star float, Ly_star float, Lz_star float,\
                  a_star float,\
                  Eax_star float, Eay_star float, Eaz_star float,\
                  b_star float,\
                  Ebx_star float, Eby_star float, Ebz_star float,\
                  c_star float,\
                  Ecx_star float, Ecy_star float, Ecz_star float,\
                  Ekin_star float, Epot_star float, lambdaE_star float,\
                  m_dm float\
                  );')
## \fn setup_snapshot()
# create table snapshot

def get_nstart():
    cmd = "SELECT nstart FROM sim WHERE active=1;"
    if(len(my.sql(cmd))<1):
        return 1
    nstart = my.sql(cmd)[0][0]
    return nstart
## \fn get_nstart()
# get start output number from active simulation

def get_lma():
    cmd = "SELECT lma FROM sim WHERE active=1;"
    if(len(my.sql(cmd))<1):
        return 19
    lma = my.sql(cmd)[0][0]
    return lma


def set_nstart(st):
    cmd = "UPDATE sim SET nstart="+str(st)+" WHERE active=1;"
    my.sql(cmd)
    return


def set_lma(st):
    cmd = "UPDATE sim SET lma="+str(st)+" WHERE active=1;"
    my.sql(cmd)
    return


def get_nstop():
    cmd = "SELECT nstop FROM sim WHERE active=1;"
    if(len(my.sql(cmd))<1):
        return 1
    nstop = my.sql(cmd)[0][0]
    return nstop


def set_nstop(st):
    cmd = "UPDATE sim SET nstop="+str(st)+" WHERE active=1;"
    my.sql(cmd)
    return

def get_active_sim():
    cmd = "SELECT name FROM sim WHERE active=1;"
    # print("len = ",len(my.sql(cmd))
    if(len(my.sql(cmd))<1):
        return " "
    sim = my.sql(cmd)[0][0]
    return sim

def set_active_sim(st):
    cmd = "UPDATE sim SET active=1 WHERE name like '"+st+"';"
    my.sql(cmd)
    cmd = "UPDATE sim SET active=0 WHERE name not like '"+st+"';"
    my.sql(cmd)
    return

def get_range():
    nstart = get_nstart()
    nstop  = get_nstop()
    return nstart,nstop

def get_active():
    sim = get_active_sim()
    nstart,nstop = get_range()
    return sim,nstart,nstop

def fill_sim(name,nsnap,dmonly):
    simdir = "/library/home/psteger/sci/sim/"+name
    print(dmonly)
    cmd  = "INSERT INTO sim(name,active,dir,nsnap,nstart,nstop,dmonly) "
    cmd += "VALUES ("
    # default values 1, nsnap for nstart, nstop
    cmd += "'"+name+"','1','"+simdir+"','"+str(nsnap)+"','"+str(1)+"','"+str(nsnap)+"',"+str(dmonly)+") "
    cmd += "ON DUPLICATE KEY UPDATE atime=now(),active='1',dmonly="+str(dmonly)+";"
    print(cmd)
    my.sql(cmd)
    cmd = "UPDATE sim SET active='0' WHERE name not like '"+name+"';"
    print(cmd)
    my.sql(cmd)
    return

def print_sim(name):
	cmd = "SELECT * FROM sim WHERE name='"+name+"';"
	out = my.sql(cmd)
	print(out)

def fill_halo(snap, buhid, val, co, cu):
    cmd  = "INSERT INTO halo(id, snap, buhid) "
    cmd += " VALUES ('"+my.md5(snap, buhid)
    cmd += "','"+str(snap)+"','"+str(buhid)+"')";
    cmd += " ON DUPLICATE KEY UPDATE atime=now();"
    my.sqlcu(cmd,co,cu)
    fmhires = float(val[38-1])
    mvir = float(val[4-1]);
    if(is_dmonly()):
        m_gas = 0.0
        m_star = 0.0
    else:
        m_gas = float(val[45-1])
        m_star = 0.0
        m_star=float(val[65-1])
    m_dm  = mvir-m_gas-m_star
    cmd  = "UPDATE halo SET "
    cmd += "npart='"      +val[5-1]+"', "
    cmd += "fmhires='"    +str(fmhires)+"', "
    cmd += "xc='"+val[6-1]+"', yc='"+val[7-1]+"', zc='"+val[8-1]+"', "
    cmd += "vxc='"+val[9-1]+"',vyc='"+val[10-1]+"',vzc='"+val[11-1]+"', "
    cmd += "mvir='"+val[4-1]+"',  rvir='"+val[12-1]+"', "
    cmd += "vmax='"+val[17-1]+"', rmax='"+val[13-1]+"', "
    cmd += "sigv='"+val[19-1]+"', lambda='"+val[20-1]+"', "
    cmd += "lx='"+val[22-1]+"',ly='"+val[23-1]+"',lz='"+val[24-1]+"', "
    cmd += "a='1', eax='"+val[27-1]+"',eay='"+val[28-1]+"',eaz='"+val[29-1]+"', "
    cmd += "b='"+val[25-1]+"', ebx='"+val[30-1]+"',eby='"+val[31-1]+"',ebz='"+val[32-1]+"', "
    cmd += "c='"+val[26-1]+"', ecx='"+val[33-1]+"',ecy='"+val[34-1]+"',ecz='"+val[35-1]+"', "
    cmd += "ovdens='"+val[36-1]+"', "
    cmd += "redge='"+val[13-1]+"', "
    cmd += "nbins='"+val[37-1]+"', "
    cmd += "ekin='"+val[39-1]+"',epot='"+val[40-1]+"', "
#    cmd += "mbp_offset='"+val[16-1]+"', com_offset='"+val[17-1]+"', "
    cmd += "r2='"+val[14-1]+"', "
    cmd += "lambdae='"+val[21-1]+"', "
    cmd += "v_esc='"+val[18-1]+"', "
    cmd += "phi0='"+val[42-1]+"', "
    if(not is_dmonly()):
        cmd += "n_gas='"+val[44-1]+"', m_gas='"+val[45-1]+"', "
        cmd += "lambda_gas='"+val[46-1]+"', "
        cmd += "lx_gas='"+val[48-1]+"',ly_gas='"+val[49-1]+"',lz_gas='"+val[50-1]+"', "
        cmd += "a_gas='1', eax_gas='"+val[53-1]+"',eay_gas='"+val[54-1]+"',eaz_gas='"+val[55-1]+"', "
        cmd += "b_gas='"+val[51-1]+"', ebx_gas='"+val[56-1]+"',eby_gas='"+val[57-1]+"',ebz_gas='"+val[58-1]+"', "
        cmd += "c_gas='"+val[52-1]+"', ecx_gas='"+val[59-1]+"',ecy_gas='"+val[60-1]+"',ecz_gas='"+val[61-1]+"', "
        cmd += "ekin_gas='"+val[62-1]+"', epot_gas='"+val[63-1]+"', "
        cmd += "lambdae_gas='"+val[47-1]+"', "
        cmd += "n_star='"+val[64-1]+"', "
        cmd += "m_star='"+val[65-1]+"', "
        cmd += "lambda_star='"+val[66-1]+"', "
        cmd += "lx_star='"+val[68-1]+"',ly_star='"+val[69-1]+"',lz_star='"+val[70-1]+"', "
        cmd += "a_star='1', eax_star='"+val[73-1]+"',eay_star='"+val[74-1]+"',eaz_star='"+val[75-1]+"', "
        cmd += "b_star='"+val[71-1]+"', ebx_star='"+val[76-1]+"',eby_star='"+val[77-1]+"',ebz_star='"+val[78-1]+"', "
        cmd += "c_star='"+val[72-1]+"', ecx_star='"+val[79-1]+"',ecy_star='"+val[80-1]+"',ecz_star='"+val[81-1]+"', "
        cmd += "ekin_star='"+val[82-1]+"', epot_star='"+val[83-1]+"', "
        cmd += "lambdae_star='"+val[67-1]+"', "
        cmd += "m_dm='"+str(m_dm)+"', "
    cmd += "sim='"+get_active_sim()+"', "
    cmd += "atime=now() "
    cmd += "WHERE id='"+my.md5(snap,buhid)+"';"
    my.sqlcu(cmd,co,cu)

def fill_snapshot(snap):
    a=str(my.read_a(snap))
    print('a read: ', a)
    z=str(my.read_z(snap))
    print('z read: ', z)
    nhalo=str(get_nhalo(snap))
    print('nhalo:', nhalo)
    cmd  = "INSERT INTO snapshot(snap,a,z,nhalo,sim,atime) "
    cmd += " VALUES ('"
    cmd += str(snap)+"','"+a+"','"+z+"','"+nhalo+"','"+get_active_sim()
    cmd += "', now())"
    cmd += " ON DUPLICATE KEY UPDATE atime=now();"
    my.sql(cmd)
    my.sql("ALTER TABLE snapshot ORDER BY snap DESC;")

def exists_snap(snap):
    c = my.sql("SELECT snap FROM snapshot WHERE snap="+str(snap)+";")
    if(len(c)==0):
        return False;
    else:
        return True;

def get_z(snap):
    c = my.sql("SELECT z FROM snapshot WHERE snap="+str(snap)+";")
    return c[0][0];

def set_ss(snap,hid,xs,ys,zs,rs):
    cmd  = "UPDATE halo SET "
    cmd += "xs="+str(xs)+", "
    cmd += "ys="+str(ys)+", "
    cmd += "zs="+str(zs)+", "
    cmd += "rs="+str(rs)+", "
    cmd += "atime=now() WHERE snap="+str(snap)+" and hid ="+str(hid)+";"
    my.sql(cmd)
    return

def set_ss_stars(snap,hid,xs_star,ys_star,zs_star,rs_star):
    cmd  = "UPDATE halo SET "
    cmd += "xs_star="+str(xs_star)+", "
    cmd += "ys_star="+str(ys_star)+", "
    cmd += "zs_star="+str(zs_star)+", "
    cmd += "rs_star="+str(rs_star)+", "
    cmd += "atime=now() WHERE snap="+str(snap)+" and hid ="+str(hid)+";"
    my.sql(cmd)

def physical_xcm(snap, h):
    actsim = my.sql("SELECT name FROM sim WHERE active=1;")
    actsim = actsim[0][0]
    sh  = str(h);
    cmd  = "UPDATE halo SET xc=xc/"+sh+", yc=yc/"+sh+", zc=zc/"+sh
    cmd += ", rvir=rvir/"+sh+", atime=now() WHERE snap="+str(snap)+" and sim like '"+actsim+"';"
    my.sql(cmd)
    return
    ## \fn physical_xcm(snap, h)
    # convert all values in DB pertaining to snapshot snap from physical AHF units to ramses code units [0,1]
    # @param snap snapshot id
    # @param h = H0/100

def exclude(snap, npart, fmhires, bndry):
    cmd  = "DELETE from halo "
    cmd += "WHERE snap="+str(snap)+" and (npart < "+str(npart)
    cmd += " or fmhires <= "+str(fmhires)
    cmd += " or xc < "+str(bndry)+" or xc > "+str(1-bndry)
    cmd += " or yc < "+str(bndry)+" or yc > "+str(1-bndry)
    cmd += " or zc < "+str(bndry)+" or zc > "+str(1-bndry)+");"
    my.sql(cmd)
    my.sql("ALTER TABLE halo ORDER BY mvir DESC;")
    my.sql("UPDATE halo SET hid=0 WHERE id!='' AND snap="+str(snap)+";")
    old_buhid=0;
    for i in range(get_nhalo(snap)):
        my.sql("UPDATE halo SET hid="+str(i+1)
               +" WHERE id!='' and hid>="+str(i)
               +" AND buhid!="+str(old_buhid)
               +" AND snap="+str(snap)+";")
        old_buhid=my.sql("SELECT buhid from halo WHERE snap="+str(snap)+" AND hid>"+str(i)
                         +" limit 1;")[0][0]
    my.sql("ANALYZE TABLE halo;")
    return
    ## \fn exclude(snap, npart, fmhires, bndry)
    # delete all entries which have
    # @param snap snapshot id
    # @param npart less than [] particles
    # @param fmhires a fraction of less than []% in high-resolution particles
    # @param bndry and which lie closer than []kpc to the boundaries

def set_rhalf_star(snap,hid,r):
    my.sql("UPDATE halo SET rhalf_star="+str(r)+
           " WHERE snap="+str(snap)+
           " AND hid="+str(hid)+";")
    return

def getxyzmr(snap, typ):
    actsim=my.sql("SELECT name FROM sim WHERE active=1;")
    actsim=actsim[0][0]
    if(typ==1):
        xc=my.sql("SELECT xc FROM halo WHERE snap="+str(snap)+" AND sim LIKE '"+actsim+"';")
        yc=my.sql("SELECT yc FROM halo WHERE snap="+str(snap)+" AND sim LIKE '"+actsim+"';")
        zc=my.sql("SELECT zc FROM halo WHERE snap="+str(snap)+" AND sim LIKE '"+actsim+"';")
    if(typ==2):
        xc=my.sql("SELECT xs FROM halo WHERE snap="+str(snap)+" AND sim LIKE '"+actsim+"';")
        yc=my.sql("SELECT ys FROM halo WHERE snap="+str(snap)+" AND sim LIKE '"+actsim+"';")
        zc=my.sql("SELECT zs FROM halo WHERE snap="+str(snap)+" AND sim LIKE '"+actsim+"';")
    mvir=my.sql("SELECT mvir FROM halo WHERE snap="+str(snap)+" AND sim LIKE '"+actsim+"';")
    rvir=my.sql("SELECT rvir FROM halo WHERE snap="+str(snap)+" AND sim LIKE '"+actsim+"';")
    x=[];y=[];z=[];r=[];m=[]
    for i in range(len(xc)):
        x.append(xc[i][0])
        y.append(yc[i][0])
        z.append(zc[i][0])
        r.append(rvir[i][0])
        m.append(mvir[i][0])
    x=np.array(x)
    y=np.array(y)
    z=np.array(z)
    r=np.array(r)
    m=np.array(m)
    order = m.argsort()
    x = x[order]
    y = y[order]
    z = z[order]
    r = r[order]
    m = m[order]
    return x[::-1],y[::-1],z[::-1],m[::-1],r[::-1]

def get_M_star(snap):
    mstar=my.sql("SELECT M_star FROM halo WHERE snap="+str(snap)+";")
    mvir=my.sql("SELECT mvir FROM halo WHERE snap="+str(snap)+";")
    ms=[]; mv=[]
    for i in range(len(mstar)):
        ms.append(mstar[i][0])
        mv.append(mvir[i][0])
    mv = np.array(mv)
    ms = np.array(ms)
    order = mv.argsort()
    mv = mv[order]
    ms = ms[order]
    return ms[::-1]

def get_M_dm(snap):
    mvir=my.sql("SELECT mvir FROM halo WHERE snap="+str(snap)+";")
    mdm = my.sql("SELECT mvir-M_star-M_gas FROM halo WHERE snap="+str(snap)+";")
    md=[]; mv = []
    for i in range(len(mdm)):
        md.append(mdm[i][0])
        mv.append(mvir[i][0])
    md = np.array(md)
    mv = np.array(mv)
    order = mv.argsort()
    mv = mv[order]
    md = md[order]
    return md[::-1]

def get_is_sub(snap):
    issub=my.sql("SELECT hosthid FROM halo WHERE snap="+str(snap)+";")
    mvir=my.sql("SELECT mvir from halo WHERE snap="+str(snap)+";")
    iss=[]; mv =[]
    for i in range(len(issub)):
        iss.append(issub[i][0]>0)
        mv.append(mvir[i][0])
    mv = np.array(mv)
    iss = np.array(iss)
    order = mv.argsort()
    iss = iss[order]
    return iss[::-1]

def get_rhalf_star(snap):
    r2s=my.sql("SELECT rhalf_star FROM halo WHERE snap="+str(snap)+";")
    mvir=my.sql("SELECT mvir FROM halo WHERE snap="+str(snap)+";")
    r=[];m=[]
    for i in range(len(r2s)):
        r.append(r2s[i][0])
        m.append(mvir[i][0])
    r = np.array(r)
    m = np.array(m)
    order = m.argsort()
    r = r[order]
    return r[::-1]

def getxyzrsnap(snap):
    wc="FROM snapshot WHERE snap="+str(snap)+";"
    xm=my.sql("SELECT xs  "+wc)
    ym=my.sql("SELECT ys  "+wc)
    zm=my.sql("SELECT zs  "+wc)
    rm=my.sql("SELECT rs "+wc)
    sm=my.sql("SELECT snap "+wc)
    x=[];y=[];z=[];r=[];s=[]
    for i in range(len(xm)):
        x.append(xm[i][0])
        y.append(ym[i][0])
        z.append(zm[i][0])
        r.append(rm[i][0])
        s.append(sm[i][0])
    return x,y,z,r,s

def mt_xyzrsnap(start,stop):
    wc="FROM snapshot WHERE snap>="+str(start)+" AND snap<="+str(stop)+";"
    xm=my.sql("SELECT xm  "+wc)
    ym=my.sql("SELECT ym  "+wc)
    zm=my.sql("SELECT zm  "+wc)
    rm=my.sql("SELECT rm "+wc)
    sm=my.sql("SELECT snap "+wc)
    x=[];y=[];z=[];r=[];s=[];
    for i in range(len(xm)):
        x.append(xm[i][0])
        y.append(ym[i][0])
        z.append(zm[i][0])
        r.append(rm[i][0])
        s.append(sm[i][0])
    return x,y,z,r,s

def mt_xyzrsnap_stars(start,stop):
    wc=" FROM snapshot WHERE snap>="+str(start)+" AND snap<="+str(stop)+";"
    xm=my.sql("SELECT xms "+wc)
    ym=my.sql("SELECT yms "+wc)
    zm=my.sql("SELECT zms "+wc)
    rm=my.sql("SELECT rms "+wc)
    sm=my.sql("SELECT snap "+wc)
    x=[];y=[];z=[];r=[];s=[];
    for i in range(len(xm)):
        x.append(xm[i][0])
        y.append(ym[i][0])
        z.append(zm[i][0])
        r.append(rm[i][0])
        s.append(sm[i][0])
    return x,y,z,r,s

def getxyzrstars_hid(snap,hid):
    wc=" FROM halo WHERE snap="+str(snap)+" AND hid="+str(hid)+";"
    xm=my.sql("SELECT xs_star "+wc)
    ym=my.sql("SELECT ys_star "+wc)
    zm=my.sql("SELECT zs_star "+wc)
    mvir=my.sql("SELECT mvir "+wc)
    rm=my.sql("SELECT rs_star "+wc)
    xms=[];yms=[];zms=[];mvirs=[];rms=[]
    for i in range(len(mvir)):
        xms.append(xm[i][0])
        yms.append(ym[i][0])
        zms.append(zm[i][0])
        mvirs.append(mvir[i][0])
        rms.append(rm[i][0])
    xms=np.array(xms); yms=np.array(yms); zms=np.array(zms); mvirs=np.array(mvirs);rms=np.array(rms)
    order = mvirs.argsort()
    xms=xms[order]; yms=yms[order]; zms=zms[order]; rms=rms[order]; mvirs=mvirs[order]
    return xms[-1],yms[-1],zms[-1],rms[-1]

def getxyzrstars(snap,typ):
    if(typ==1):
        xs=my.sql("SELECT xc FROM halo WHERE snap="+str(snap)+";")
        ys=my.sql("SELECT yc FROM halo WHERE snap="+str(snap)+";")
        zs=my.sql("SELECT zc FROM halo WHERE snap="+str(snap)+";")
        rs=my.sql("SELECT rvir FROM halo WHERE snap="+str(snap)+";")
        mv=my.sql("SELECT mvir FROM halo WHERE snap="+str(snap)+";")
    if(typ==2):
        xs=my.sql("SELECT xs_star FROM halo WHERE snap="+str(snap)+";")
        ys=my.sql("SELECT ys_star FROM halo WHERE snap="+str(snap)+";")
        zs=my.sql("SELECT zs_star FROM halo WHERE snap="+str(snap)+";")
        rs=my.sql("SELECT rs_star FROM halo WHERE snap="+str(snap)+";")
        mv=my.sql("SELECT mvir FROM halo WHERE snap="+str(snap)+";")
    x=[]; y=[]; z=[]; r=[]; m=[]
    for i in range(len(xs)):
        x.append(xs[i][0])
        y.append(ys[i][0])
        z.append(zs[i][0])
        r.append(rs[i][0])
        m.append(mv[i][0])
    x=np.array(x);y=np.array(y);z=np.array(z);r=np.array(r);m=np.array(m)
    order = m.argsort()
    x=x[order]; y=y[order]; z=z[order]; r=r[order]; m=m[order]
    return x[::-1],y[::-1],z[::-1],r[::-1]

def get_nhalo(snap):
    nhalo=my.sql("SELECT atime FROM halo WHERE snap="+str(snap)+";")
    return len(nhalo)

def getNmax():
    d=my.sql("SELECT nsnap FROM sim WHERE active=1;")
    return d[0][0]

def is_dmonly():
    dmonly=my.sql("SELECT dmonly FROM sim WHERE active=1;")
    return dmonly[0][0]

def d(nsnap):
    d = my.sql("SELECT dir FROM sim WHERE active=1;")
    d = d[0][0]+"/output_"+str(nsnap).zfill(5)+"/"
    return d

def simdir():
    d = my.sql("SELECT dir FROM sim WHERE active=1;")
    d = d[0][0]+"/"
    return d

def snap_exists(n):
    exists=my.file_exists(d(n)+"/info_"+str(n).zfill(5)+".txt");
    return exists

def get_sims():
    cmd="select name from sim;"
    o=my.sql(cmd)
    if(len(o)<1):
        o=[['find...']]
    p=np.array(o)
    q=np.transpose(p)[0]
    return q
