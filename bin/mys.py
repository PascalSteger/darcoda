#!/usr/bin/python
# -*- coding: utf-8 -*-

import initialize as my

def clear():
    clear_sim();
    clear_snapshot();
    clear_halo();

def clear_sim():
    my.sql('DROP table if exists sim;')

def clear_snapshot():
    my.sql('DROP table if exists snapshot;')

def clear_halo():
    my.sql('DROP table if exists halo;')

def setup():
    # sample statement
    print 'Server version:'
    cmd = "select version();"
    my.sqlout(cmd)

    #set up tables
    setup_sim();
    setup_snapshot();
    setup_halo();

def setup_sim():
    my.sql('CREATE TABLE if not exists sim (\
                  id INT NOT NULL PRIMARY KEY AUTO_INCREMENT,\
                  dir VARCHAR(255) NOT NULL,\
                  nsnap INT NOT NULL,\
                  nstart INT NOT NULL,\
                  nstop INT NOT NULL,\
                  atime TIMESTAMP DEFAULT NOW()\
                  );')

def setup_snapshot():
    my.sql('CREATE TABLE if not exists snapshot (\
                  snap INT NOT NULL PRIMARY KEY,\
                  a float, z float,\
                  nhalo int,\
                  xm float, ym float, zm float, rm float,\
                  xms float, yms float, zms float, rms float,\
                  atime TIMESTAMP DEFAULT NOW()\
                  );')

def setup_halo():
    my.sql('CREATE TABLE if not exists halo (\
                  id varchar(33) NOT NULL PRIMARY KEY,\
                  snap INT NOT NULL,\
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
                  mvir FLOAT, rvir FLOAT,\
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
                  Ekin_star float, Epot_star float, lambdaE_star float\
                  );')

    #my.sqlout(cu,'describe halo')
    #fill table
    #my.sql(cu,'truncate table halo')
    #my.sql(cu,"INSERT INTO halo(buhid,snap,mvir) VALUES('1','270','134.03')")
    #my.sql(cu,"INSERT INTO halo(buhid,snap,mvir) VALUES('20','270','4.03')")
    #my.sql(cu,"INSERT INTO halo(buhid,snap,mvir) VALUES('1','269','1.402')")

    #read out one property
    #my.sqlout(cu,"select id,snap,mvir,rvir from halo where snap=270")

def fill_sim(simdir,nsnap,nstart,nstop):
    cmd  = "INSERT INTO sim(dir,nsnap,nstart,nstop) "
    cmd += "VALUES ("
    cmd += "'"+simdir+"',"+str(nsnap)+","+str(nstart)+","+str(nstop)+");"
    my.sql(cmd)

def fill_halo(snap,buhid,val):
    cmd  = "INSERT INTO halo(id,snap,buhid) "
    cmd += " VALUES ('"+my.md5(snap,buhid)
    cmd += "','"+str(snap)+"','"+str(buhid)+"')";
    cmd += " ON DUPLICATE KEY UPDATE atime=now();"
    my.sql(cmd)

    cmd  = "UPDATE halo SET "
    cmd += "npart='"      +val[0]+"', "
    cmd += "fmhires='"    +val[1]+"', "
    cmd += "xc='"+val[2]+"', yc='"+val[3]+"', zc='"+val[4]+"', "
    cmd += "vxc='"+val[5]+"',vyc='"+val[6]+"',vzc='"+val[7]+"', "
    cmd += "mvir='"+val[8]+"',  rvir='"+val[9]+"', "
    cmd += "vmax='"+val[10]+"', rmax='"+val[11]+"', "
    cmd += "sigv='"+val[12]+"', "
    cmd += "lambda='"+val[13]+"', "
    cmd += "lx='"+val[14]+"',ly='"+val[15]+"',lz='"+val[16]+"', "
    cmd += "a='"+val[17]+"', "
    cmd += "eax='"+val[18]+"',eay='"+val[19]+"',eaz='"+val[20]+"', "
    cmd += "b='"+val[21]+"', "
    cmd += "ebx='"+val[22]+"',eby='"+val[23]+"',ebz='"+val[24]+"', "
    cmd += "c='"+val[25]+"', "
    cmd += "ecx='"+val[26]+"',ecy='"+val[27]+"',ecz='"+val[28]+"', "
    cmd += "ovdens='"+val[29]+"', "
    cmd += "redge='"+val[30]+"', "
    cmd += "nbins='"+val[31]+"', "
    cmd += "ekin='"+val[32]+"',epot='"+val[33]+"', "
    cmd += "mbp_offset='"+val[34]+"', "
    cmd += "com_offset='"+val[35]+"', "
    cmd += "r2='"+val[36]+"', "
    cmd += "lambdae='"+val[37]+"', "
    cmd += "v_esc='"+val[38]+"', "
    cmd += "phi0='"+val[39]+"', "
    cmd += "n_gas='"+val[40]+"', m_gas='"+val[41]+"', "
    cmd += "lambda_gas='"+val[42]+"', "
    cmd += "lx_gas='"+val[43]+"',ly_gas='"+val[44]+"',lz_gas='"+val[45]+"', "
    cmd += "a_gas='"+val[46]+"', "
    cmd += "eax_gas='"+val[47]+"',eay_gas='"+val[48]+"',eaz_gas='"+val[49]+"', "
    cmd += "b_gas='"+val[50]+"', "
    cmd += "ebx_gas='"+val[51]+"', "
    cmd += "eby_gas='"+val[52]+"', "
    cmd += "ebz_gas='"+val[53]+"', "
    cmd += "c_gas='"+val[54]+"', "
    cmd += "ecx_gas='"+val[55]+"', "
    cmd += "ecy_gas='"+val[56]+"', "
    cmd += "ecz_gas='"+val[57]+"', "
    cmd += "ekin_gas='"+val[58]+"', "
    cmd += "epot_gas='"+val[59]+"', "
    cmd += "lambdae_gas='"+val[60]+"', "
    cmd += "n_star='"+val[61]+"', "
    cmd += "m_star='"+val[62]+"', "
    cmd += "lambda_star='"+val[63]+"', "
    cmd += "lx_star='"+val[64]+"', "
    cmd += "ly_star='"+val[65]+"', "
    cmd += "lz_star='"+val[66]+"', "
    cmd += "a_star='"+val[67]+"', "
    cmd += "eax_star='"+val[68]+"', "
    cmd += "eay_star='"+val[69]+"', "
    cmd += "eaz_star='"+val[70]+"', "
    cmd += "b_star='"+val[71]+"', "
    cmd += "ebx_star='"+val[72]+"', "
    cmd += "eby_star='"+val[73]+"', "
    cmd += "ebz_star='"+val[74]+"', "
    cmd += "c_star='"+val[75]+"', "
    cmd += "ecx_star='"+val[76]+"', "
    cmd += "ecy_star='"+val[77]+"', "
    cmd += "ecz_star='"+val[78]+"', "
    cmd += "ekin_star='"+val[79]+"', "
    cmd += "epot_star='"+val[80]+"', "
    cmd += "lambdae_star='"+val[81]+"', "
    cmd += "atime=now() "
    cmd += "WHERE id='"+my.md5(snap,buhid)+"';"
    my.sql(cmd)

def fill_snapshot(snap):
    a=str(my.read_a(snap))
    print 'a read: ',a
    z=str(my.read_z(snap))
    print 'z read: ',z
    nhalo=str(get_nhalo(snap))
    print 'nhalo:',nhalo

    cmd  = "INSERT INTO snapshot(snap,a,z,nhalo,atime) "
    cmd += " VALUES ('"
    cmd += str(snap)+"','"+a+"','"+z+"','"+nhalo
    cmd += "', now())"
    cmd += " ON DUPLICATE KEY UPDATE atime=now();"
    my.sql(cmd)

    my.sql("alter table snapshot order by snap desc;")



def exists_snap(snap):
    c = my.sql("SELECT snap from snapshot where snap="+str(snap)+";")

    if(len(c)==0):
        return False;
    else:
        return True;

def get_z(snap):
    c = my.sql("select z from snapshot where snap="+str(snap)+";")
    return c[0][0];


def set_ss(snap,hid,xs,ys,zs,rs):
    cmd  = "UPDATE halo SET "
    cmd += "xs="+str(xs)+", "
    cmd += "ys="+str(ys)+", "
    cmd += "zs="+str(zs)+", "
    cmd += "rs="+str(rs)+", "
    cmd += "atime=now() WHERE snap="+str(snap)+" and hid ="+str(hid)+";"
    my.sql(cmd)

def set_ss_stars(snap,hid,xs_star,ys_star,zs_star,rs_star):
    cmd  = "UPDATE halo SET "
    cmd += "xs_star="+str(xs_star)+", "
    cmd += "ys_star="+str(ys_star)+", "
    cmd += "zs_star="+str(zs_star)+", "
    cmd += "rs_star="+str(rs_star)+", "
    cmd += "atime=now() WHERE snap="+str(snap)+" and hid ="+str(hid)+";"
    my.sql(cmd)

    
def physical_xcm(snap,h):
    sh  = str(h);     sh3 = str(1000*h)
    cmd  = "UPDATE halo SET xc=xc/"+sh+", yc=yc/"+sh+", zc=zc/"+sh
    cmd += ", rvir=rvir/"+sh3+", atime=now() WHERE snap="+str(snap)+";"
    my.sql(cmd)

def exclude(snap,npart,fmhires,bndry):
    cmd  = "DELETE from halo "
    cmd += "WHERE snap="+str(snap)+" and (npart < "+str(npart)
    cmd += " or fmhires < "+str(fmhires)
    cmd += " or xc < "+str(bndry)+" or xc > "+str(1-bndry)
    cmd += " or yc < "+str(bndry)+" or yc > "+str(1-bndry)
    cmd += " or zc < "+str(bndry)+" or zc > "+str(1-bndry)+");"
    my.sql(cmd)

    my.sql("alter table halo order by mvir desc;")
    my.sql("update halo set hid=0 where id!='' and snap="+str(snap)+";")
    old_buhid=0;
    for i in range(get_nhalo(snap)):
        my.sql("update halo set hid="+str(i+1)
               +" where id!='' and hid>="+str(i)
               +" and buhid!="+str(old_buhid)
               +" and snap="+str(snap)+";")
        old_buhid=my.sql("select buhid from halo where snap="+str(snap)+" and hid>"+str(i)
                         +" limit 1;")[0][0]
#        print old_buhid
    my.sql("analyze table halo;")

def set_rhalf_star(snap,hid,r):
    my.sql("update halo set rhalf_star="+str(r)+
           " where snap="+str(snap)+
           " and hid="+str(hid)+";")

def getxyzmr(snap,typ):

    if(typ==1):
        xc=my.sql("select xc from halo where snap="+str(snap)+";")
        yc=my.sql("select yc from halo where snap="+str(snap)+";")
        zc=my.sql("select zc from halo where snap="+str(snap)+";")
    if(typ==2):
        xc=my.sql("select xs from halo where snap="+str(snap)+";")
        yc=my.sql("select ys from halo where snap="+str(snap)+";")
        zc=my.sql("select zs from halo where snap="+str(snap)+";")
    mvir=my.sql("select mvir from halo where snap="+str(snap)+";")
    rvir=my.sql("select rvir from halo where snap="+str(snap)+";")

    x=[];y=[];z=[];r=[];m=[]
    for i in range(len(xc)):
        x.append(xc[i][0])
        y.append(yc[i][0])
        z.append(zc[i][0])
        r.append(rvir[i][0])
        m.append(mvir[i][0])
        
    return x,y,z,m,r

def get_M_star(snap):

    mstar=my.sql("select M_star from halo where snap="+str(snap)+";")
    m=[]
    for i in range(len(mstar)):
        m.append(mstar[i][0])
    return m

def get_M_dm(snap):

    mstar=my.sql("select mvir-M_star-M_gas from halo where snap="+str(snap)+";")
    m=[]
    for i in range(len(mstar)):
        m.append(mstar[i][0])
        
    return m

def get_is_sub(snap):

    issub=my.sql("select hosthid from halo where snap="+str(snap)+";")
    iss=[]
    for i in range(len(issub)):
        iss.append(issub[i][0]>0)
    return iss

def get_rhalf_star(snap):

    r2s=my.sql("select rhalf_star from halo where snap="+str(snap)+";")
    r=[]
    for i in range(len(r2s)):
        r.append(r2s[i][0])
        
    return r

def getxyzrsnap(start,stop):

    wc="from snapshot where snap>="+str(start)+" and snap<="+str(stop)+";"
    xm=my.sql("select xm  "+wc)
    ym=my.sql("select ym  "+wc)
    zm=my.sql("select zm  "+wc)
    rm=my.sql("select rm "+wc)
    sm=my.sql("select snap "+wc)

    x=[];y=[];z=[];r=[];s=[];sn=[]
    for i in range(len(xm)):
        x.append(xm[i][0])
        y.append(ym[i][0])
        z.append(zm[i][0])
        r.append(rm[i][0])
        s.append(sm[i][0])
        
    return x,y,z,r,s


def getxyzrsnap_stars(start,stop):

    wc=" from halo where snap>="+str(start)+" and snap<="+str(stop)+";"
    xm=my.sql("select xs_star "+wc)
    ym=my.sql("select ys_star "+wc)
    zm=my.sql("select zs_star "+wc)
    rm=my.sql("select rs_star "+wc)
    sm=my.sql("select snap "+wc)

    x=[];y=[];z=[];r=[];s=[];sn=[]
    for i in range(len(xm)):
        x.append(xm[i][0])
        y.append(ym[i][0])
        z.append(zm[i][0])
        r.append(rm[i][0])
        s.append(sm[i][0])
        
    return x,y,z,r,s

def getxyzrstars_hid(snap,hid):

    wc=" from halo where snap="+str(snap)+" and hid="+str(hid)+";"
    xm=my.sql("select xs_star "+wc)
    ym=my.sql("select ys_star "+wc)
    zm=my.sql("select zs_star "+wc)
    rm=my.sql("select rs_star "+wc)

    return xm[0][0],ym[0][0],zm[0][0],rm[0][0]

def getxyzrstars(snap,typ):

    if(typ==1):
        xs=my.sql("select xc from halo where snap="+str(snap)+";")
        ys=my.sql("select yc from halo where snap="+str(snap)+";")
        zs=my.sql("select zc from halo where snap="+str(snap)+";")
        rs=my.sql("select rvir from halo where snap="+str(snap)+";")
    if(typ==2):
        xs=my.sql("select xs_star from halo where snap="+str(snap)+";")
        ys=my.sql("select ys_star from halo where snap="+str(snap)+";")
        zs=my.sql("select zs_star from halo where snap="+str(snap)+";")
        rs=my.sql("select rs_star from halo where snap="+str(snap)+";")

    x=[];y=[];z=[];r=[];
    for i in range(len(xs)):
        x.append(xs[i][0])
        y.append(ys[i][0])
        z.append(zs[i][0])
        r.append(rs[i][0])
        
    return x,y,z,r

def get_nhalo(snap):
    nhalo=my.sql("select atime from halo where snap="+str(snap))
#    print len(nhalo)
    return len(nhalo)

def d(nsnap):
    d = my.sql("select dir from sim where id=1")
#    print d[0][0]
#    d = d[0][0]+
    d = "/scratch/psteger/sim/sim_aboley/output_"+str(nsnap).zfill(5)+"/"
    d = "/scratch/psteger/sim/nec_20111220/output_"+str(nsnap).zfill(5)+"/"
    return d
