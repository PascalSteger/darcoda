#!/usr/bin/python

import os
import sys
from numpy import *

nstart=int(sys.argv[1])
nstop =int(sys.argv[2])

# read in halos
all_new=[]; all_old=[]; old=[]; new=[];
s="{"
for i in range(nstop-nstart+1):
	nc = nstart+i
	infile = "mt/"+str(nc)+"_idx"
#	print infile
	try:
		in1=open(infile,"r")
	except(IOError), e:
#		print "no file:",nc
		all_new.append(new)
		all_old.append(new)
		continue

	old=[]; new=[]
	for line in in1:
		val = line.split()
		old.append(int(val[0]))
		new.append(int(val[1]))
	all_new.append(new)
	all_old.append(old)
	in1.close()

tb_id=[]; tb_x=[]; tb_y=[]; tb_z=[]; tb_r=[]; tb_m=[]
# get halo of interest in latest snapshot
tb_now = 0
tb_id.append(tb_now)
# trace it backwards
for i in range(nstop-nstart+1):
	# get prospective progenitors
	j = nstop-nstart-i
#	print nstart+j
	prog = all_new[j]
#	print "prog",prog
#	print "len(all_old)",len(all_old[j]),", len(prog)",len(prog)
	pos_prog=[]
	for k in range(len(prog)):
		if(int(prog[k])==tb_now):
			pos_prog.append(all_old[j][k])

	infile = "output_"+str(nstart+j).zfill(5)+"/halof"
	try:
		in1=open(infile,"r")
	except(IOError), e:
#		print "backtrace: no file:",j
		tb_x.append(0)
		tb_y.append(0)
		tb_z.append(0)
		tb_m.append(-1.0)
		tb_r.append(0)
		continue
	
	# read x/y/z/m/rvir
	id=[];x=[];y=[];z=[];m=[];r=[]
	for line in in1:
		val = line.split()
		id.append(int(val[0]))
		m.append(float(val[1]))
		x.append(float(val[2]))
		y.append(float(val[3]))
		z.append(float(val[4]))
		r.append(float(val[5]))
	in1.close()
	#print len(x)

	# take most massive one
	mmax = -1.0; imax=0
#	print "pos_prog",pos_prog
#	print "len(m)",len(m)
	for k in range(len(pos_prog)):
		if(m[pos_prog[k]]>mmax):
#			print "if"
			mmax=m[pos_prog[k]]
			imax=id[pos_prog[k]]-1
	tb_now = imax
	#print "imax",imax
	# for later analysis: copy respective dm_*.dat file to dm.dat
	dmdat = "/scratch/psteger/sim_aboley/output_"+str(nstart+j).zfill(5)+"/dm/dm.dat"
	dmimaxdat = "/scratch/psteger/sim_aboley/output_"+str(nstart+j).zfill(5)+"/dm/dm_"+str(imax+1)+".dat"
	os.system("rm "+dmdat)
	os.system("ln -s "+dmimaxdat + " " + dmdat)
	
	# fill in x/y/z/m/rvir
	if(len(x)==0):
		tb_x.append(0)
		tb_y.append(0)
		tb_z.append(0)
		tb_m.append(-1.0)
		tb_r.append(0)
	else:
		tb_x.append(x[imax])
		tb_y.append(y[imax])
		tb_z.append(z[imax])
		tb_m.append(m[imax])
		tb_r.append(r[imax])

# if all is ok (all halos exist, stop here)
# TODO: handle missing snapshots
#print len(tb_x)
for i in range(nstop-nstart+1):
	j=nstop-nstart-i
#	print j
#	print tb_x[j],tb_y[j],tb_z[j],tb_r[j]

#print " fill in interpolated x/y/z/m/rvir"
nc = nstop-nstart
#print "nc ",nc, tb_x[nc], tb_m[nc]
if(tb_m[nc]<0):
	tb_x[nc]=0.5; tb_y[nc]=0.5; tb_z[nc]=0.5; tb_r[nc]=0.0; tb_m[nc]=1;
for i in range(nstop-nstart+1):
	#print tb_m[i]
	if(tb_m[i]<0):
#		print "no info at i=",i
		stop = 0
		for j in range(nstop-nstart+1):
#			print "i+j<len: ",i,j,len(tb_m)
			if(not tb_m[i+j]<0):
				stop = j+1
				break
		for j in range(stop):
			#print nci,ncj,stop,len(tb_x)
			tb_x[i+j]=tb_x[i-1]-(tb_x[i-1]-tb_x[i-1+stop])*(j+1)/stop
			tb_y[i+j]=tb_y[i-1]-(tb_y[i-1]-tb_y[i-1+stop])*(j+1)/stop
			tb_z[i+j]=tb_z[i-1]-(tb_z[i-1]-tb_y[i-1+stop])*(j+1)/stop
			tb_m[i+j]=tb_m[i-1]-(tb_m[i-1]-tb_m[i-1+stop])*(j+1)/stop
			tb_r[i+j]=tb_r[i-1]-(tb_r[i-1]-tb_r[i-1+stop])*(j+1)/stop
		i = i+stop-1

for i in range(nstop-nstart+1):
	j = nstop-nstart-i
	print tb_x[j],tb_y[j],tb_z[j],tb_r[j]
