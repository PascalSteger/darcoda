#!/usr/bin/python

import os
import sys
from numpy import *

# read in halos
all_new=[]; all_old=[]; old=[]; new=[];
s="{"
for i in range(270):
	infile = "mt/"+str(i)+"_idx"
	try:
		in1=open(infile,"r")
	except(IOError), e:
		print "no file"
		for j in range(len(old)):
			s = s + str(i) + "." + new[j]+"->"
			s = s + str(i+1)+"."+ new[j]+","
		continue

	old=[]; new=[]
	for line in in1:
		val = line.split()
		old.append(val[0])
		new.append(val[1])
		s = s + str(i) + "." + val[0]+"->"
		s = s + str(i+1)+"."+ val[1]+","
	all_new.append(new)
	all_old.append(old)
	in1.close()
s = s + "0->0}"
print s	



