#!/usr/bin/python2.5

import sys
import os

if len(sys.argv)!=3:
	print "Incorrect number of arguments"

f=open(sys.argv[1],"r")
g=open(sys.argv[2],"r")

find_these_ids=[]
for line in f:
	if line[0]=="#" or line[1]=="#": continue
	val=line.split()
	find_these_ids.append(int(val[4]))

f.close()

n=len(find_these_ids)
vmax=max(find_these_ids)
vmin=min(find_these_ids)
print "#",n
for line in g:
	if line[0]=="#" or line[1]=="#": continue
	val=line.split()
	id=int(val[4])
	if vmin <= id <= vmax:
		i=find_these_ids.count(id)
		if i>0:
			i=find_these_ids.index(id)
			print line.replace("\n","")
			find_these_ids.pop(i)
g.close()
