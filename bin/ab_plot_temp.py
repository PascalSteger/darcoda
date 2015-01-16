#!/usr/bin/python
'''plot temperature'''

import os

# file-dependent variables.
istart=270
istop =270

directory="MAPS/"
prefix1="p_diry"
prefix2="g_diry"
prefix3="tk_diry"

RM=4
mp=1.67e-24
kb=1.3801e-16
USECENTERING=0

for i in xrange(istart,istop):
	number=repr(i).zfill(5)
	filename1=directory+prefix1+number+".dat"         
	filename2=directory+prefix2+number+".dat"         
	filename3=directory+prefix3+number+".png"
	command = "./ab_map_tk.py "+filename1+" "+filename2+" "+number+" "+filename3+" 1"
	print command 
	os.system(command)
        

