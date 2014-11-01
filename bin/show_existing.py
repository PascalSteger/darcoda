#!/usr/bin/python
# handling missing snapshots: copy old maps to missing ones

import os

for i in range(270):
	i5 = str(i).zfill(5)
        try:
            the_file = open("vid/map_"+i5, "r")
        except(IOError), e:
            print "Unable to open the file ", i, ", copying last snapshot\n", e
            os.system("cp vid/map_"+str(i-1).zfill(5)+" vid/map_"+i5)
