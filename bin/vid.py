#!/usr/bin/python

import os
import sys
import threading

nproc = 48
os.nice(8)

semaphore = threading.Semaphore(nproc)
def run_command(cmd):
    with semaphore:
        os.system(cmd)

os.system("mkdir vid")
center = open("mt/pos_sss","r")
i=0
for line in center:
    i=i+1
    val = line.split()
    x = float(val[0])
    y = float(val[1])
    z = float(val[2])
    r = float(val[3])
    r = 0.001

    # no following, simply 0.5, box of 10kpc
    #x = 0.5; y=0.5; z=0.5; r=0.01

    # for dm density:
    threading.Thread(target=run_command, args=("part2map -dir z -nx 512 -ny 512 -inp output_"+str(i).zfill(5)+" -out vid/map_d_"+str(i).zfill(5)+" -xmi "+str(x-r)+" -xma "+str(x+r)+" -ymi "+str(y-r)+" -yma "+str(y+r)+" -zmi "+str(z-r)+" -zma "+str(z+r), )).start()

    #for star density
    threading.Thread(target=run_command, args=("part2map -str true -dir z -nx 512 -ny 512 -inp output_"+str(i).zfill(5)+" -out vid/map_s_"+str(i).zfill(5)+" -xmi "+str(x-r)+" -xma "+str(x+r)+" -ymi "+str(y-r)+" -yma "+str(y+r)+" -zmi "+str(z-r)+" -zma "+str(z+r), )).start()

    # for gas density
    threading.Thread(target=run_command, args=("amr2map -dir z -nx 512 -ny 512 -lma 15 -typ 1 -inp output_"+str(i).zfill(5)+" -out vid/map_g_"+str(i).zfill(5)+" -xmi "+str(x-r)+" -xma "+str(x+r)+" -ymi "+str(y-r)+" -yma "+str(y+r)+" -zmi "+str(z-r)+" -zma "+str(z+r), )).start()

    # for spherical dm density
    threading.Thread(target=run_command, args=("get_sphere_dm_8 -inp output_"+str(i).zfill(5)+"/ -xc "+str(x)+" -yc "+str(y)+" -zc "+str(z)+" -r "+str(r)+" > vid/part_"+str(i).zfill(5)+" && "+"get_prof.py vid/part_"+str(i).zfill(5)+" "+str(i)+" "+str(x)+" "+str(y)+" "+str(z)+" "+str(r)+">vid/prof_"+str(i).zfill(5), )).start()

# fill non-existing snapshots with previous snapshots
#os.system("show_existing.py")

#for i in range(270):
#    s = str(i+1).zfill(5)
#    threading.Thread(target=run_command, args=("map2img.py vid/map_"+s+" -o vid/map_"+s+".png -l -c jet", )).start()
