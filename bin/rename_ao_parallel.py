#!/usr/bin/python

import os

for i in range(270):
    cdc = "cd output_"+str(i+1).zfill(5)+"/ && "
    os.system( cdc + "mv *_halos halos" )
    os.system( cdc + "mv *_centres centres" )
    os.system( cdc + "mv *_particles particles" )
    os.system( cdc + "mv *_profiles profiles" )
    os.system( cdc + "mv *_substructure substructure" )
    os.system( cdc + "rm halo" )
