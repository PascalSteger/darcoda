#!/usr/bin/python
import fileinput
import os

for i in range(270):
    cdc   = "output_"+str(i+1).zfill(5)+"/"
    if(not os.path.exists(cdc+"halos")):
        print "file does not exist!"
        continue
    halos = open(cdc+"halos","r");
    halo  = open(cdc+"halo_unsorted","w")
    h = 0.719
    for line in halos:
        if(line[0]=="#"):
            continue
        val = line.split()
        con =  str(float(val[9-1])/h)
        con += " "
        con += str(float(val[3-1])/h)
        con += " "
        con += str(float(val[4-1])/h)
        con += " "
        con += str(float(val[5-1])/h)
        con += " "
        con += str(float(val[10-1])/h/1000.)
        con += "\n"
        halo.write(con)
    halo.close();    halos.close()

    os.system("sort -nr "+cdc+"halo_unsorted > "+cdc+"halom")
    halom = open(cdc+"halom","r");
    halo  = open(cdc+"halo","w")
    i = 0
    for line in halom:
        i = i + 1
        halo.write(str(i)+" "+line)
    halom.close();
    halo.close()
