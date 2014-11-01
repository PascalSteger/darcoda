#!/usr/bin/python

import os
import sys
import threading
semaphore = threading.Semaphore(48)

def run_command(cmd):
    with semaphore:
        os.system(cmd)
                
for root,dirs,files in os.walk('.'):
    for d in dirs:
        print d

        os.system('mkdir '+d+'/amr')
        os.system('mv '+d+'/amr_* '+d+'/amr')

        os.system('mkdir '+d+'/grav')
        os.system('mv '+d+'/grav_* '+d+'/grav')
        
        os.system('mkdir '+d+'/hydro')
        os.system('mv '+d+'/hydro_* '+d+'/hydro')

        os.system('mkdir '+d+'/part')
        os.system('mv '+d+'/part_* '+d+'/part')

#        threading.Thread(target=run_command, args=("mkdir "+d+"/r2g", )).start()
        
#    for f in files:
#        print os.path.splitext(os.path.join(root,f))
