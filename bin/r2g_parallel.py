#!/usr/bin/python

import os
import sys
import threading
semaphore = threading.Semaphore(48)

def run_command(cmd):
    with semaphore:
        os.system(cmd)
                
for d in range(270):
    threading.Thread(target=run_command, args=("cd output_"+str(d).zfill(5)+" && r2g -i ../output_"+str(d).zfill(5), )).start()
