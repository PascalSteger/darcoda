#!/usr/bin/python
# sample script: do 10 times sleeping in parallel

import os
import sys
import threading

nproc = 48
njob  =128

semaphore = threading.Semaphore(nproc)

def run_command(cmd):
    with semaphore:
        os.system(cmd)
                
for d in range(njob):
    threading.Thread(target=run_command, args=("sleep 10", )).start()
