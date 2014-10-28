#!/usr/bin/env python3

#Submit_Lisa_Job.py
#
#20 October 2014, Hamish Silverwood
#
# This programme copies the gravlite code to a holding area, produces an appopriate
# PDB file, and then submits the PDB file to the Lisa cluster. The PDB file copies
# code from the holding area, rather than the main programs folder. Multiple copies
# of the code are held in the holding area under different time stamps, each with
# different PDB files.

from numpy import *
import os, datetime, shutil
from subprocess import call
import pdb

nodes=1
cores=16
ppn=16
walltime='00:00:15:00'

gravlite_path = os.path.abspath('../')
holding_stack_path = gravlite_path + '/holding_stack/'

#Check for holding area folder, create one if none exists
if os.path.isdir(holding_stack_path) != True:
    os.mkdir(holding_stack_path)

#Generate holding number
holding_number = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
#os.mkdir(holding_stack_path + holding_number)

#Copy code to holding area with timestamp
shutil.copytree(gravlite_path + '/programs', holding_stack_path + holding_number + '/programs/')

#Create PBS file
pbs_filename = holding_stack_path + holding_number + '/programs/PBS_LisaSubmit_' + holding_number + '.pbs'
pbs_file = open(pbs_filename, 'w')
pbs_file.writelines('#!/bin/bash'+'\n')
pbs_file.writelines('#PBS -S /bin/bash'+'\n')
pbs_file.writelines('#PBS -lnodes=' + str(nodes) + ':cores' + str(cores) + ':ppn=' + str(ppn) + ',walltime=' + walltime+'\n')
pbs_file.writelines('# Copying program files to scratch'+'\n')
pbs_file.writelines('#module load openmpi/intel'+'\n')
pbs_file.writelines('export RUNDIR=$TMPDIR'+'\n')
pbs_file.writelines('cd $TMPDIR'+'\n')
pbs_file.writelines('mkdir -p $TMPDIR/darcoda/gravlite/programs'+'\n')
pbs_file.writelines('mkdir $TMPDIR/darcoda/gravlite/programs/HoldingNumberWas_' + holding_number)
pbs_file.writelines('cp -r $HOME/LoDaM/darcoda/gravlite/holding_stack/' + str(holding_number) +'/programs/* $TMPDIR/darcoda/gravlite/programs/'+'\n')
pbs_file.writelines('cd $TMPDIR/darcoda/gravlite/programs'+'\n')
pbs_file.writelines('# Calculate run time for gravlite, less than wall time to allow for data to be'+'\n')
pbs_file.writelines('# copied back, allow [transft] seconds for transfer.'+'\n')
pbs_file.writelines('echo PBS_WALLTIME = $PBS_WALLTIME'+'\n')
pbs_file.writelines('transft=180'+'\n')
pbs_file.writelines('echo Transfer time = $transft'+'\n')
pbs_file.writelines('runtime=$(expr $PBS_WALLTIME - $transft)'+'\n')
pbs_file.writelines('echo gravlite runtime = $runtime'+'\n')
pbs_file.writelines('python3 gravlite.py& PID=$!; sleep $runtime; kill $PID'+'\n')
pbs_file.writelines('echo gravlite killed, transfering data'+'\n')
pbs_file.writelines('cp -r $TMPDIR/darcoda/gravlite/DTdiscmock/0/* $HOME/LoDaM/darcoda/gravlite/DTdiscmock/0/'+'\n')
pbs_file.writelines('echo Data transfered, job finished'+'\n')
pbs_file.close()

#Submit PDB file
pdb.set_trace()
os.chdir(holding_stack_path + holding_number + '/programs/')
os.system('qsub ' + pbs_filename)
