#!/usr/bin/env python3

#Submit_Lisa_Job.py
#
#20 October 2014, Hamish Silverwood
#
# This programme copies the gravimage code to a holding area, produces an appopriate
# PDB file, and then submits the PDB file to the Lisa cluster. The PDB file copies
# code from the holding area, rather than the main programs folder. Multiple copies
# of the code are held in the holding area under different time stamps, each with
# different PDB files.

from numpy import *
import os, datetime, shutil
from subprocess import call
import pdb

nodes=1
cores='16'
ppn=1
walltime='03:00:00:00'

gravimage_path = os.path.abspath('../')
holding_stack_path = gravimage_path + '/holding_stack/'
investigation = 'simplenu'

#NOT YET IMPLEMENTED
resume_run = False
resume_ts = '201507141033'

#Check for holding area folder, create one if none exists
if os.path.isdir(holding_stack_path) != True:
    os.mkdir(holding_stack_path)

#Generate holding number
holding_number = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
#os.mkdir(holding_stack_path + holding_number)

#Copy code to holding area with timestamp
if resume_run:
    #shutil.copytree(gravimage_path +'/DTsimplenu/0/' + resume_ts + '/programs', holding_stack_path + holding_number + '/programs/')
    shutil.copytree(gravimage_path +'/DTsimplenu/0/' + resume_ts, holding_stack_path + holding_number + '/')
else:
    shutil.copytree(gravimage_path + '/programs', holding_stack_path + holding_number + '/programs/')

#Remove __pycache__ from the holding area
shutil.rmtree(holding_stack_path + holding_number + '/programs/__pycache__', ignore_errors=True)

#Create PBS file
pbs_filename = holding_stack_path + holding_number + '/programs/PBS_LisaSubmit_' + holding_number + '.pbs'
pbs_file = open(pbs_filename, 'w')
pbs_file.writelines('#!/bin/bash'+'\n')
pbs_file.writelines('#PBS -S /bin/bash'+'\n')
if cores == 'any':
    pbs_file.writelines('#PBS -lnodes=' + str(nodes) + ':ppn=' + str(ppn) + ',walltime=' + walltime+'\n')
else:
    pbs_file.writelines('#PBS -lnodes=' + str(nodes) + ':cores' + str(cores) + ':ppn=' + str(ppn) + ',walltime=' + walltime+'\n')
pbs_file.writelines('# Copying program files to scratch'+'\n')
pbs_file.writelines('#module load openmpi/intel'+'\n')
pbs_file.writelines('export RUNDIR=$TMPDIR'+'\n')
pbs_file.writelines('cd $TMPDIR'+'\n')
pbs_file.writelines('mkdir -p $TMPDIR/darcoda/gravimage/programs'+'\n')
pbs_file.writelines('mkdir $TMPDIR/darcoda/gravimage/programs/HoldingNumberWas_' + holding_number + '\n')
pbs_file.writelines('cp -r $HOME/LoDaM/darcoda/gravimage/holding_stack/' + str(holding_number) +'/programs/* $TMPDIR/darcoda/gravimage/programs/'+'\n')
pbs_file.writelines('cd $TMPDIR/darcoda/gravimage/programs'+'\n')
pbs_file.writelines('# Calculate run time for gravimage, less than wall time to allow for data to be'+'\n')
pbs_file.writelines('# copied back, allow [transft] seconds for transfer.'+'\n')
pbs_file.writelines('echo PBS_WALLTIME = $PBS_WALLTIME'+'\n')
pbs_file.writelines('transft=180'+'\n')
pbs_file.writelines('echo Transfer time = $transft'+'\n')
pbs_file.writelines('runtime=$(expr $PBS_WALLTIME - $transft)'+'\n')
pbs_file.writelines('echo gravimage runtime = $runtime'+'\n')

pbs_file.writelines('echo Python3 version used:' + '\n')
pbs_file.writelines('which python3' + '\n')

pbs_file.writelines('echo Loading LISA version' + '\n')
pbs_file.writelines('module load python/3.4.2' + '\n')

pbs_file.writelines('echo Python3 version used:' + '\n')
pbs_file.writelines('which python3' + '\n')

pbs_file.writelines('module load openmpi/intel' + '\n')

pbs_file.writelines('echo Contents of programs folder:' + '\n')
pbs_file.writelines('ls -l -a' + '\n')

pbs_file.writelines('echo Contents of DT folder:' + '\n')
pbs_file.writelines('ls -R -l $TMPDIR/darcoda/gravimage/DT' + '\n')

pbs_file.writelines('python3 gravimage.py --investigation ' + investigation + '& PID=$!; sleep $runtime; kill $PID'+'\n')
#pbs_file.writelines('mpiexec -n 2 python3 gravimage.py --investigation simplenu' + '\n')
#pbs_file.writelines('mpiexec -n 16 python3 gravimage.py --investigation ' + investigation + '& PID=$!; sleep $runtime; kill $PID'+'\n')
pbs_file.writelines('echo gravimage killed, transfering data'+'\n')
pbs_file.writelines('pwd'+'\n')
pbs_file.writelines('ls -l -a $TMPDIR/darcoda/gravimage/DT' + investigation +'/0/*'+'\n')
pbs_file.writelines('cp -r $TMPDIR/darcoda/gravimage/DT' + investigation +'/0/* $HOME/LoDaM/darcoda/gravimage/DT'+ investigation +'/0/'+'\n')
pbs_file.writelines('echo Data transfered, job finished'+'\n')
pbs_file.close()

#Submit PDB file
#pdb.set_trace()
os.chdir(holding_stack_path + holding_number + '/programs/')
os.system('qsub ' + pbs_filename)
