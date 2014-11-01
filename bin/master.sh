#!/bin/bash
#PBS -M psteger@phys.ethz.ch
#PBS -N torque_run
#PBS -l walltime=08:00:00
#PBS -l nodes=1:ppn=32
#PBS -l mem=16gb
#PBS -l vmem=48gb
cd $PBS_O_WORKDIR

export PATH=$PATH:/scratch/psteger/tools/bin/:/scratch/psteger/software/bin/

/scratch/psteger/tools/bin/master.py $1