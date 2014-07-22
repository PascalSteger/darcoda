#!/bin/bash

#PBS -N gravlite
#PBS -o gravlite.out
#PBS -e gravlite.err
#PBS -l walltime=00:03:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=2gb
#PBS -l vmem=8gb
#PBS -m abe
#PBS -M psteger@phys.ethz.ch
## cd $PBS_O_WORKDIR

mpirun -output-filename chains/example ./gravlite_pymulti.py
