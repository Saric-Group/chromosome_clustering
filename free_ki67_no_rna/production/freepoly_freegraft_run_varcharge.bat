#!/bin/bash
#
#SBATCH -J FG_run
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=4G 
#SBATCH --time=170:00:00 
#SBATCH --mail-user=valeriosorichetti@gmail.com 
#SBATCH --mail-type=END
#SBATCH --export=NONE 

module load openmpi/4.1.4

mkdir dump dumplin restart
printf $RANDOM > lseed.dat
printf $RANDOM > vseed.dat

time mpirun /nfs/scistore15/saricgrp/vsoriche/lammps_23062022/src/lmp_mpi -in freepoly_freegraft_run_varcharge.in 
