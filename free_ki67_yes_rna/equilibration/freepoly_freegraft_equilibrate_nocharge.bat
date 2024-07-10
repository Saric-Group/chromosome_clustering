#!/bin/bash
#
#SBATCH -J FPFG_eq
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=2G 
#SBATCH --time=170:00:00 
#SBATCH --mail-user=valeriosorichetti@gmail.com 
#SBATCH --mail-type=END
#SBATCH --export=NONE 

module load lammps/20220623b

mkdir dump dumplin restart
rm -f freepoly*lammpsdata
python create_freepoly_freegraft_rw.py
cp freepoly*lammpsdata data
printf $RANDOM > lseed.dat
printf $RANDOM > vseed.dat

time mpirun lmp_mpi -in freepoly_freegraft_equilibrate_nocharge.in 
