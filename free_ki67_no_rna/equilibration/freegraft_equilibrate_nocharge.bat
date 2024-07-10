#!/bin/bash
#
#SBATCH -J freegraft_eq
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=4G 
#SBATCH --time=170:00:00 
#SBATCH --mail-user=valeriosorichetti@gmail.com 
#SBATCH --mail-type=END
#SBATCH --export=NONE 

module load lammps/20220623b

mkdir dump dumplin restart
rm -f freepoly*lammpsdata
python create_freegraft_rw.py
cp freepoly*lammpsdata data
printf $RANDOM > lseed.dat
printf $RANDOM > vseed.dat

time mpirun lmp_mpi -in freegraft_equilibrate_nocharge.in 
