#!/bin/bash
#
#SBATCH -J 2S_eq_nocharge
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=2G 
#SBATCH --time=120:00:00 
#SBATCH --mail-user=valeriosorichetti@gmail.com 
#SBATCH --mail-type=END
#SBATCH --export=NONE 

module load lammps/20220623b

mkdir dump dumplin restart
rm -f surfaces*
python create_2surf_freepoly_uniform_b5.py
cp surfaces* data
printf $RANDOM > lseed.dat
printf $RANDOM > vseed.dat

time mpirun lmp_mpi -in 2surf_freepoly_uniform_b5_equilibrate_nocharge.in 
