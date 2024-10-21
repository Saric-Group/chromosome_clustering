# chromosome_clustering
LAMMPS and Python codes for simulations of Ki-67 and RNA mediated chromosome clustering.

The simulation model is described in detail in:

_A liquid-like coat mediates chromosome clustering during mitotic exit_
Alberto Hernandez-Armendariz, Valerio Sorichetti, Yuki Hayashi, Zuzana Koskova, Andreas Brunner, Jan Ellenberg, Anđela Šarić, Sara Cuylen-Haering, [Molecular Cell, Volume 84, Issue 17, (2024)](https://doi.org/10.1016/j.molcel.2024.07.022) 

## What this repository contains

LAMMPS codes to simulate Ki-67 polymers interacting with RNA polymers in three different geometries:

- Half cylinder grafted with Ki-67 polymers with and wihtout free RNA polymers.
- Parallel surfaces grafted with Ki-67 polymers with and without free RNA polymers.
- Free Ki-67 polymers with and without free RNA polymers.

## How to run the code

Each folder contains an "equilibration" sub-folder and a "production" sub-folder. The equilibration sub-folder contains a python script that can be used to generate the initial configuration. The initial configuration can then be equilibrated by running the LAMMPS script contained in the equilibration sub-folder. During equilibration, all charges are set to zero.

Once equilibration is completed, the restart file produced can be used to perform the production run. During production, the charge distribution of the Ki-67 molecules is gradually changed using the sigmoidal function described in bioRxiv 2023.11.23.568026. 

The code was originally wrote for and run with LAMMPS stable version 23.06.2022

## How to cite this material

When citing this material, reference the original work and the following DOI: [10.5281/zenodo.12706119](https://doi.org/10.5281/zenodo.13961010) (v2).
