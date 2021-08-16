#!/bin/bash
#SBATCH --mem 60000
#SBATCH -J num-evol
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --ntasks-per-core=1
Rscript main-kirk-2.R

