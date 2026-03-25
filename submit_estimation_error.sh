#!/bin/bash -l
#PBS -N hdGLS_sim
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=1:mem=50G
#PBS -j oe

module purge
module add tools/dev
module load R/4.2.1-foss-2022a

cd $PBS_O_WORKDIR

R --file=$PBS_O_WORKDIR/HPC_estimation_error_simulation.R
