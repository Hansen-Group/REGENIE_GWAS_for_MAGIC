#!/bin/bash
# Set SLURM options:
#SBATCH --nodes=1
#SBATCH --mem=256G
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=05:00:00
#SBATCH --partition=standardqueue
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load --auto R/4.3.3

COHORT="cohort_name"   
Rscript extract_cov.R $COHORT 
