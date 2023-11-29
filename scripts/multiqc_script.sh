#!/bin/bash

# This script performs MultiQC on the different results of a data quality control step (with FastQC, Qualimap or samtools).

#SBATCH --output=/mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/multiqc/slurm-%j.out
#SBATCH --error=/mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/multiqc/slurm-%j.err

# First argument of the script defines the directory containing all the results we want to use as input of MultiQC
INPUT_QC_DIR=$1

# loading the module
module load multiqc

multiqc $1/* -o $1/
