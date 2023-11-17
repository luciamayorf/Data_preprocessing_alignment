#!/bin/bash

# This script performs MultiQC on the different results of a data quality control step (with FastQC, Qualimap or samtools).

#SBATCH --output=logs/multiqc/sample1_fastq_pair_1.out
#SBATCH --error=logs/multiqc/sample1_fastq_pair_1.err

# First argument of the script defines the directory containing all the results we want to use as input of MultiQC
INPUT_QC_DIR=$1

# loading the module
module load multiqc

multiqc $1/* -o $1/