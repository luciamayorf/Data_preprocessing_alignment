#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/fastqc/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/fastqc/slurm-%j.out

# This script performs a fastqc of a given list of files.

# Example of the fastqc command:
#   fastqc -o /path/to/pair1/fastqc -t 6 /path/to/pair1/pair1_r1.fastq.gz /path/to/pair1/pair1_r2.fastq.gz

# Usage of the script simultaneously for several bam files:
#    for input_fastq in $(cat fastq_list.txt); do 
#        job_id=(sbatch -c 6 --mem=4GB -t 03:00:00 fastqs_fastQC.sh <${input_fastq}> | awk '{print $4}')
#        echo "${job_id} ${input_bam}" >> /logs/fastqc/job_ids_fastqc.txt
#    done

# Load the fastqc module:
module load cesga/2020 fastqc/0.11.9
module load multiqc

# Define the input bam file:
fastq=${1}

# Define the output directory:
output_directory=$(dirname ${fastq})

# Create the output directory if it doesn't exist:
mkdir -p ${output_directory}/fastqc

# Run FastQC in the pair of reads:
fastqc -o ${output_directory}/fastqc -t 6 ${fastq}
