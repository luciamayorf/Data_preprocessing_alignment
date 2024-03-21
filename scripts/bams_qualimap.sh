#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/qualimap/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/qualimap/slurm-%j.out

# This script performs a qualimap quality QC of a given list of BAM files.

# Example of the qualimap command:
#   qualimap bamqc sample1.bam --java-mem-size=19G -outfile alignment1_qualimap.html -outformat html -outdir /path/to/bams/qualimap_output/alignment1_qualimap

# Usage of the script simultaneously for several bam files:
#    for input_bam in $(ls /path/to/bams/*.bam); do 
#        job_id=(sbatch -c 10 --mem=20GB -t 06:00:00 bams_qualimap.sh <${input_bam}> | awk '{print $4}')
#        echo "${job_id} ${input_bam}" >> /logs/qualimap/job_ids_qualimap.txt
#    done

# Load the qualimap module:
module load qualimap

# Define the input bam file:
bam=${1}

# Define the output directory:
output_directory=$(dirname ${bam})

# Define the name of the bam file:
bam_name=$(basename ${bam} .bam)

# Create the output directory if it doesn't exist:
mkdir -p ${output_directory}/qualimap_output

# Run qualimap:
qualimap bamqc -bam ${bam} --skip-duplicated --java-mem-size=19G -outfile ${bam_name}_qualimap.html -outformat html -outdir ${output_directory}/qualimap_output/${bam_name}_qualimap
