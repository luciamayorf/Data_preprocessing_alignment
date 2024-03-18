#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/fastp/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/fastp/slurm-%j.out

# This script performs a fastp trimming of a given list of FASTQ files.

# Example of the fastp command:
#   fastp \
# -i path/to/pair1/pair1_r1.fastq.gz -I path/to/pair1/pair1_r2.fastq.gz \
# -o path/to/pair1/fastp/pair1_r1.fastp.fastq.gz -O path/to/pair1/fastp/pair1_r2.fastp.fastq.gz \
# -h path/to/pair1/fastp/pair1_fastp.html -j path/to/pair1/fastp/pair1_fastp.json \
# --unpaired1 path/to/pair1/fastp/pair1_unpaired.fastq.gz --unpaired2 path/to/pair1/fastp/pair1_unpaired.fastq.gz \
# --failed_out path/to/pair1/fastp/pair1_failed.fastq.gz \
# --dont_overwrite \
# --trim_poly_g \
# --trim_poly_x \
# --correction \
# --detect_adapter_for_pe \
# --thread 6

# Usage of the script simultaneously for several bam files:
#    for input_fastq in $(cat fastq_list.txt); do 
#        job_id=(sbatch -c 6 --mem=4GB -t 03:00:00 fastp_trimming.sh <${input_fastq}> | awk '{print $4}')
#        echo "${job_id} ${input_bam}" >> /logs/fastp/job_ids_fastp.txt
#    done

# Load the fastp module:
module load cesga/2020 gcccore/system fastp/0.22.0

# Define the input bam file:
fastq=${1}

# Define the fastq filename:
fastq_basename=$(basename ${fastq})

# Define the fastq directory:
fastq_directory=$(dirname ${fastq})

# Create the fastp output directory if it doesn't exist:
mkdir -p ${fastq_directory}/fastp

# Run fastp in the pair of reads:
fastp \
    -i ${fastq}_1.fastq.gz -I ${fastq}_2.fastq.gz \
    -o ${fastq_directory}/fastp/${fastq_basename}_1.fastp.fastq.gz -O ${fastq_directory}/fastp/${fastq_basename}_2.fastp.fastq.gz \
    -h ${fastq_directory}/fastp/${fastq_basename}_fastp.html -j ${fastq_directory}/fastp/${fastq_basename}_fastp.json \
    --unpaired1 ${fastq_directory}/fastp/${fastq_basename}_unpaired.fastq.gz --unpaired2 ${fastq_directory}/fastp/${fastq_basename}_unpaired.fastq.gz \
    --failed_out ${fastq_directory}/fastp/${fastq_basename}_failed.fastq.gz \
    --dont_overwrite \
    --trim_poly_g \
    --trim_poly_x \
    --correction \
    --detect_adapter_for_pe \
    --thread 6