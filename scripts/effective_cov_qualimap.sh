#!/bin/bash
#SBATCH --output=/mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/qualimap/slurm-%j.out
#SBATCH --error=/mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/qualimap/slurm-%j.err

# This script calculates the effective coverage per sample from the Qualimap output.
# The effective coverage is estimate has -ln(1 - f_cov), where f_cov is the fraction of genomic positions coveref by at least one read.

# Usage of the script: bash effective_cov_qualimap.sh <input_qualimap_directory>

# Define the input directory
input_dir=${1}

# find all the "genome_fraction_coverage.txt" files 
qualimap_list=$(find ${1} -name "genome_fraction_coverage.txt") 

# create the output table header
echo -e "sample\teffective_coverage" > ${input_dir}/effective_coverage.txt

# iterate over the files
for file in ${qualimap_list}; do

    # get the sample name
    sample=$(echo "${file}" | sed 's/.*\/\(.*\)_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner_qualimap\/raw_data_qualimapReport\/genome_fraction_coverage.txt/\1/')

    # get the fraction of genomic positions covered by at least one read
    f_cov=$(echo "$(sed -n '2p' ${file} | cut -f2) / 100" | bc -l)
    
    # calculate the effective coverage
    effective_cov=$(echo "-l(1 - $f_cov)" | bc -l)
    
    # print the sample name and the effective coverage
    echo -e "${sample}\t${effective_cov}" >> ${input_dir}/effective_coverage.txt

done
