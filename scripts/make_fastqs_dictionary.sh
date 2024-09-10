#!/bin/bash

# In this script I will create a dictionary for Enrico's python package (congenomics_fastq_align) to generate the aligments scripts with the python script (run_alignment.py) 
# The dictionary will have the following format:
# ---
# sample_dict:
#  sample1:
#    fastq1:
#      ['path/to/fastq1', 'file_1.fastq', 'file_2.fastq']
#    fastq2:
#      ['path/to/fastq2', 'file2_1.fastq', 'file2_2.fastq']
#  sample2:
#    ['path/to/fastq2', 'file3_1.fastq', 'file3_2.fastq']

# The script will take as input a list of fastq files that contains the path to each fastq and the prefix of each pair of files (without the suffix) in the first column, and the final sample name in the second column, tab-separated. 

# Specify the paths to the fastq list and output dictionary folder
FASTQ_LIST='/path/to/fastq_list.txt'
DICTIONARY='/path/to/output/test.yml'

# Print the first line with "sample_dict:"
echo "---\nsample_dict:" > "${DICTIONARY}"

# Initialize an empty array to store processed sample names
declare -a PROCESSED_SAMPLES=()

# Read the FASTQ list file line by line
while read line ; do
	# get the fastq ID from the fastq list:
	FQ=($(echo "${line}" | cut -f1))

	# get the fastq code:
	CODEFQ=($(basename "${FQ}"))

	# get the fastq uniq code (the lane where it was sequenced). This variable definition might change a lot depending on the format of the fastq name:
    	CODEFQ_ID=($(echo "${CODEFQ}" | rev | cut -f2- -d'_' | rev))

	# get the final sample name:
	SAMPLE_NAME=($(echo "${line}" | cut -f2))

	# generate the fastp variables:
	FASTP1=($(echo "${CODEFQ}_1.fastp.fastq.gz"))
	FASTP2=($(echo "${CODEFQ}_2.fastp.fastq.gz"))

	# get the input faspt fastq path:
	FASTP_PATH=($(dirname "${FQ}"))
	FASTP_PATH="${FASTP_PATH}/fastp"

	# Check if the sample name has been processed before and add the fastq right after it otherwise
    if [[ " ${PROCESSED_SAMPLES[@]} " =~ " ${SAMPLE_NAME} " ]]; then
        # If the sample name has been processed before, omit the sample name
        echo -e "    ${CODEFQ_ID}:\n      ['${FASTP_PATH}', '${FASTP1}', '${FASTP2}']" >> "${DICTIONARY}"
    else
        # If the sample name has not been processed before, print the sample name and the rest of the variables
        echo -e "  ${SAMPLE_NAME}:\n    ${CODEFQ_ID}:\n      ['${FASTP_PATH}', '${FASTP1}', '${FASTP2}']" >> "${DICTIONARY}"
        # Add the sample name to the array of processed sample names
        PROCESSED_SAMPLES+=("${SAMPLE_NAME}")
    fi

done <  ${FASTQ_LIST}

echo -e "\n\nreference_fasta: \"/mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.scaffolds.revcomp.scaffolds.fa\"" >> "${DICTIONARY}"
echo -e "\noutput_folder: \"/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/pool_epil_all\"" >> "${DICTIONARY}"
echo -e "\nthreads: 20" >> "${DICTIONARY}"
echo -e "\nalignment_name: \"mLynPar1.2_ref\"" >> "${DICTIONARY}"
echo -e "\ncall_bwa: \"bwa\"" >> "${DICTIONARY}"
echo -e "\ncall_samtools: \"samtools\"" >> "${DICTIONARY}"
echo -e "\ncall_picard: \"java -jar \$EBROOTPICARD/picard.jar\"" >> "${DICTIONARY}"
echo -e "\ncall_gatk: \"java -jar \$EBROOTGATK/GenomeAnalysisTK.jar\"" >> "${DICTIONARY}"
