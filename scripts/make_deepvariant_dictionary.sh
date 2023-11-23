#!/bin/bash

# In this script I will create a dictionary for the python script (make_deepvariant_scripts.py) that writes the scripts for the deepvariant pipeline.
# The dictionary will have the following format:
# sample_dict:
#  sample1:
#   ['path/to/bam1', 'file1.bam', 'male']
# sample2:
#   ['path/to/bam2', 'file2.bam', 'female']

## CAUTION: I need to be carefull with the sample list and hidden characters (you can see them with sed -n l filename.txt). 
	# I need to remove the hidden characters from the sample list before running this script, in this case I used sed -i 's/male\t/male/g' fastq_samples_list.txt. 


# The script will take as input a list of fastq files and a folder with the bam files. It will then search for the bam files that correspond to the fastq files and create the dictionary.

# Specify the paths to the fastq list and bam folder.
FASTQ_LIST='/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/FASTQ_files/novogene_lp_sept2023/fastq_samples_list.txt'
BAM_FOLDER='/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams'
DICTIONARY='all_bams_calling.yml'


while read line ; do
	# get the sample name from the fastq list:
	SAMPLE_NAMEFQ=($(echo "${line}" | cut -f1))

	#get the final sample name:
	SAMPLE_NAME=($(echo "${line}" | cut -f2))

	# get the individuals' sex:
	SEX=($(echo "${line}" | cut -f3))

	# get the bam file name and path:
	BAM_INFO=$(find "${BAM_FOLDER}" -type f -name "${SAMPLE_NAMEFQ}_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.bam")

	# get the bam file name:
	BAM=($(basename "${BAM_INFO}"))

	# get the bam file path:
	BAM_PATH=($(dirname "${BAM_INFO}"))

	# Now we will arrange the information in the dictionary format:
	dict_entry="${SAMPLE_NAME}:\n   ['${BAM_PATH}', '${BAM}', '${SEX}']"

	# Now we will print the dictionary entry:
	echo -e "${dict_entry}" >> ${DICTIONARY}

	echo "${BAM}"
done < <(sed 's/_EKDN[0-9]*-1A_[A-Z0-9]*_L[0-9]*_[0-9]*.fastp.fq.gz//g' ${FASTQ_LIST} | uniq)
