# Data_preprocessing_alignment_v2

This repository is an inspiration from Enrico's [Lynx_demography](https://github.com/Enricobazzi/Lynx_Demography) repository.

In this repository, I keep the scripts used for sequencing data quality control ([fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [multiQC](https://multiqc.info/docs/)), trimming ([fastp](https://github.com/OpenGene/fastp)), alignment (Enrico's [congenomics_fastq_align](https://github.com/Enricobazzi/congenomics_fastq_align) python package) and quality control of the aligment ([qualimap](http://qualimap.conesalab.org/doc_html/analysis.html)).

** The scripts are only tested with one sample of the ones sequenced in previous projects, but they seem to work.

## 1. Raw data quality control

First, we do a first FastQC analysis of the FASTQ files received with FastQC and multiQC to summarize the results. I run the script [fastqs_fastQC.sh](https://github.com/luciamayorf/Data_preprocessing_alignment_v2/blob/main/scripts/fastqs_fastQC.sh) <input_fastq> to obtain the fastQC of each pair of reads.

Note: the sample list here have a very specific format, with the following columns: path/to/fastq(witought fastq_suffix not distinguising between read1 and read2), sample name, sex.

```bash
for input_fastq in $(cut -f1 /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/FASTQ_files/fastq_paths_samples_list_old_sequences.txt); do
  job_id=$(sbatch -c 6 --mem=4GB -t 03:00:00 /home/csic/eye/lmf/scripts/Data_preprocessing_alignment/fastqs_fastQC.sh ${input_fastq} | awk '{print $4}')
  echo "${job_id} ${input_fastq}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/fastqc/job_ids_fastqc_old_sequences.txt
done
```

Then we summarize the results with multiqc, running the script [multiqc_script.sh](https://github.com/luciamayorf/Data_preprocessing_alignment_v2/blob/main/scripts/multiqc_script.sh) <fastq_output_directory>

```bash
sbatch -t 00:15:00 -c 10 --mem 5GB /home/csic/eye/lmf/scripts/Data_preprocessing_alignment/multiqc_script.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/FASTQ_files/genome_proyect_backup/fastq_genome_project/fastqc
```
