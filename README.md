# Data_preprocessing_alignment
In this repository I keep the scripts used for sequencing data quality control (fastQC and multiQC), trimming (fastp) and alignment (Enrico's [congenomics_fastq_align](https://github.com/Enricobazzi/congenomics_fastq_align) python package).

It is mainly designed for the processing of low-coverage libraries generated with NovaSeqX. 

## Data quality control

First, we do a first FastQC analysis of the raw data received with FastQC, and summarize the results obtained with MultiQC.

```
python scripts/rawdataqc_scripts.py config/all_rawreads_fastqs.yml
```




## Reads trimming and quality control

We do a first trimming with fastp, followed by a quality control with FastQC and MultiQC of the trimmed reads.





## Alignment and quality control

Finally, we align the reads to the reference genome using [congenomics_fastq_align](https://github.com/Enricobazzi/congenomics_fastq_align) to generate the scripts. After that, we perform a quality control of the aligment with qualimap, samtools and MultiQC.
