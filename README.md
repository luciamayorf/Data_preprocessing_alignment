# Data_preprocessing_alignment

This repository is an inspiration from Enrico's [Lynx_demography](https://github.com/Enricobazzi/Lynx_Demography) repository.

In this repository I keep the scripts used for sequencing data quality control (fastQC and multiQC), trimming (fastp) and alignment (Enrico's [congenomics_fastq_align](https://github.com/Enricobazzi/congenomics_fastq_align) python package).

It is mainly designed for the processing of low-coverage libraries generated with NovaSeqX. 

## Raw data quality control

First, we do a first FastQC analysis of the FASTQ files received with FastQC.

```
python scripts/make_fastqc_scripts.py config/all_rawreads_fastqs.yml
```

Once obtained the scripts, we submit them to CESGA ft3:
```
for sh in $(ls scripts/rawreads_fastqc/*.sh)
 do
  echo "sbatch $sh"
  sbatch -t 03:00:00 -c 20 --mem 5GB $sh
done
```
OJO, TENGO QUE CAMBIAR DONDE SE GUARDAN LOS SCRIPTS EN FUNCIÓN DE CADA LIBRERÍA, PARA NO CORRERLOS TODOS (LOS YA HECHOS PREVIAMENTE) OTRA VEZ:


Finally, we summarize the obtained results running the MultiQC script.
```
sbatch -t 01:00:00 -c 20 --mem 5GB multiqc_script.sh <path/to/fastqc/output>
```





## Reads trimming and quality control

We do a first trimming with fastp, followed by a quality control with FastQC and MultiQC of the trimmed reads.





## Alignment and quality control

Finally, we align the reads to the reference genome using [congenomics_fastq_align](https://github.com/Enricobazzi/congenomics_fastq_align) to generate the scripts. After that, we perform a quality control of the aligment with qualimap, samtools and MultiQC.
