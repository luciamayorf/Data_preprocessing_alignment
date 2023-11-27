# Data_preprocessing_alignment_variantcalling

This repository is an inspiration from Enrico's [Lynx_demography](https://github.com/Enricobazzi/Lynx_Demography) repository.

In this repository I keep the scripts used for sequencing data quality control (fastQC and multiQC), trimming (fastp), alignment (Enrico's [congenomics_fastq_align](https://github.com/Enricobazzi/congenomics_fastq_align) python package) and variant calling (DeepVariant).
 

## 1. Raw data quality control

First, we do a first FastQC analysis of the FASTQ files received with FastQC.

```
python scripts/make_fastqc_scripts.py config/all_rawreads_fastqs.yml
```

Once obtained the scripts, we submit them to CESGA ft3 (change the output directory accordingly to the run):
```
for sh in $(ls scripts/rawreads_fastqc/*.sh)
 do
  echo "sbatch $sh"
  sbatch -t 03:00:00 -c 20 --mem 5GB $sh
done
```

Finally, we summarize the obtained results running the MultiQC script.
```
sbatch -t 01:00:00 -c 20 --mem 5GB multiqc_script.sh <path/to/fastqc/output>
```


## 2. Reads trimming and quality control

We do a first trimming with fastp, followed by a quality control with FastQC and MultiQC of the trimmed reads.

### 2.1. Fastp

We run fastp
```
python scripts/make_fastp_scripts.py config/all_rawreads_fastqs.yml
```
  here I still need to change the path to the adapter fasta

```
for sh in $(ls scripts/fastp/rawreads_fastqc/*.sh)
 do
  echo "sbatch $sh"
  sbatch $sh
done
```

### 2.2. Trimming quality control

Then we assess the quality of the trimming:

```
python scripts/make_fastqc_scripts.py config/all_fastp_fastqs.yml
```

```
for sh in $(ls scripts/fastp/fastqc/*.sh)
 do
  echo "sbatch $sh"
  sbatch $sh
done
```

Finally, we summarize the obtained results running the MultiQC script.
```
sbatch -t 01:00:00 -c 20 --mem 5GB multiqc_script.sh <path/to/fastp_fastqc/output>
```

ME QUEDO AQUÍ

## Alignment and quality control

Finally, we align the reads to the reference genome using [congenomics_fastq_align](https://github.com/Enricobazzi/congenomics_fastq_align) to generate the scripts. After that, we perform a quality control of the aligment with qualimap, samtools and MultiQC.


```

```



## Variant calling

This part is only applicable to the high coverage dataset (50 individuals at ~30X). We use the WGS model from [DeepVariant](https://github.com/google/deepvariant) to call variants in our dataset.

We first use a bash script to create the samples bams dictionary:

```
bash ./scripts/make_deepvariant_dictionary.sh
```

Then we generate a deepvariant script per sample:
```
python scripts/make_deepvariant_scripts.py config/all_bams_calling.yml
```




