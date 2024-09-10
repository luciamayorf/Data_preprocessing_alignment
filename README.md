# Data_preprocessing_alignment_v2

This repository is an inspiration from Enrico's [Lynx_demography](https://github.com/Enricobazzi/Lynx_Demography) repository.

In this repository, I keep the scripts used for sequencing data quality control ([fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [multiQC](https://multiqc.info/docs/)), trimming ([fastp](https://github.com/OpenGene/fastp)), alignment (Enrico's [congenomics_fastq_align](https://github.com/Enricobazzi/congenomics_fastq_align) python package) and quality control of the aligment ([qualimap](http://qualimap.conesalab.org/doc_html/analysis.html)).

** The scripts are only tested with one sample of the ones sequenced in previous projects, but they seem to work.

---

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

---

## 2. Reads trimming and quality control

We do a first trimming with fastp, followed by a quality control with FastQC and MultiQC of the trimmed reads.

We run fastp with the script [fastp_trimming.sh](https://github.com/luciamayorf/Data_preprocessing_alignment_v2/blob/main/scripts/fastp_trimming.sh) <input_fastq>:
```bash
for input_fastq in $(cut -f1 /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/FASTQ_files/fastq_paths_samples_list_old_sequences.txt); do
  job_id=$(sbatch -c 6 --mem=6GB -t 02:00:00 /home/csic/eye/lmf/scripts/Data_preprocessing_alignment/fastp_trimming.sh ${input_fastq} | awk '{print $4}')
  echo "${job_id} ${input_fastq}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/fastp/job_ids_fastp_old_sequences.txt
done
```

### Trimming quality control

Then we assess the quality of the trimming with the script [fastp_fastqc.sh](https://github.com/luciamayorf/Data_preprocessing_alignment_v2/blob/main/scripts/fastp_fastqc.sh) <input_fastq>:

```bash
for input_fastq in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/FASTQ_files/genome_proyect_backup/fastq_genome_project/fastp/*.fastp.fastq.gz); do
  job_id=$(sbatch -c 6 --mem=4GB -t 03:00:00 /home/csic/eye/lmf/scripts/Data_preprocessing_alignment/fastp_fastqc.sh ${input_fastq} | awk '{print $4}')
  echo "${job_id} ${input_fastq}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/fastqc/job_ids_fastqc_old_sequences.txt
done
```

Finally, we summarize the obtained results running [multiqc_script.sh](https://github.com/luciamayorf/Data_preprocessing_alignment_v2/blob/main/scripts/multiqc_script.sh) <fastq_output_directory>.

```bash
sbatch -t 00:15:00 -c 10 --mem 5GB /home/csic/eye/lmf/scripts/Data_preprocessing_alignment/multiqc_script.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/FASTQ_files/genome_proyect_backup/fastq_genome_project/fastp/fastqc
```

---

## 3. Alignment and quality control

After that, we perform a quality control of the aligment with qualimap, samtools and MultiQC.

We align the reads to the reference genome using Enrico's package [congenomics_fastq_align](https://github.com/Enricobazzi/congenomics_fastq_align) to generate the scripts. I downloaded the whole package and installed it in my $HOME in CESGA's ft3. I work from directory "/home/csic/eye/lmf/alignments/congenomics_fastq_align-0.1.0"

For that, we need to generate a sample dictionary using the custom script [make_fastqs_dictionary.sh](https://github.com/luciamayorf/Data_preprocessing_alignment_v2/blob/main/scripts/make_fastqs_dictionary.sh). This script requires a list where the first column contains the path and the pair of fastq prefix (no .fq.gz) and the second column contains the final sample name, tab-separated (Careful with the CODE_IDFQ variable definition, it highly depends on the fastq name format!!!).

> /path/to/FASTQ_files/LYNX_06_08/C5TMUACXX_2_1nf     c_lp_sm_0134
/path/to/FASTQ_files/LYNX_06_08/C5TN1ACXX_7_1nf     c_lp_sm_0134
/path/to/FASTQ_files/LYNX_06_08/C5TMUACXX_2_2nf     c_lp_do_0141


We start from a list of the samples, in this case, I obtain it from the second column of fastq sample list. Here, the YAML template not only contains a sample dictionary, but also the paths to the reference genome, the output folder, and the modules that need to be loaded to run the script.

```bash
for i in $(cut -f2 /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/FASTQ_files/fastq_paths_samples_list_old_sequences.txt | sort -u); do
  python ./run_alignment.py --sample ${i} --config template_alignment_old_sequences.yml --test
done

mv *_aligner.sh /home/csic/eye/lmf/scripts/alignment_sh/old_sequences/
```

After generating those scripts, we need to add some lines to the script to load the modules and launch the jobs.
```bash
for script in /home/csic/eye/lmf/scripts/alignment_sh/old_sequences/*_mLynPar1.2_ref_aligner.sh; do
    if [ -f "$script" ]; then
        sed -i '2i\#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/alignments/slurm-%j.err\n#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/alignments/slurm-%j.out\n' "$script"
        sed -i '6i\module load bwa\nmodule load samtools\nmodule load picard\nmodule load gatk/3.7-0-gcfedb67\n' "$script"
    fi
done
```

Then, we can launch the jobs:
```bash
for script in /home/csic/eye/lmf/scripts/alignment_sh/old_sequences/*_mLynPar1.2_ref_aligner.sh; do
  echo "sbatch of ${script}"
  sbatch -t 05:00:00 -c 20 --mem 25GB ${script}    
done
```

### Alignment quality control

We will use Qualimap bamqc with the script [bams_qualimap.sh](https://github.com/luciamayorf/Data_preprocessing_alignment_v2/blob/main/scripts/bams_qualimap.sh) <input_bam>. 

```{bash}
for input_bam in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/old_sequences/*_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.bam); do
  job_id=$(sbatch -c 10 --mem=20GB -t 06:00:00 /home/csic/eye/lmf/scripts/Data_preprocessing_alignment/bams_qualimap.sh ${input_bam} | awk '{print $4}')
  echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/qualimap/job_ids_qualimap_old_sequences.txt
done

sbatch -t 00:30:00 -c 10 --mem 5GB /home/csic/eye/lmf/scripts/Data_preprocessing_alignment/multiqc_script.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/old_sequences/qualimap_output
```

#### Effective coverage estimation

Using one of the outputs of Qualimap, we want to estimate the effective coverage, defined as the -ln (1 - f_cov) by [Steward et al. (2024)](https://www.biorxiv.org/content/10.1101/2024.01.30.578044v1), where f_cov is the fraction of the reference genome covered by at least one read. We will use the custom script [effective_cov_qualimap.sh](https://github.com/luciamayorf/Data_preprocessing_alignment_v2/blob/main/scripts/effective_cov_qualimap.sh) <input_qualimap_directory>

```{bash}
sbatch -t 00:10:00 --mem 1GB /home/csic/eye/lmf/scripts/alignmentQC_sh/effective_cov_qualimap.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/pool_epil_all/qualimap_outpu
```


---
