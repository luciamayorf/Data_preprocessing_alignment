# Data_preprocessing_alignment_variantcalling

This repository is an inspiration from Enrico's [Lynx_demography](https://github.com/Enricobazzi/Lynx_Demography) repository.

In this repository, I keep the scripts used for sequencing data quality control ([fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [multiQC](https://multiqc.info/docs/)), trimming ([fastp](https://github.com/OpenGene/fastp)), alignment (Enrico's [congenomics_fastq_align](https://github.com/Enricobazzi/congenomics_fastq_align) python package), variant calling ([DeepVariant](https://github.com/google/deepvariant)) and gVCF merging ([GLnexus](https://github.com/dnanexus-rnd/GLnexus)).

---

## 1. Raw data quality control

First, we do a first FastQC analysis of the FASTQ files received with FastQC.

I run all the scripts from the package directoy in CESGA ft3:/home/csic/eye/lmf/Data_preprocessing_alignment.

```py
python scripts/make_fastqc_scripts.py config/all_rawreads_fastqs.yml
```

Once obtained the scripts, we submit them to CESGA ft3 (change the output directory accordingly to the run):

```bash
for sh in $(ls scripts/rawreads_fastqc/*.sh)
 do
  echo "sbatch $sh"
  sbatch -t 03:00:00 -c 20 --mem 5GB $sh
done
```

PRUEBO ESTA PARTE, FUNCIONA OK.

Finally, we summarize the obtained results running the MultiQC script.

```bash
sbatch -t 01:00:00 -c 20 --mem 5GB multiqc_script.sh <path/to/fastqc/output>
```

---

## 2. Reads trimming and quality control

We do a first trimming with fastp, followed by a quality control with FastQC and MultiQC of the trimmed reads.

### 2.1. Fastp

We run fastp:

```py
python scripts/make_fastp_scripts.py config/all_rawreads_fastqs.yml
```
  
  here I still need to change the path to the adapter fasta

```bash
for sh in $(ls scripts/fastp/rawreads_fastqc/*.sh)
 do
  echo "sbatch $sh"
  sbatch $sh
done
```

### 2.2. Trimming quality control

Then we assess the quality of the trimming:

```py
python scripts/make_fastqc_scripts.py config/all_fastp_fastqs.yml
```

```bash
for sh in $(ls scripts/fastp/fastqc/*.sh)
 do
  echo "sbatch $sh"
  sbatch $sh
done
```

Finally, we summarize the obtained results running the MultiQC script.

```bash
sbatch -t 01:00:00 -c 20 --mem 5GB multiqc_script.sh <path/to/fastp_fastqc/output>
```

---

## 3. Alignment and quality control

 After that, we perform a quality control of the aligment with qualimap, samtools and MultiQC.

### 3.1. Alignment

We align the reads to the reference genome using Enrico's [congenomics_fastq_align](https://github.com/Enricobazzi/congenomics_fastq_align) to generate the scripts. I downloaded the whole package and installed it in my $HOME in CESGA's ft3.

We start from a list of the samples, in this case, I obtain it from the first column of fastq sample list. Here, the YAML template not only contains a sample dictionary, but also the paths to the reference genome, the output folder, and the modules that need to be loaded to run the script.

```bash
for i in $(cut -f1 /path/to/fastq/sample/list.txt | cut -d'_' -f1,2 | sort -u); do
  python /home/csic/eye/lmf/alignments/congenomics_fastq_align-0.1.0/run_alignment.py --sample ${i} --config template_novogene_oct23.yml --test
done 
```

After generating those scripts, we need to add some lines to the script to load the modules and launch the jobs.

```bash
for script in /path/to/aligment/scripts/*_mLynPar1.2_ref_aligner.sh; do
    if [ -f "$script" ]; then
        sed -i '2i\#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/slurm-%j.err\n#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/slurm-%j.out\n' "$script"
        sed -i '6i\module load bwa\nmodule load samtools\nmodule load picard\nmodule load gatk/3.7-0-gcfedb67\n' "$script"
    fi
done
```

Then, we can launch the jobs:

```bash
for script in /path/to/aligment/scripts/*_mLynPar1.2_ref_aligner.sh; do
  echo "sbatch of ${script}"
  sbatch -t 05:00:00 -c 20 --mem 25GB ${script}    
done
```

### 3.2. Mapping quality control

We will use Qualimap bamqc. We first create a script for each bam file:

```py
python scripts/make_qualimap_scripts.py config/all_bams_qualimap_and_calling.yml
```

Then we launch the scripts:

```bash
for script in /path/to/qualimap/scripts/*_qualimap.sh; do
  echo "sbatch of ${script}"
  sbatch ${script}    
done
```

ME QUEDO AQUÍ. NO HE PROBADO EL FUNCIONAMIENTO DE NINGÚN SCRIPT EXCEPTO EL DE DEEPVARIANT.

---

## 4. Variant calling

This part is only applicable to the high coverage dataset (50 individuals at ~30X). We use the WGS model from [DeepVariant](https://github.com/google/deepvariant) to call variants in our dataset.

For the pseudoatosomal regions, we will establish a standard PAR1 region of 7 Mb in the beginning of the X and Y chromosomes, according to the existing bibliography: between [6 Mb](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5155386/)-[6.5Mb](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5522595/) and [10 Mb](https://onlinelibrary.wiley.com/doi/full/10.1111/eva.13302).

We first use a bash script to create the samples bams dictionary:

```bash
bash make_deepvariant_dictionary.sh
```

Then we generate a deepvariant script per sample:

```py
python scripts/make_deepvariant_scripts.py config/all_bams_qualimap_and_calling.yml
```

Finally, we run the scripts:

```bash
for script in /home/csic/eye/lmf/scripts/deepvariant/novogene_lp_sept23/*_gvcf.sh; do
  echo "sbatch of ${script}"
  sbatch ${script}    
done
```

---

## 5. gVCF merging

We will use GLnexus, which converts multiple VCF files to a single BCF file, which then needs to be converted to a VCF.

Initial test, should be something like this:

```bash
module load cesga/2020 glnexus/1.4.1
module load cesga/2020 bzip2/1.0.8
module load cesga/2020 samtools/1.14

glnexus_cli --config DeepVariantWGS /in/*.g.vcf.gz | bcftools view - | bgzip -@ 4 -c > output.vcf.gz
```

MIRAR LO DEL UNFILTERED MODO DE GLNEXUS-DV. MIRAR LOS INPUTS Y OUTPUTS DE TODOS LOS SCRITPS.
