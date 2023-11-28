import argparse
import yaml
import os

'''
This script generates bash scripts for running DeepVariant on multiple BAMs.

Usage:
python make_deepvariant_scripts.py <path/to/config_file.yaml>

The config file should be a YAML file with the following structure:
 - The 'sample_dict' key should be a dictionary of dictionaries.
   - The first level of keys should be the sample names.
     - The values for the second level of keys should be a list of three strings:
       the path to the folder containing the BAM files, the name of the BAM file and the sex of the individual.
       
Here's an example of what the YAML file should look like:
sample_dict:
  sample1:
    ['path/to/bam1', 'file1.bam', 'male']
  sample2:
    ['path/to/bam2', 'file2.bam', 'female']

The script will generate a bash script for each bam in the sample_dict.
Here's an example of the bash script that would be generated for a bam file:

#!/bin/bash
#SBATCH --job-name=alignment1tyu8jik_gvcf
#SBATCH --output=logs/deepvariant/alignment1_gvcf.out
#SBATCH --error=logs/deepvariant/alignment1_gvcf.err
#SBATCH --time=6:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:a100:1

module load cesga/2020 deepvariant

deepvariant_gpu --model_type=WGS --ref=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.scaffolds.revcomp.scaffolds.fa \
--reads=/path/to/bam --output_gvcf=/path/to/output/output.g.vcf.gz --num_shards=32
'''

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('config_file', help='Path to YAML config file')
args = parser.parse_args()

# Load the YAML file
with open(args.config_file, 'r') as f:
    data = yaml.safe_load(f)


# Loop through each sample bam in the YAML file
for sample in data['sample_dict']:
    # Get the containing folder the bam file name:
    containing_folder = data['sample_dict'][sample][0]
    bam = os.path.join(containing_folder, data['sample_dict'][sample][1])
    bam_name = data['sample_dict'][sample][1]
    bam_name_root = os.path.splitext(bam_name)[0]                                   ###### NO LO ESTÁ LEYENDO BIEN, NO ESTÁ QUITANDO LA EXTENSIÓN .bam, COGE TODA LA RUTA. REVISAR TB EL OTRO SCRIPT.
        
    # Generate the deepvariant command
    command = ' '.join([
      'deepvariant_gpu --model_type=WGS',
      '--ref=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/',
      f'--reads={bam}',
      f'--output_gvcf=mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/{sample}_mLynPar1.2_ref.g.vcf.gz',
      '--num_shards=32'
    ])
       
    # Check the sex of the sample and add --haploid-contigs flag if the sample is male:
    sex_flag = ''

    if data['sample_dict'][sample][2] == "male":
      sex_flag = ''.join([
        '--haploid-contigs="mLynPar1.2_ChrX,mLynPar1.2_ChrY,mLynPar1.2_ChrY_unloc_1,mLynPar1.2_ChrY_unloc_2,mLynPar1.2_ChrY_unloc_3,',
        'mLynPar1.2_ChrY_unloc_4,mLynPar1.2_ChrY_unloc_5,mLynPar1.2_ChrY_unloc_6,mLynPar1.2_ChrY_unloc_7,mLynPar1.2_ChrY_unloc_8,',
        'mLynPar1.2_ChrY_unloc_9,mLynPar1.2_ChrY_unloc_10,mLynPar1.2_ChrY_unloc_11,mLynPar1.2_ChrY_unloc_12"',
        ' --par_regions_bed="/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.chrX_PAR_sexChr.bed"'
      ])

    # Generate the bash script
    script = f'''\
#!/bin/bash
#SBATCH --job-name={sample}_{bam_name_root}_gvcf
#SBATCH --output=logs/deepvariant/{sample}_{bam_name_root}_gvcf.out			##DUDA: ¿RUTA ABSOLUTA O RELATIVA?
#SBATCH --error=logs/deepvariant/{sample}_{bam_name_root}_gvcf.err
#SBATCH --time=6:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:a100:1

module load cesga/2020 deepvariant/1.6.0

{command} {sex_flag}
'''

        
    # Write the bash script to a file
    script_file = f'scripts/deepvariant/{sample}_gvcf.sh'
    with open(script_file, 'w') as f:
      f.write(script)
            
            
            
            #### TENGO QUE AÑADIR LA RUTA A LAS REGIONES PSEUDOATOSÓMICAS (crear el archivo mLynPar1.2_par.bed, ver las regiones que acumulan más cobertura en el cromosoma X).