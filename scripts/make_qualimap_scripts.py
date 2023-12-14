import argparse
import yaml
import os

'''
This script generates bash scripts for running Qualimap bamqc on multiple bam files.

Usage:
python make_fastqc_scripts.py <path/to/config_file.yaml>

The config file should be a YAML file with the following structure:
 - The 'sample_dict' key should be a dictionary of dictionaries.
   - The first level of keys should be the sample names.
     - The values for the second level of keys should be a list of three strings:
       the path to the folder containing the BAM files, the name of the BAM file and the sex of the individual (the sex won't be used here, but it be useful to reuse the dictionary for variant calling).

Here's an example of what the YAML file should look like:
Here's an example of what the YAML file should look like:
sample_dict:
  sample1:
    ['path/to/bam1', 'file1.bam', 'male']
  sample2:
    ['path/to/bam2', 'file2.bam', 'female']

The script will generate a bash script for each bam file the sample_dict.
Here's an example of the bash script that would be generated for a bam file:

#!/bin/bash
#SBATCH --job-name=alignment1_qualimap
#SBATCH --output=logs/qualimap/alignment1_qualimap.out
#SBATCH --error=logs/qualimap/alignment1_qualimap.err
#SBATCH --time=6:00:00
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=10

module load qualimap

qualimap bamqc sample1.bam --java-mem-size=10G -outfile alignment1_qualimap.html -outformat html -outdir /path/to/bams/qualimap_output/alignment1_qualimap
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
    bam_name_root = os.path.splitext(bam_name)[0]

    # Generate the qualimap command
    command = f'qualimap bamqc -bam {bam} --java-mem-size=10G -outfile {bam_name_root}_qualimap.html -outformat html -outdir {containing_folder}/qualimap_output/{bam_name_root}_qualimap'
        
    # Generate the bash script
    script = f'''\
#!/bin/bash
#SBATCH --job-name={bam_name_root}_qualimap
#SBATCH --output=logs/qualimap/{bam_name_root}_qualimap.out
#SBATCH --error=logs/qualimap/{bam_name_root}_qualimap.err
#SBATCH --time=6:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=10

module load qualimap

{command}
'''
        
    # Write the bash script to a file
    script_file = f'scripts/qualimap/{bam_name_root}_qualimap.sh'
    with open(script_file, 'w') as f:
        f.write(script)

