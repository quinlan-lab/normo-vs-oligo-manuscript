import os
import subprocess
import pandas as pd

## Directory of complete bam files
bam_dir = "/scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/output/bams/"
output_bam_dir = "/scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/output/bams/downsampled_bams"

## Read in manifest
manifest_file = "/scratch/ucgd/lustre-labs/quinlan/u1240855/normo-vs-oligo/data/sample_metadata.txt"
manifest_df = pd.read_csv(manifest_file, sep="\t", header=0)

## Get sperm and blood samples
blood_samples = manifest_df.loc[(manifest_df.genomex_id == "19610R") & (manifest_df.tissue == "Blood_DNA"), "sample"].tolist() 
sperm_samples = manifest_df.loc[(manifest_df.genomex_id == "19610R") & (manifest_df.tissue == "Sperm_DNA"), "sample"].tolist() 

## Function to create downsampled BAM file
def downsample(sample, fraction, bam_dir=bam_dir, output_bam_dir=output_bam_dir):
    input_file = os.path.join(bam_dir, f"{sample}.mapped.bam")
    output_file = os.path.join(output_bam_dir, f"{sample}_downsampled.bam")
    log_dir = "/scratch/ucgd/lustre-labs/quinlan/u1240855/normo-vs-oligo/logs"
    slurm_script = f"""#!/bin/bash
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --job-name=downsample_{sample}
#SBATCH --output={log_dir}/downsample_{sample}.out
#SBATCH --error={log_dir}/downsample_{sample}.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4

module load samtools
echo "Downsampling {input_file} by {fraction}"
samtools view -bs {fraction} {input_file} | samtools sort -@ 4 -o {output_file}
samtools index {output_file}
"""

    ## Write downsampling shell script
    script_dir = "/scratch/ucgd/lustre-labs/quinlan/u1240855/normo-vs-oligo/scripts/downsample_bams_shell_scripts/"
    with open(os.path.join(script_dir, f"downsample_{sample}.sh"), "w") as file:
        file.write(slurm_script)
        
    ## Submit downsampling job
    print(os.path.join(script_dir, f"downsample_{sample}.sh"))
    subprocess.run(["sbatch", os.path.join(script_dir, f"downsample_{sample}.sh")])

# Iterate over blood samples and downsample each one
for sample in blood_samples:
    downsample(sample, 0.5)

# Iterate over sperm samples and downsample each one
for sample in sperm_samples:
    downsample(sample, 0.33)
