#!/bin/bash
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --job-name=downsample_19610X3
#SBATCH --output=/scratch/ucgd/lustre-labs/quinlan/u1240855/normo-vs-oligo/logs/downsample_19610X3.out
#SBATCH --error=/scratch/ucgd/lustre-labs/quinlan/u1240855/normo-vs-oligo/logs/downsample_19610X3.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4

module load samtools
echo "Downsampling /scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/output/bams/19610X3.mapped.bam by 0.33"
samtools view -bs 0.33 /scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/output/bams/19610X3.mapped.bam | samtools sort -@ 4 -o /scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/output/bams/downsampled_bams/19610X3_downsampled.bam
samtools index /scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/output/bams/downsampled_bams/19610X3_downsampled.bam
