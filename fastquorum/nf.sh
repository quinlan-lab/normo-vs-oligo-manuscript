## Load libraries

module load jdk/11 ## load Java version 11
module load singularity

## Create conda environment and install nextflow (v23.10.1)
# conda create -n fastquorum_env
# conda activate fastquorum_env
# conda install -c bioconda nextflow

#################
##### SETUP #####
#################

set -euo pipefail

export TOWER_WORKSPACE_ID=52086791413503
export TOWER_ACCESS_TOKEN=eyJ0aWQiOiA5MTEwfS5kMmE0YzBhMzliNGMxYjUwMTJlZWE1ZDJjZTIwYjc2NGZkMDYwOTdh

ref38=/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta
ref38_fai=/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta.fai
experiment="19610R"
sample_sheet=/scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/sample_sheets/${experiment}_sample_sheet.csv

#############################################################
##### RUN NF-CORE/FASTQUORUM ON ALL SAMPLES IN PARALLEL #####
#############################################################

## aab72e6473 (https://github.com/nf-core/fastquorum/commit/aab72e6473ae3362a7ed00b51af09149d0b2015b)
nextflow run -r aab72e6473 nf-core/fastquorum \
    -with-tower \
    -profile redwood \
    -c ./nextflow-$experiment.config \
    -with-trace \
    -resume \
    --input $sample_sheet \
    --outdir output-$experiment \
    --mode 'ht' \
    --fasta $ref38 \
    --fasta_fai $ref38_fai \
    --duplex_seq \
    --groupreadsbyumi_strategy 'Paired' \
    --groupreadsbyumi_edits 1 \
    --call_min_reads '3 3 3' \
    --call_min_baseq 20 \
    --filter_min_reads '3 3 3' \
    --filter_min_baseq 45 \
    --filter_max_base_error_rate 0.2

cat <<EOF

The raw FASTQ pair contains the read structure 8M1S+T which breaks down into:

8M: 8bp of molecular index (UMI)
1S: 1bp of skip sequence
+T: A variable number of template sequence based on the Illumina kit that was used


Input/output options
  --input                      [string]  Path to comma-separated file containing information about the samples in the experiment.
  --outdir                     [string]  The output directory where the results will be saved. You have to use absolute paths to storage on Cloud
                                         infrastructure.
  --email                      [string]  Email address for completion summary.
  --multiqc_title              [string]  MultiQC report title. Printed as page header, used for filename if not otherwise specified.

Main options
  --mode                       [string]  The pipeline mode to use, either 'rd' for R&D or 'ht' for High-Throughput
  --duplex_seq                 [boolean] Enable when the input is duplex sequenecing.

Read grouping options
  --groupreadsbyumi_strategy   [string]  Grouping strategy
  --groupreadsbyumi_edits      [integer] Maximum number of edits

Consensus reads options
  --call_min_reads             [string]  Minimum reads to call a consensus
  --call_min_baseq             [integer] Minimum input base quality

Consensus filtering options
  --filter_min_reads           [string]  Minimum reads to keep a consensus
  --filter_min_baseq           [integer] Minimum consensus base quality
  --filter_max_base_error_rate [number]  The maximum error rate for a single consensus base

Reference genome options
  --genome                     [string]  Name of iGenomes reference.
  --fasta                      [string]  Path to FASTA genome file.
  --fasta_faix                 [string]  Path to FASTA reference index.
EOF
