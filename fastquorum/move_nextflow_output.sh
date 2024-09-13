#!/bin/bash

module load samtools

# List of project IDs
# projects=("19610R" "20070R" "20963R")
projects=("19959R")

for project in "${projects[@]}"; do

    echo "Copying fastquorum output for $project..."

    if [ "$project" == "19959R" ]; then
        source_dir=/scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/nextflow-$project/output-$project-all
    else
        source_dir=/scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/nextflow-$project/output-$project
    fi    

    ## Define destination dir for copied bam files
    dest_dir=/scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/output

    ## Define the consensus bam dir
    bam_dir=$source_dir/filtering/align_consensus_bam

    for bam_file in "$bam_dir"/*/*.mapped.bam; do

        sample=$(basename $bam_file | cut -d'.' -f1)

        # Submit the job with sbatch
        sbatch --account=quinlan-rw \
            --partition=quinlan-shared-rw \
            --output=/scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/logs/copy_${sample}.out \
            --error=/scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/logs/copy_${sample}.err \
            --wrap="\
        # Copy the bam file
        echo \"Copying BAM file for $sample...\"
        cp \"$bam_file\" \"$dest_dir/bams\"

        # Index the bam file
        echo \"Indexing $bam_file\"
        samtools index \"$bam_file\"

        # Copy .bai file
        echo \"Copying $bam_file.bai to $dest_dir/bams\"
        cp \"$bam_file.bai\" \"$dest_dir/bams\""

        # ## Copy the bam filecd
        # echo "Copying BAM file for $sample..."
        # cp $bam_file $dest_dir/bams

        # ## Index the bam file
        # echo "Indexing $bam_file"
        # samtools index $bam_file

        # ## Copy .bai file
        # echo "Copying  $bam_file.bai to $dest_dir/bams"
        # cp $bam_file.bai $dest_dir/bams

    done

    ## Copy metrics files
    echo "Copying duplex metrics files from $source_dir to $dest_dir"
    cp -r $source_dir/metrics/duplex_seq/*/* $dest_dir/metrics

    ## Copy family size files
    echo "Copying family size files from $source_dir to $dest_dir"
    cp -r $source_dir/grouping/groupreadsbyumi/*/* $dest_dir/grouping

done