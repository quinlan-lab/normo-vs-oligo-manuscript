#!/bin/bash

## Define arguments based on flags
while getopts "f:s:b:i:u:p:r:t:o:d:" opt; do
    case $opt in
        f) fasta=$OPTARG ;;
        s) sample=$OPTARG ;;
        b) bam_file=$OPTARG ;;
        i) probe_bed_file=$OPTARG ;;
        u) subset_probe_bed_file=$OPTARG ;;
        p) probe=$OPTARG ;;
        r) filter=$OPTARG ;;
        t) threshold=$OPTARG ;;
        o) output_file=$OPTARG ;;
        d) output_dir=$OPTARG ;;
        \?) echo "Invalid option: -$OPTARG" >&2 ;;
    esac
done

# --------------------------------------------------------------------------------

## Run pbr
echo "Running per-base read pileup (pbr; https://github.com/brentp/pbr)"
echo "Fasta: $fasta"
echo "Sample: $sample"
echo "BAM file: $bam_file"
echo "Probe bed file: $probe_bed_file"
echo "Subset probe bed file: $subset_probe_bed_file"
echo "Probe: $probe"
echo "Filter: $filter"
echo "Threshold: $threshold"
echo "Output file: $output_file"

pbr="/scratch/ucgd/lustre-work/quinlan/u1240855/DuplexSeqAnalysis/Snakemake/scripts/pbr"

# --------------------------------------------------------------------------------

## Determine which probe bed file to use
intervals=$probe_bed_file
if [[ "$probe" == "subset_intervals" ]]; then
    echo "Using subset probe interval set..."
    intervals=${subset_probe_bed_file}
fi

## Make temp file
tempfile=$(mktemp)

# ---------- APPLY DIFFERENT FILTER COMBINATIONS ----------

## Apply all filters
if [[ "$filter" == "all_filters" ]]; then

    if [[ "$threshold" == "baseline" ]]; then
        $pbr --threads 8 --fasta $fasta $bam_file \
        "return string_count(read.sequence, 'N') < 0.05 * read.length and read.distance_from_5prime >= 15 and read.distance_from_3prime >= 15 and bit32.band(read.flags, bit32.bor(4, 256, 512, 1024)) == 0" \
        --max-depth 100000 \
        -p "return pile.n / pile.depth < 0.05" \
        --bedfile $intervals > $tempfile
    fi

    if [[ "$threshold" == "strict" ]]; then
        $pbr --threads 8 --fasta $fasta $bam_file \
        "return string_count(read.sequence, 'N') < 0.025 * read.length and read.distance_from_5prime >= 15 and read.distance_from_3prime >= 15 and bit32.band(read.flags, bit32.bor(4, 256, 512, 1024)) == 0" \
        --max-depth 100000 \
        -p "return pile.n / pile.depth < 0.025" \
        --bedfile $intervals > $tempfile
    fi

    if [[ "$threshold" == "lenient" ]]; then
        $pbr --threads 8 --fasta $fasta $bam_file \
        "return string_count(read.sequence, 'N') < 0.10 * read.length and read.distance_from_5prime >= 10 and read.distance_from_3prime >= 10 and bit32.band(read.flags, bit32.bor(4, 256, 512, 1024)) == 0" \
        --max-depth 100000 \
        -p "return pile.n / pile.depth < 0.10" \
        --bedfile $intervals > $tempfile
    fi
fi

# --------------------------------------------------------------------------------

## Apply read_N_prop filter
if [[ "$filter" == "read_N_prop" ]]; then

    if [[ "$threshold" == "baseline" ]]; then
        $pbr --threads 8 --fasta $fasta $bam_file \
            "return string_count(read.sequence, 'N') < 0.05 * read.length" \
            --max-depth 100000 \
            --bedfile $intervals > $tempfile
    fi

    if [[ "$threshold" == "strict" ]]; then
        $pbr --threads 8 --fasta $fasta $bam_file \
            "return string_count(read.sequence, 'N') < 0.025 * read.length" \
            --max-depth 100000 \
            --bedfile $intervals > $tempfile
    fi 

    if [[ "$threshold" == "lenient" ]]; then
        $pbr --threads 8 --fasta $fasta $bam_file \
            "return string_count(read.sequence, 'N') < 0.10 * read.length" \
            --max-depth 100000 \
            --bedfile $intervals > $tempfile
    fi     
fi

## Apply site N frequency filter
if [[ "$filter" == "site_N_frequency" ]]; then

    if [[ "$threshold" == "baseline" ]]; then
        $pbr --threads 8 --fasta $fasta $bam_file \
        "return true" \
        -p "return pile.n / pile.depth < 0.05" \
        --max-depth 100000 \
        --bedfile $intervals > $tempfile
    fi

    if [[ "$threshold" == "strict" ]]; then
        $pbr --threads 8 --fasta $fasta $bam_file \
        "return true" \
        -p "return pile.n / pile.depth < 0.025" \
        --max-depth 100000 \
        --bedfile $intervals > $tempfile
    fi  

    if [[ "$threshold" == "lenient" ]]; then
        $pbr --threads 8 --fasta $fasta $bam_file \
        "return true" \
        -p "return pile.n / pile.depth < 0.10" \
        --max-depth 100000 \
        --bedfile $intervals > $tempfile
    fi          
fi 

## Filter based on read position
if [[ "$filter" == "read_position" ]]; then

    if [[ "$threshold" == "baseline" ]]; then
        $pbr --threads 8 --fasta $fasta $bam_file \
        "return read.distance_from_5prime >= 15 and read.distance_from_3prime >= 15" \
        --max-depth 100000 \
        --bedfile $intervals > $tempfile
    fi

    if [[ "$threshold" == "strict" ]]; then
        $pbr --threads 8 --fasta $fasta $bam_file \
        "return read.distance_from_5prime >= 30 and read.distance_from_3prime >= 30" \
        --max-depth 100000 \
        --bedfile $intervals > $tempfile
    fi

    if [[ "$threshold" == "lenient" ]]; then
        $pbr --threads 8 --fasta $fasta $bam_file \
        "return read.distance_from_5prime >= 5 and read.distance_from_3prime >= 5" \
        --max-depth 100000 \
        --bedfile $intervals > $tempfile
    fi    
fi 

# --------------------------------------------------------------------------------

## Get pbr depth when no filter filters are applied
if [[ "$filter" == "none_filters" ]]; then
    output_file=$output_dir/coverage/pbr_files/$sample.$probe.$filter.$threshold.nucleotide_count.annotated.bed.gz
    $pbr --threads 8 --fasta $fasta $bam_file \
    "return true" \
    --max-depth 100000 \
    --bedfile $intervals > $tempfile
fi

# --------------------------------------------------------------------------------

echo "Creating bed file and appending sample id"
cat $tempfile | tail -n +2 | awk -v OFS="\t" '{print $1,$2,$2+1,$3,$4,$5,$6,$7,$8,$9}' | awk -v sample="$sample" -v filter="$probe_$filter_$threshold" '{print $0, "\t", sample, "\t", filter}' | bgzip > $output_file
tabix -p bed ${output_file}
