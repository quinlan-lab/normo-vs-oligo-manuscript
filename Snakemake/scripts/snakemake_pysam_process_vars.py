## Load libraries

import sys
import pysam
import pandas as pd
from pyfaidx import Fasta
from cyvcf2 import VCF
from itertools import chain
from collections import OrderedDict

sys.path.append("/scratch/ucgd/lustre-labs/quinlan/u1240855/DuplexSeqAnalysis/Snakemake/scripts")
import snakemake_pysam_generics as pysam_generics

# ---------- Setup ----------

sys.stdout = sys.stderr = open(snakemake.log[0], "wt")

sample: str = snakemake.wildcards.sample
pae_variant_file: str = snakemake.input.pae_variant_file
pbr_file: str = snakemake.input.sample_unfiltered_pbr_file
bam_file: str = snakemake.input.bam_file
vcf_file: str = snakemake.input.vcf_file[0]
fasta_file: str = snakemake.input.fasta_file
mnv_indel_coords_file: str = snakemake.input.mnv_indel_coords_file[0]
probe_bed_file: str = snakemake.input.probe_bed_file
output_file: str = snakemake.output.output_file
max_pysam_depth: int = snakemake.params.max_pysam_depth
flank_length: int = snakemake.params.flank_length

# sample = "19610X23"
# pae_variant_file = "/scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/workflow_process_duplex_bam/data/pae_variants_hg38.txt"
# pbr_file = "/scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/workflow_process_duplex_bam/output/v17/coverage/pbr_files/19610X23.all_intervals.none_filters.baseline.nucleotide_count.annotated.bed.gz"
# bam_file = "/scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/output/bams/downsampled_bams/19610X23.downsampled.bam"
# vcf_file = "/scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/workflow_process_duplex_bam/output/v17/fulcrum/19610X23.vardict.annotated.in_probe_intervals.vcf.gz"
# fasta_file = "/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta"
# mnv_indel_coords_file = "/scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/workflow_process_duplex_bam/output/v17/tmp/19610X23.mnv_indel_coords.txt"
# probe_bed_file = "/scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/workflow_process_duplex_bam/data/combined_probe_intervals.bed"
# output_file = "/scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/workflow_process_duplex_bam/output/v17/tmp/19610X23.pysam_processed.txt"
# max_pysam_depth = 25
# flank_length = 3

print(f"Sample: {sample}")
print(f"PAE variants file: {pae_variant_file}")
print(f"PBR file: {pbr_file}")
print(f"BAM file: {bam_file}")
print(f"VCF file: {vcf_file}")
print(f"Reference genome: {fasta_file}")
print(f"MNV and indel coordinates file: {mnv_indel_coords_file}")
print(f"Probe bed file: {probe_bed_file}")
print(f"Output file: {output_file}")
print(f"Max pysam depth: {max_pysam_depth}")
print(f"Flank length: {flank_length}")

## Read in files
print("Reading in the PAE variant file")
pae_var_df = pd.read_csv(pae_variant_file, sep="\t", header=0)
print("Reading in the PBR file")
pbr_df = pd.read_csv(pbr_file, sep="\t", header=None, compression='gzip')
pbr_df.columns = ["chr", "start", "stop", "ref_allele", "pbr_depth",
                  "a_count", "c_count", "g_count", "t_count", "n_count", 
                  "sample", "filter_threshold"]
print("Reading in the probe bed file")
probe_bed_df = pd.read_csv(probe_bed_file, sep="\t", header=None)
probe_bed_df.columns = ["chr", "start", "stop", "gene"]
print("Reading in the BAM file")
bam = pysam.AlignmentFile(bam_file, "rb")
print("Reading in the reference genome")
ref_genome = Fasta(fasta_file)
print("Reading in the VCF file")
vcf = VCF(vcf_file)
print("Reading in the MNV and indel coordinates file")
mnv_indel_coords = pd.read_csv(mnv_indel_coords_file, sep="\t", header=0)

## Define PAE genes
pae_genes = ['BRAF', 'RET', 'FGFR3', 'FGFR2', 'HRAS', 'PTPN11', 'KRAS', 'MAP2K1', 'MAP2K2', 'SOS1', 'RAF1', 'CBL']

# ---------- Function to get read information for specific variant ----------

## Function to initialize variant dictionary
def init_var_dict(var, ref_genome, probe_bed_df, pae_var_df, pae_genes, flank_length, no_call_rate_cutoff, mnv_indel_coords): 
    
    var_dict = OrderedDict()

    ## Get variant coordinates    
    chrom = var.CHROM    
    start = int(var.start)
    stop = int(var.end)
    ref = var.REF
    alt = var.ALT[0]
    coord = chrom + ":" + str(start) + "-" + str(stop) + "-" + ref + ">" + str(alt)
    
    ## Get the mutation type
    variation_type = pysam_generics.check_variation_type(variant = var)
    if variation_type != "snv": 
        print(f"Variant {coord} is not an SNV... It is a {variation_type}... Skipping...")
        return None
    else: 
        print(f"Variant {coord} is an SNV... Proceeding...")
        
    ## Check if pbr depth was measured 
    if var.INFO.get('pbr_depth') is None:
        print(f"Variant {coord} has no pbr depth information... Skipping...")
        return None

    ## At this point, the mutation should be an SNV
    assert var.INFO.get('pbr_ref') == var.REF
    assert (len(ref) == 1) and (len(alt) == 1)
    assert ref != alt
    var_dict['chr'], var_dict['start'], var_dict['stop'], var_dict['ref'], var_dict['alt'], var_dict['var_coord'] = chrom, start, stop, ref, alt, coord

    ## Determine gene and panel for the variant
    gene, is_pae_gene = pysam_generics.get_gene_name(variant=var, probe_bed_file=probe_bed_df, pae_genes=pae_genes)
    if (gene == "neutral_region"):
        panel = "mutagenesis"
    else: 
        panel = "custom"
    var_dict['gene'] = gene
    var_dict['is_pae_gene'] = is_pae_gene
    var_dict['is_pae_var'] = pysam_generics.get_pae_info(var_chrom=chrom, var_start=start, var_stop=stop, gene=gene, pae_var_df=pae_var_df)
    var_dict['panel'] = panel

    ## Get variant information from vardict
    vardict_alt_count = int(var.INFO.get('VD'))
    vardict_total_depth = int(var.INFO.get('DP'))    
    vardict_af = float(var.INFO.get('AF'))
    vardict_hiaf = float(var.INFO.get('HIAF'))
    var_dict['vardict_alt_count'], var_dict['vardict_total_depth'], var_dict['vardict_af'], var_dict['vardict_hiaf'] = vardict_alt_count, vardict_total_depth, vardict_af, vardict_hiaf

    ## Mutation subtype and trinucleotide context
    context = pysam_generics.get_trinucleotide_context(variant=var, ref_genome=ref_genome)
    collapsed_context = pysam_generics.get_collapsed_trinucleotide_context(context=context)
    is_cpg_site = pysam_generics.get_is_cpg(collapsed_context=collapsed_context)
    mutation_subtype_uncollapsed, mutation_subtype_collapsed = pysam_generics.get_mutation_subtype(variant=var, is_cpg_site=is_cpg_site)  
    var_dict['subtype_uncollapsed'], var_dict['subtype_collapsed'] = mutation_subtype_uncollapsed, mutation_subtype_collapsed
    var_dict['context'], var_dict['collapsed_context'], var_dict['is_cpg_site'] = context, collapsed_context, is_cpg_site
    
    ## Get the clinical annotation
    var_dict['clin_sig'] = str(var.INFO.get('CLNSIG'))
    var_dict['clin_disease'] = str(var.INFO.get('CLNDN'))

    ## Get filter information from vardict
    var_dict['vardict_filter'] = str(var.FILTER)
    
    ## Get twinstrand filter information
    left_hp_filter, left_seq = pysam_generics.check_flanking_hp(variant=var, alt_allele=alt, ref_genome=ref_genome, flank_length=flank_length, direction="left") ## Homopolymer filter (left)
    right_hp_filter, right_seq = pysam_generics.check_flanking_hp(variant=var, alt_allele=alt, ref_genome=ref_genome, flank_length=flank_length, direction="right")   ## Homopolymer filter (right)
    var_dict['left_flanking_seq'] = left_seq
    var_dict['right_flanking_seq'] = right_seq
    overlap_mnv_indel_filter = pysam_generics.overlap_clonal_mnv_indel(variant=var, mnv_indel_coords=mnv_indel_coords)    ## Overlapping clonal MNV or indel filter
    no_call_rate_filter = pysam_generics.check_no_call_rate(variant=var, sample_mean_no_call_rate_cutoff=no_call_rate_cutoff)     ## No call rate filter

    # Initialize an empty list
    twinstrand_filter = []

    # Append corresponding strings for True values
    if left_hp_filter:
        twinstrand_filter.append("left_hp_filter")
    if right_hp_filter:
        twinstrand_filter.append("right_hp_filter")
    if overlap_mnv_indel_filter:
        twinstrand_filter.append("overlap_mnv_indel_filter")
    if no_call_rate_filter:
        twinstrand_filter.append("no_call_rate_filter")
    
    # Join the list elements with a ";" separator
    twinstrand_filter = ";".join(twinstrand_filter)
    var_dict['twinstrand_filter'] = twinstrand_filter

    ## Get the pbr nucleotide counts
    var_dict['a_count'] = int(var.INFO.get('a_count'))
    var_dict['c_count'] = int(var.INFO.get('c_count'))
    var_dict['g_count'] = int(var.INFO.get('g_count'))
    var_dict['t_count'] = int(var.INFO.get('t_count'))
    var_dict['n_count'] = int(var.INFO.get('n_count'))
    var_dict['pbr_depth'] = int(var.INFO.get('pbr_depth'))
    var_dict['site_N_freq'] = var_dict['n_count'] / var_dict['pbr_depth']
    
    ## Get homopolymer information
    var_dict['hp_chr'] = str(var.INFO.get('hp_chr'))
    var_dict['hp_start'] = int(var.INFO.get('hp_start'))
    var_dict['hp_stop'] = int(var.INFO.get('hp_stop'))
    var_dict['hp_nucleotide'] = str(var.INFO.get('hp_nucleotide'))
    var_dict['hp_length'] = int(var.INFO.get('hp_length'))
    var_dict['hp_distance'] = int(var.INFO.get('hp_distance'))
    
    ## Get fraguracy error information
    fraguracy_error_count = var.INFO.get('fraguracy_errors')
    if fraguracy_error_count is None: 
        var_dict['fraguracy_error_count'] = 0
    else:
        var_dict['fraguracy_error_count'] = int(fraguracy_error_count)
    
    ## Initialize empty dictionaries with read information
    empty_read_dict = {'read_name': "no read", 'read_pos': None, 'read_orientation': None, 'read_length': None, \
                        'region_N_count': None, 'region_len': None, 'region_N_prop': None,
                                'read_N_prop': None, 'read_CIGAR': None, 'read_N_ct': None, 'strand': None, \
                                    'read_is_del': None, 'read_is_refskip': None, 'pair_status': None, \
                                        'mapping_status': None, 'duplicate_status': None, \
                                            'primary_status': None, 'supplementary_status': None, \
                                                'aD': None, 'bD': None, 'cD': None, \
                                                    'aE': None, 'bE': None, 'cE': None}
    
    final_var_dict = OrderedDict({**var_dict, **empty_read_dict})
    
    return final_var_dict

## Use pysam to extract information from variant supporting reads 
def process_vcf_variant(var, bam, ref_genome, 
                        probe_bed_df, pae_var_df, pae_genes,
                        max_pysam_depth, flank_length, no_call_rate_cutoff, mnv_indel_coords): 
    
    var_list = []
    tmp_var_list = []
    read_count = 0
    var_found_in_read = False

    var_dict = init_var_dict(var=var, 
                             ref_genome=ref_genome, 
                             probe_bed_df=probe_bed_df, 
                             pae_var_df=pae_var_df, 
                             pae_genes=pae_genes,
                             flank_length=flank_length, 
                             no_call_rate_cutoff=no_call_rate_cutoff, 
                             mnv_indel_coords=mnv_indel_coords)
    
    ## Skips if the variant is not an SNV
    if var_dict == None: 
        return None

    ## Gather all reads that overlap the variant position
    for pileupcolumn in bam.pileup(var_dict['chr'], var_dict['start'], var_dict['stop'], 
                                   ignore_overlaps=False, max_depth=999999, stepper='nofilter', 
                                   min_base_quality = 0): 
        
        ## Only examine pileup of reads that overlap the variant position
        if pileupcolumn.pos != var_dict['start']:
            continue
        
        print(f"Pileup column at position {pileupcolumn.pos}...")

        ## Iterate over reads at the pileup site
        for pileupread in pileupcolumn.pileups: 
            
            ################################################################
            ##### Get information if the read contains the alt allele #####
            ################################################################
            
            ## Get the position of the alternate site on the read 
            read_pos = pileupread.query_position
            
            ## Check to verify if a read position was determined, None is reported if ths position occurs in a soft-clipped region? 
            if (read_pos == None): 
                # print("read position is None")
                continue
            
            ## Get the read allele at the variant position
            read_allele = pileupread.alignment.query_sequence[read_pos]  

            ## Only look at reads with the variant in question (ignores reference and N alleles)
            if read_allele != var_dict['alt']: 
                continue
            
            ############################################################################                
            ##### Perform filters to determine if variant is found in a valid read #####
            ############################################################################

            ## Skip read if there is a deletion at the position
            if (pileupread.is_del != True) and (pileupread.is_del == True): 
                print(read_pos + " position in the read is a deletion or skipped in the reference genome.")
                continue

            ###############################################################
            ##### At this point, the variant is found in a valid read #####
            ###############################################################
            
            pysam_var_dict = var_dict.copy()     
            # print(f"New read, here is the initial read name: {pysam_var_dict['read_name']}")
            assert read_allele == var_dict['alt']
            var_found_in_read = True

            ## Tally up number of reads supporting variant
            read_count += 1
                                            
            ## Get the read sequence
            read_seq = pileupread.alignment.query_sequence

            ## Get N count in flank regions of the variant in the read
            region_N_ct, region_len, region_N_prop = \
                pysam_generics.get_Ncount_flank(read_seq=read_seq, read_pos=read_pos)

            ## Count the number of N's in the read
            N_ct = pysam_generics.get_Ncount(read_seq = read_seq)
            
            ## Get read information
            read_name = pileupread.alignment.query_name ## Get the read name
            read_pair_orientation = pysam_generics.extract_readpair(read = pileupread) # Get read pair orientation
            read_length = pileupread.alignment.query_alignment_length ## Get the length of the read
            N_prop = N_ct / read_length
            read_CIGAR = pileupread.alignment.cigarstring ## Get the read CIGAR string
            strandedness = pysam_generics.extract_strandedness(read = pileupread) ## Get the strandedness 
            read_isdel = pileupread.is_del
            read_isrefskip = pileupread.is_refskip
            pair_status = pysam_generics.extract_pair_status(read = pileupread) ## Get pair status
            mapping_status = pysam_generics.extract_mapping_status(read = pileupread) ## Get mapping status
            duplicate_status = pysam_generics.extract_duplicate_status(read = pileupread) ## Get duplicate status
            primary_status = pysam_generics.extract_primary_status(read = pileupread) ## Get primary/secondary status
            supplementary_status = pysam_generics.extract_supplementary_status(read = pileupread) ## Get supplementary statys
        
            pysam_var_dict['read_name'] = read_name
            pysam_var_dict['read_pos'] = read_pos + 1 ## read_pos is 0-base, convert to 1-base
            pysam_var_dict['read_orientation'] = read_pair_orientation
            pysam_var_dict['read_length'] = read_length
            pysam_var_dict['region_N_count'] = region_N_ct
            pysam_var_dict['region_len'] = region_len
            pysam_var_dict['region_N_prop'] = region_N_prop
            pysam_var_dict['read_N_prop'] = N_prop
            pysam_var_dict['read_CIGAR'] = read_CIGAR
            pysam_var_dict['strand'] = strandedness
            pysam_var_dict['read_is_del'] = read_isdel
            pysam_var_dict['read_is_refskip'] = read_isrefskip
            pysam_var_dict['pair_status'] = pair_status
            pysam_var_dict['mapping_status'] = mapping_status
            pysam_var_dict['duplicate_status'] = duplicate_status
            pysam_var_dict['primary_status'] = primary_status
            pysam_var_dict['supplementary_status'] = supplementary_status
            
            ## Get read tags
            aD, bD, cD, aE, bE, cE = pysam_generics.extract_tags(read = pileupread)
            pysam_var_dict.update({'aD': aD, 'bD': bD, 'cD': cD, 'aE': aE, 'bE': bE, 'cE': cE})
            
            ## Append dictionary to the temporary var list
            tmp_var_list.append(pysam_var_dict)
        
            ## Only include mutations with pysam depth less than user defined threshold                
            if (read_count > max_pysam_depth): 
                print(f"Variant {var_dict['var_coord']} appears inherited with pysam pileup {len(tmp_var_list)}")
                break
        
        if var_found_in_read == False: 
            print("No reads with variant " + var_dict['var_coord'] + " reported... Returning dictionary with empty read information...")
            assert read_count == 0 and var_dict['read_name'] == "no read"
            var_dict['pysam_depth'] = 0
            var_dict['var_source'] = "error"
            var_list.append(var_dict)
            
        elif (var_found_in_read == True) and (read_count > max_pysam_depth): 
            print(f"Variant {var_dict['var_coord']} had pysam depth of {read_count} which exceed the user defined threshold of {max_pysam_depth}...")
            assert len(tmp_var_list) == read_count
            var_dict['read_name'] = "read count exceeded max depth"
            var_dict['pysam_depth'] = read_count
            var_dict['var_source'] = "inherited"
            var_list.append(var_dict)
            
        else: 
            print("Variant " + var_dict['var_coord'] + " found in the BAM file...")
            assert (len(tmp_var_list) == read_count) and (read_count >= 1) and (read_count <= max_pysam_depth)
            tmp_var_list = [{**d, 'pysam_depth': read_count, 'var_source': "likely_denovo"} for d in tmp_var_list]
            var_list = tmp_var_list
                    
    return var_list
    
# ---------- Calculate ALT allele count from pysam read pileups ----------

def get_pysam_ac(var_list): 
    
    new_var_source = None
    
    ## If no pysam pileup 
    if var_list[0]['read_name'] == "no read": 
        assert (len(var_list) == 1) and (var_list[0]['var_source'] == "error")
        print(f"No pysam pileup detected for {var_list[0]['var_coord']}... Reporting an alt_count_v2 value of 0...")
        read_num = var_list[0]['pysam_depth']
        
    ## If mutation was inherited
    elif var_list[0]['read_name'] == "read count exceeded max depth": 
        assert (len(var_list) == 1) and (var_list[0]['var_source'] == "inherited")
        read_num = var_list[0]['pysam_depth']
        print(f"The variant {var_list[0]['var_coord']} is likely inherited... Reporting an alt_count_v2 value of {str(read_num)} which is the pysam-pileup allele count...")

    ## If likely de novo mutation was called
    else: 
        
        # [print(v['var_source']) for v in var_list]
        # print(all(v['var_source'] == "likely_denovo" for v in var_list))
        reads = [read['read_name'] for read in var_list if (\
            (read['supplementary_status'] == "not supplementary") & 
            (read['pair_status'] == "proper pair") & 
            (read['duplicate_status'] == "not duplicate") & 
            (read['mapping_status'] == "mapped") & 
            (read['primary_status'] == "primary"))]
        
        assert all(x not in reads for x in ['no read', 'read count exceeded max depth'])
        assert all(v['var_source'] == "likely_denovo" for v in var_list)
        
        
        ## TODO: Fix this. This is always returning a value of 1. 
        read_num = len(set(reads))
        
        ## If the reads failed --> change var_source
        if (read_num == 0): 
            print(f"The variant {var_list[0]['var_coord']} is was marked as likely denovo + had pysam pileup but supporting reads failed qc... Changing the var_source to error")
            new_var_source = "error"
    
    return read_num, new_var_source
    
# ---------- Calculate the sample's mean no call rate from the pbr output ----------
"""
Iterate through each site in the pbr_file and calculate the no call rate (n_count / depth)
"""

def get_sample_no_call_rate_cutoff(pbr_df):
    
        print("Calculating the sample's mean no call rate from the pbr output...")
        
        ## Get the sample's no call rate
        pbr_df['site_no_call_rate'] = pbr_df['n_count'] / pbr_df['pbr_depth']
        
        ## Calculate mean no call rate
        mean_sample_no_call_rate = pbr_df['site_no_call_rate'].mean()
        std_dev_sample_no_call_rate = pbr_df['site_no_call_rate'].std()
        
        ## Print the sample mean no call rate
        print(f"The sample's mean no call rate is {mean_sample_no_call_rate} with a standard deviation of {std_dev_sample_no_call_rate}...")
        
        return mean_sample_no_call_rate + std_dev_sample_no_call_rate

# ---------- GET VARIANT INFORMATION + FILTER ----------

if __name__ == "__main__":
        
        ## Get the mean no call rate for the sample
        no_call_rate_cutoff = get_sample_no_call_rate_cutoff(pbr_df=pbr_df)
        
        ## Get the pysam read information for each variant in the VCF
        sample_var_list = []

        ## Iterate through each variant in the VCF
        for var in vcf: 
            
            # if (var.CHROM != "chr4") or (var.start != 1806722) or (var.end != 1806723): 
            #     continue
            
            ## Get read information for the variant
            var_list = process_vcf_variant(var=var, 
                                           bam=bam, 
                                           ref_genome=ref_genome,
                                           probe_bed_df=probe_bed_df,
                                           pae_var_df=pae_var_df,
                                           pae_genes=pae_genes,
                                           max_pysam_depth=max_pysam_depth, 
                                           flank_length=flank_length,
                                           no_call_rate_cutoff=no_call_rate_cutoff, 
                                           mnv_indel_coords=mnv_indel_coords)
            
            if var_list == None: 
                continue
            
            ## Calculate pysam allele count
            alt_count_v2, new_var_source = get_pysam_ac(var_list=var_list)
            # print(var_list)
            # print(alt_count_v2)
            
            ## Update list of variant dictionaries with adjusted alt count value
            final_var_list = []
            for v in var_list:     
                if new_var_source is not None: 
                    v['var_source'] = new_var_source
                v['alt_count_v2'] = alt_count_v2
                v['AF_v2'] = v['alt_count_v2'] / v['pbr_depth']
                final_var_list.append(v)
            
            sample_var_list.append(final_var_list)
    
        ## Combine information from all variants 
        var_df = pd.DataFrame(list(chain.from_iterable(sample_var_list)))
        var_df.to_csv(output_file, sep="\t", index=False, header=True)

