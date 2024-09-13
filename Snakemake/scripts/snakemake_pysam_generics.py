# ---------- FILTERS TO APPLY TO VARIANTS ----------

## EXCLUDE NON-SNVS
def check_variation_type(variant):
    ref_len = len(variant.REF)
    alt_len = len(variant.ALT[0])
    
    if ref_len == 1 and alt_len == 1: 
        variation_type = "snv"
    elif ref_len > alt_len: 
        variation_type = "del"
    elif ref_len < alt_len:
        variation_type = "ins"
    elif ref_len > 1 and alt_len > 1 and ref_len == alt_len: 
        variation_type = "mnv"
    else: 
        raise ValueError("Unknown variation type for: {variant.REF}>{variant.ALT[0]}")
    return variation_type

## EXCLUDE SNVS WHERE EITHER FLANK IS A HOMOPOLYMER + ALT IS THE HOMOPOLYMER BASE
def check_flanking_hp(variant, alt_allele, ref_genome, flank_length, direction):
    
    ## Subtract 1 because pyfaidx uses 0-based positions
    variant_pos = variant.end - 1
    
    ## Get the 3 nucleotides to the left of the position
    if direction == "left":
        hp = ref_genome[variant.CHROM][variant_pos - flank_length:variant_pos].seq
    elif direction == "right": 
        hp = ref_genome[variant.CHROM][variant_pos + 1:variant_pos + flank_length + 1].seq
    else: 
         raise ValueError("Flanking homopolymer direction must be 'left' or 'right'")
   
    ## Check if the sequence is a homopolymer
    is_hp = all(x == hp[0] for x in hp)

    ## Return true if flanking sequence is a homopolymer and the alt allele is the homopolymer base
    if is_hp and alt_allele == hp[0]:
        return True, hp
    else: 
        return False, hp

## EXCLUDE VARIANTS OVERLAPPING CLONAL MNVs OR WITHIN 10BP OF CLONAL INDELS
def overlap_clonal_mnv_indel(variant, mnv_indel_coords):
    clonal_mnvs = mnv_indel_coords.loc[(mnv_indel_coords.chr == variant.CHROM) &
                                       (mnv_indel_coords.start == variant.start) &
                                       (mnv_indel_coords.stop == variant.end) &
                                       (mnv_indel_coords.vardict_depth > 1) & 
                                       (mnv_indel_coords.variation_type == 'mnv')].shape[0]
    
    clonal_indels = mnv_indel_coords.loc[(mnv_indel_coords.chr == variant.CHROM) &
                                        (mnv_indel_coords.start - 10 <= variant.start) &
                                        (mnv_indel_coords.stop + 10 >= variant.end) &
                                        (mnv_indel_coords.vardict_depth > 1) &
                                        (mnv_indel_coords.variation_type == 'indel')].shape[0]
                                    
    if (clonal_mnvs > 0) or (clonal_indels > 0): 
        return True
    else:     
        return False

## EXCLUDE VARIANTS >1 SD AWAY FORM MEAN NO CALL RATE ACROSS ALL SEQUENCED SITES IN THE SAMPLE
def check_no_call_rate(variant, sample_mean_no_call_rate_cutoff):
    print(f"Variant: {variant.CHROM}:{variant.start}-{variant.end}_{variant.REF}>{variant.ALT[0]}")
    print(f"No calls: {variant.INFO.get('n_count')}")
    print(f"PBR depth: {variant.INFO.get('pbr_depth')}")
    var_n_freq = int(variant.INFO.get('n_count')) / int(variant.INFO.get('pbr_depth'))
    
    if var_n_freq > sample_mean_no_call_rate_cutoff: 
        return True
    else: 
        return False

# ---------- FUNCTIONS TO ANNOTATE VARIANT ENTRIES ----------

## GET MUTATION SUBTYPE FOR SNVS
def get_mutation_subtype(variant, is_cpg_site):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    ref, alt = variant.REF, variant.ALT[0]
    mutation_subtype_uncollapsed = f"{ref}>{alt}"
    
    ## Collapse the mutation type to C or T
    if ref in ['A', 'G']: 
        ref = complement[ref]
        alt = complement[alt]
    
    if ref == 'C' and alt == 'T' and is_cpg_site: 
        mutation_subtype_collapsed = "CpG>TpG"
    else: 
        mutation_subtype_collapsed = f"{ref}>{alt}"
    return mutation_subtype_uncollapsed, mutation_subtype_collapsed

## GET THE TRINUCLEOTIDE CONTEXT FOR THE SNV
def get_trinucleotide_context(variant, ref_genome):
    """
    Get the trinucleotide context for a given variant.
    """
    
    # Subtract 1 from the variant_pos because pysam uses 0-based positions
    variant_pos = variant.end - 1

    # Get the base before, at, and after the variant position
    context = ref_genome[variant.CHROM][variant_pos - 1:variant_pos + 2].seq
    
   
    return context.upper()

## GET THE COLLAPSED TRINUCLEOTIDE CONTEXT FOR SNVS 
def get_collapsed_trinucleotide_context(context):
    """
    Get the collapsed trinucleotide context. 
    If the central nucleotide is A or G, return the reverse complement of the trinucleotide context.
    For example, a CGA context would be converted to TCG --> and would be marked as a CpG site in the get_is_cpg function.
    """
    
    ## Define the complement of each nucleotide
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    
    ## Check if the central nucleotide is A or G
    if context[1] in ['A', 'G']: 
        collapsed_context = "".join(complement[base] for base in reversed(context))
    else: 
        collapsed_context = context
    
    return collapsed_context
        
        
## GET IS_CPG STATUS
def get_is_cpg(collapsed_context):
    """
    If the central nucleotide is a C and the last nucleotide is a G, return True.
    """
    if collapsed_context[1] == 'C' and collapsed_context[2] == 'G': 
        return True
    else: 
        return False

## Get the gene name for a variant
def get_gene_name(variant, probe_bed_file, pae_genes):
    
    ## Get the gene name for the variant based on the interval it falls under
    gene = probe_bed_file.loc[(probe_bed_file.chr == variant.CHROM) & 
                       (variant.start >= probe_bed_file.start) &
                       (variant.end <= probe_bed_file.stop), 'gene'].values[0]
    
    ## Check if the gene is a known PAE gene
    if gene in pae_genes: 
        is_pae_gene = True  
    else: 
        is_pae_gene = False
    
    return gene, is_pae_gene

## GET PAE INFORMATION FOR SNVS (GENE, PAE_EFFECT)
def get_pae_info(var_chrom, var_start, var_stop, gene, pae_var_df):

    ## Check if variant is in the pae_var_df
    
    if ((pae_var_df['chr'] == var_chrom) & 
        (pae_var_df['start'] == var_start) & 
        (pae_var_df['stop'] == var_stop) &
        (pae_var_df['gene'] == gene)).any():
        is_pae_var = True
    else: 
        is_pae_var = False
        
    return is_pae_var

# ---------- Functions to extract information from pysam reads ----------

## Function to extract tag information from a read
def extract_tags(read): 
    aD = read.alignment.get_tag("aD")
    bD = read.alignment.get_tag("bD")
    cD = read.alignment.get_tag("cD")

    aE = read.alignment.get_tag("aE")
    bE = read.alignment.get_tag("bE")
    cE = read.alignment.get_tag("cE")

    return aD, bD, cD, aE, bE, cE

## Function to extract strandedness from a read
def extract_strandedness(read): 
    if read.alignment.is_reverse == True: 
        strandedness = "reverse"

    else: 
        strandedness = "forward"

    return strandedness

## Function to determine if read is paired/in a proper pair
def extract_pair_status(read):
    if read.alignment.is_proper_pair == True: 
        pair_status = "proper pair"
    else: 
        pair_status = "not proper pair"

    return pair_status 

## Functon to determine if read is mapped
def extract_mapping_status(read): 
    if read.alignment.is_unmapped == True: 
        mapping_status = "not_mapped"
    else: 
        mapping_status = "mapped"
    
    return mapping_status

## Function to determine if read is a duplicate
def extract_duplicate_status(read): 
    if read.alignment.is_duplicate == True: 
        duplicate_status = "duplicate"
    else: 
        duplicate_status = "not duplicate"

    return duplicate_status


## Function to extract primary/secondary status
def extract_primary_status(read): 
    if read.alignment.is_secondary == True: 
        primary_status = "secondary"

    else: 
        primary_status = "primary"

    return primary_status

## Function to extract supplementary status
def extract_supplementary_status(read): 
    if read.alignment.is_supplementary == True: 
        supplementary_status = "supplementary"
    else: 
        supplementary_status = "not supplementary"

    return supplementary_status

## Function to extract read pair orientation
def extract_readpair(read): 
    if read.alignment.is_read1 == True: 
        read_pair_orientation = "read1"

    else: 
        read_pair_orientation = "read2"

    return read_pair_orientation

## Function to count # of N's in a read
def get_Ncount(read_seq): 
    N_ct = read_seq.count('N')

    return N_ct

## Function to get proportion of N's in a read 
def get_Ncount_flank(read_seq, read_pos, flank=10): 
    left_pos = read_pos - flank
    right_pos = read_pos + flank

    if (left_pos < 0): 
        left_pos = 0

    if (right_pos > len(read_seq)): 
        right_pos = len(read_seq)

    region_seq = read_seq[left_pos:right_pos+1]
    region_N_ct = get_Ncount(region_seq)
    region_len = len(region_seq)
    region_N_prop = region_N_ct / region_len

    return region_N_ct, region_len, region_N_prop