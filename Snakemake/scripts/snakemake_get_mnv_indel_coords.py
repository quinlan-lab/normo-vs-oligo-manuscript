## Load libraries

import pandas as pd
from cyvcf2 import VCF

# ---------- Setup ----------

sys.stdout = sys.stderr = open(snakemake.log[0], "wt")

sample: str = snakemake.wildcards.sample
vcf_file: str = snakemake.input.vcf_file[0]
output_file: str = snakemake.output.mnv_indel_coords

print(f"Sample: {sample}")
print(f"VCF file: {vcf_file}")
print(f"Output file: {output_file}")

vcf = VCF(vcf_file)

# ---------- GET COORDINATES FOR MNVs AND INDELS IN THE SAMPLE VCF ----------

"""
MNVs are alterations where the len(ref) == len(alt) and len(ref) > 1 and len(alt)> 1
TwinStrand recommendation: remove snvs that overlap a clonal MNV
Information needed: 
    - var_coord
    - clonality
    
INDELS are alterations where the len(ref) != len(alt) and len(ref) > 1 and len(alt)> 1
TwinStrand recommendation: remove snvs that are +/- 10bp away from a clonal indel
Information needed:
    - var_coord
    - clonality
"""


print("Getting coordinates for clonal MNVs and indels in the sample VCF...")

## Initialize list to store variant dictionaries
var_coord_list = []

## Itearate through each variant in the VCF
for var in vcf: 
    
    ## Get variant information
    chrom, start, stop = var.CHROM, var.start, var.end
    ref = var.REF
    alt = var.ALT[0]
    
    ## Determine variation type
    if (len(ref) == len(alt)) and (len(ref) > 1) and (len(alt) > 1): 
        variation_type = "mnv"
    elif (len(ref) != len(alt)) and ((len(ref) > 1) or (len(alt) > 1)):
        variation_type = "indel"
    elif (len(ref) == 1) and (len(alt) == 1): 
        variation_type = "snv"
        print(f"Variant: {chrom}:{start}-{stop}_{ref}>{alt}")
    else: 
        variation_type = "other"
        error_message = f"Variant {chrom}:{start}-{stop}_{ref}>{alt} is not a snv, mnv, or indel... Check this..."
        raise ValueError(error_message)
        
    ## Get the clonality of the variant
    var_depth_vardict =  int(var.INFO.get('VD'))
    
    ## Store in a dictionary
    var_dict = {'chr': chrom, 'start': start, 'stop': stop, 'ref': ref, 'alt': alt, 'variation_type': variation_type, 'vardict_depth': var_depth_vardict}
    var_coord_list.append(var_dict)
    
## Return list of variant dictionaries and write to output file
var_coord_df = pd.DataFrame(var_coord_list)
var_coord_df.to_csv(output_file, index = False, header=True, sep="\t")
    