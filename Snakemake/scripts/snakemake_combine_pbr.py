## Load libraries

import polars as pl

# ---------- Setup ----------

sys.stdout = sys.stderr = open(snakemake.log[0], "wt")

mutations_file = snakemake.input.mutations[0]

output_file = snakemake.output.combined_pbr

samples = snakemake.params.samples
filters = snakemake.params.filters
thresholds = snakemake.params.thresholds
probes = snakemake.params.probes
pbr_dir = snakemake.params.pbr_dir

print(f"PBR directory: {pbr_dir}")
print(f"Samples: {samples}")
print(f"Mutation file: {mutations_file}")
print(f"Filters: {filters}")
print(f"Thresholds: {thresholds}")
print(f"Probes: {probes}")
print(f"Output file: {output_file}")

## Get probe, filter, and threshold combinations
combinations = [f"{probe}.{applied_filter}.{threshold}" for probe in probes for applied_filter in filters for threshold in thresholds]
print(combinations)

## Read in mutations and get recurrent sites
mutations = pl.read_csv(mutations_file, separator="\t", has_header=True)
recurrent_sites = mutations.filter((pl.col("within_donor_recurrent") == True) | 
                                   (pl.col("within_cohort_recurrent") == True) |
                                   (pl.col("across_cohort_recurrent") == True)).select(["chr", "start", "stop"]).unique(subset=["chr", "start", "stop"])
recurrent_sites = recurrent_sites.with_columns(
    recurrent_status = pl.lit("recurrent")
)

print("Dataframe for recurrently mutated sites: ")
print(recurrent_sites)

## Function to get the file name 
def get_file_name(pbr_dir, sample, probe, applied_filter, threshold):
    filename = f"{pbr_dir}/{sample}.{probe}.{applied_filter}.{threshold}.nucleotide_count.annotated.bed.gz"
    print(f"Filename: {filename}")
    return filename

# ---------- Iterate through each file and calculate depth ----------

agg_pbr_list = []

## Iterate through each sample
for sample in samples: 
    
    print(f"Processing pbr output for sample: {sample}")
                
    ## Iterate through each filter combination
    for combination in combinations:
    
        ## Parse combination
        probe, applied_filter, threshold = combination.split(".")
        
        print(f"Processing {probe} {applied_filter} {threshold}")
        
        ## Get file name
        pbr_file = get_file_name(pbr_dir, sample, probe, applied_filter, threshold)
        
        ## Read in pbr file
        colnames = ["chr", "start", "stop", "ref_base", "depth", "a_count", "c_count", "g_count", "t_count", "n_count", "sample", "threshold"]
        pbr_data = pl.read_csv(pbr_file, separator="\t", has_header=False, new_columns=colnames)
        
        ## Sum the depth
        depth = pbr_data['depth'].sum()
        pbr_data = pbr_data.with_columns(
            total_depth = pl.lit(depth)
        )
        pbr_data = pbr_data.with_columns(
            probe_filter_threshold = pl.lit(combination)
        )
        
        ## Get relevant columns of interest
        sample_pbr_data = pbr_data.select(["sample", "total_depth", "threshold", "probe_filter_threshold"]).unique()            
            
        ## assert there is one row in pbr_data
        assert sample_pbr_data.shape[0] == 1, "DataFrame has more than one row"        
        
        ## Concatenate all DataFrames in the list
        agg_pbr_list.append(sample_pbr_data)

        ## Also get sample coverage after excluding recurrently mutated sites if all_filters is selected
        if (applied_filter == "all_filters"): 
            rm_rec_pbr_data = pbr_data.join(recurrent_sites, on=["chr", "start", "stop"], how="left", coalesce=True)
            rm_rec_pbr_data = rm_rec_pbr_data.with_columns(
                pl.col("recurrent_status").replace(None, "non-recurrent")
            )
            rm_rec_pbr_data = rm_rec_pbr_data.filter(pl.col("recurrent_status") == "non-recurrent")
            
            print("PBR dataframe with coverage removed at recurrently mutated sites: ")
            print(rm_rec_pbr_data)
            
            ## Sum the depth
            depth = rm_rec_pbr_data['depth'].sum()
            rm_rec_pbr_data = rm_rec_pbr_data.with_columns(
                total_depth = pl.lit(depth)
            )
            rm_rec_pbr_data = rm_rec_pbr_data.with_columns(
                probe_filter_threshold = pl.lit(combination + "_remove_recurrence")
            )

            ## Get relevant columns of interest
            exclude_rec_pbr = rm_rec_pbr_data.select(["sample", "total_depth", "threshold", "probe_filter_threshold"]).unique()
            print("Dataframe with coverage removed at recurrently mutated sites:")
            print(exclude_rec_pbr)
            
            ## assert there is one row in pbr_data
            assert exclude_rec_pbr.shape[0] > 0, "DataFrame is empty"
            assert exclude_rec_pbr.shape[0] == 1, "DataFrame has more than one row"
            
            ## Concatenate all DataFrames in the list
            agg_pbr_list.append(exclude_rec_pbr)
        
## Merge all pandas dataframes in the list
pbr_df = pl.concat(agg_pbr_list, how="vertical_relaxed")

print("Final dataframe: ")
print(pbr_df)
    
## Convert to dataframe
# pbr_df = pd.DataFrame(agg_pbr_list)
pbr_df.write_csv(output_file, separator="\t", include_header=True)
