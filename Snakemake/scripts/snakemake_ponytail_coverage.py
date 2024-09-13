## Libraries 

import polars as pl
import itertools
import os
from multiprocessing import Pool, cpu_count

# ---------- Setup ----------

sys.stdout = sys.stderr = open(snakemake.log[0], "wt")

mutations_file = snakemake.input.mutations[0]

output_file: str = snakemake.output.output_file

probes: str = snakemake.params.probes
filters: str = snakemake.params.filters
thresholds: str = snakemake.params.thresholds
samples: str = snakemake.params.samples
pbr_dir: str = snakemake.params.pbr_dir

print(f"PBR directory: {pbr_dir}")
print(f"Samples: {samples}")
print(f"Mutation file: {mutations_file}")
print(f"Probes: {probes}")
print(f"Filters: {filters}")
print(f"Thresholds: {thresholds}")
print(f"Output file: {output_file}")

## Read in mutations and get recurrent sites
mutations = pl.read_csv(mutations_file, separator="\t", has_header=True)
recurrent_sites = mutations.filter((pl.col("within_donor_recurrent") == True) | 
                                   (pl.col("within_cohort_recurrent") == True) |
                                   (pl.col("across_cohort_recurrent") == True)).select(["chr", "start", "stop"]).unique(subset=["chr", "start", "stop"])
recurrent_sites = recurrent_sites.with_columns(
    recurrent_status = pl.lit("recurrent")
)

# ---------- Functions ----------

def get_cov_prop(args): 
    
    sample, probe, applied_filter, threshold, pbr_dir, depth_cutoff = args
    
    filename = f"{pbr_dir}/{sample}.{probe}.{applied_filter}.{threshold}.nucleotide_count.annotated.bed.gz"
    # print(filename)
    if not os.path.exists(filename):
        print(f"File does not exist: {filename}")
        return []
    assert os.path.exists(filename)

    print(f"Processing ponytaill coverage data for: {filename}")
        
    ## Read in coverage dataset
    colnames = ["chr", "start", "stop", "ref_base", "depth", "a_count", "c_count", "g_count", "t_count", "n_count", "sample", "threshold"]
    coverage_df = pl.read_csv(filename, separator="\t", has_header=False, new_columns=colnames)
    
    ## Exclude recurrently mutated sites if all_filters is selected
    if (applied_filter == "all_filters"): 
        coverage_df = coverage_df.join(recurrent_sites, on=["chr", "start", "stop"], how="left", coalesce=True)
        coverage_df = coverage_df.with_columns(
            pl.col("recurrent_status").replace(None, "non-recurrent")
        )
        coverage_df = coverage_df.filter(pl.col("recurrent_status") == "non-recurrent")
            
    ## Count the total number of bases sequenced
    num_bases = coverage_df.shape[0]
    
    ## Count number of bases sequenced to depth > threshold
    num_bases_above_threshold = coverage_df.filter(pl.col("depth") >= depth_cutoff).shape[0]
    prop_bases_above_threshold = num_bases_above_threshold / num_bases
    
    sample_dict = {'sample': sample, 'depth_cutoff': depth_cutoff, 
                   'probe': probe, 'applied_filter': applied_filter, 'threshold': threshold,
                   'num_bases': num_bases, 'num_bases_above_threshold': num_bases_above_threshold,
                   'prop_bases_above_threshold': prop_bases_above_threshold}
    
    return sample_dict


# ----------------------------------------------------------------------------------------------------

if __name__ == "__main__": 
    
    ## Calculate proportion of bases covered by a certain depth value
    depth_list = list(range(0, 30000, 500))

    ## Get every sample, probe, filter, threshold combination
    print("Generating combinations for probes, filters, and thresholds")
    combinations = list(itertools.product(samples, probes, filters, thresholds, [pbr_dir], depth_list))
    # print(combinations)

    ## Perform multiproecssing
    print("Performing multiprocessing")
    with Pool(processes = 16) as pool: 
        results = pool.map(get_cov_prop, combinations)

    print(results)
    
    ## Flatten list of lists
    print("Flattening multiprocessed results")
    # flattened_results = [item for sublist in results for item in sublist]

    ## Write to output
    print("Writing results to output with polars")
    df = pl.DataFrame(results)
    df.write_csv(output_file, separator="\t", include_header=True)        

    # # Combine results
    # ponytail_df = pd.concat([pd.DataFrame(r) for r in results], ignore_index=True)
    
    # # Write to ouptut
    # ponytail_df.to_csv(output_file, sep="\t", index=False, header=True)