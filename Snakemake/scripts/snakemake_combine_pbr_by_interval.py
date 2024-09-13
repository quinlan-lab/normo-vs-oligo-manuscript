## Load libraries

import polars as pl
from multiprocessing import Pool, cpu_count

# ---------- Setup ----------

sys.stdout = sys.stderr = open(snakemake.log[0], "wt")

mutations_file = snakemake.input.mutations[0]
probe_bed_file: str = snakemake.input.probe_bed_file

output_file: str = snakemake.output.output_file
samples: str = snakemake.params.samples
pbr_dir: str = snakemake.params.pbr_dir

print(f"Probe bed file: {probe_bed_file}")
print(f"Mutation file: {mutations_file}")
print(f"Samples: {samples}")
print(f"PBR directory: {pbr_dir}")
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

## Function to get interval coverage
def get_interval_coverage(sample_pbr_df, chrom, start, stop):
    coverage = sample_pbr_df.filter((pl.col("chr") == chrom) & 
                                    (pl.col("start") >= start) &
                                    (pl.col("stop") <= stop)).select("depth").sum()
    return coverage

## Function to process coverage per interval in a sample
def process_sample(args):
        
        sample, probe_bed_file, pbr_dir, desired_filter = args
        print(f"Processing sample: {sample}... Reading in the probe bed file {probe_bed_file[0]} and using the PBR directory: {pbr_dir[0]}")
        
        ## Read in the probe bed file
        probe_cols = ['chr', 'start', 'stop', 'gene']
        probe_df = pl.read_csv(str(probe_bed_file[0]), separator="\t", has_header=False, new_columns=probe_cols)
        
        ## Specify if pre/post filter pbr file should be used
        if (desired_filter == "postfilter") | (desired_filter == "postfilter_remove_recurrence"): 
            sample_file = f"{pbr_dir[0]}/{sample}.subset_intervals.all_filters.baseline.nucleotide_count.annotated.bed.gz"
        if desired_filter == "prefilter": 
            sample_file = f"{pbr_dir[0]}/{sample}.all_intervals.none_filters.baseline.nucleotide_count.annotated.bed.gz"
            
        ## Read in sample's pbr file
        pbr_cols = ['chr', 'start', 'stop', 'ref_allele', 'depth', 'a_count', 'c_count', 'g_count', 't_count', 'n_count', 'sample', 'threshold']
        sample_pbr_df = pl.read_csv(sample_file, separator="\t", has_header=False, new_columns=pbr_cols)
        
        ## Exclude recurrently mutated sites if all_filters is selected
        if desired_filter == "postfilter_remove_recurrence": 
            
            sample_pbr_df = sample_pbr_df.join(recurrent_sites, on=["chr", "start", "stop"], how="left", coalesce=True)
            sample_pbr_df = sample_pbr_df.with_columns(
                pl.col("recurrent_status").replace(None, "non-recurrent")
            )
            sample_pbr_df = sample_pbr_df.filter(pl.col("recurrent_status") == "non-recurrent")
        
        sample_list = []
        
        ## Get the coverage for each interval
        for row in probe_df.rows(named=True):
            # if index > 15: 
            #     continue
            # print(index)
            chrom = row['chr']
            start = row['start']
            stop = row['stop']
            gene = row['gene']
            
            ## Get the coverage for the interval
            coverage = get_interval_coverage(sample_pbr_df=sample_pbr_df, chrom=chrom, start=start, stop=stop)
            
            sample_list.append({'sample': sample, 
                                'chr': chrom,
                                'start': start,
                                'stop': stop,
                                'gene': gene,
                                'coverage': coverage, 
                                'filter': desired_filter})
        return sample_list
            
# ---------- PERFORM MULTIPROCESSING ----------

## Set up arguments for multiprocessing
prefilter_args = [(sample, [probe_bed_file], [pbr_dir], "prefilter") for sample in samples]
postfilter_args = [(sample, [probe_bed_file], [pbr_dir], "postfilter") for sample in samples]
postfilter_remove_recurrence_args = [(sample, [probe_bed_file], [pbr_dir], "postfilter_remove_recurrence") for sample in samples]

with Pool(processes=8) as pool:
    prefiltered_results = pool.map(process_sample, prefilter_args)
    postfiltered_results = pool.map(process_sample, postfilter_args)
    postfiltered_remove_recurrence_results = pool.map(process_sample, postfilter_remove_recurrence_args)
    
    
## Combine results and write to output
print("writing to output file...")

prefiltered_flattened_results = [item for sublist in prefiltered_results for item in sublist]
postfiltered_flattened_results = [item for sublist in postfiltered_results for item in sublist]
postfiltered_remove_recurrence_flattened_results = [item for sublist in postfiltered_remove_recurrence_results for item in sublist]

prefiltered_interval_cov_df = pl.DataFrame(prefiltered_flattened_results)
postfiltered_interval_cov_df = pl.DataFrame(postfiltered_flattened_results)
postfiltered_remove_recurrence_interval_cov_df = pl.DataFrame(postfiltered_remove_recurrence_flattened_results)

interval_cov_df = pl.concat([prefiltered_interval_cov_df, postfiltered_interval_cov_df, postfiltered_remove_recurrence_interval_cov_df])

interval_cov_df.write_csv(output_file, separator="\t", include_header=True)
