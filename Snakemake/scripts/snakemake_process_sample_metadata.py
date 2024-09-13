import os
import pandas as pd

sys.stdout = sys.stderr = open(snakemake.log[0], "wt")

experiment_metadata_file = snakemake.input.experiment_metadata_file
clinical_metadata_file = snakemake.input.clinical_metadata_file
output_file = snakemake.output.sample_metadata
samples = snakemake.params.samples
grouping_dir = snakemake.params.grouping_dir
duplex_yield_dir = snakemake.params.duplex_yield_dir

## Read in experiment metadata
experiment_df = pd.read_csv(experiment_metadata_file, sep="\t", header=0)

## Read in clinical metadata
clinical_df = pd.read_csv(clinical_metadata_file, sep="\t", header=0)

## merge clinical data with experiment data
merged_df = pd.merge(experiment_df, clinical_df, on=["donor_id", "tissue", "age", "fertility_status", "timepoint"], how="inner")

## Get the PTFS and total reads per sample
sample_list = []
for sample in samples: 
    grouping_file = os.path.join(grouping_dir, sample + ".grouped-family-sizes.txt")
    yield_file = os.path.join(duplex_yield_dir, sample + ".duplex_seq_metrics.duplex_yield_metrics.txt")
    
    grouping_df = pd.read_csv(grouping_file, sep="\t", header=0)
    yield_df = pd.read_csv(yield_file, sep="\t", header=0)
    
    ## Get the sample ptfs
    grouping_df['total_reads_per_family'] = grouping_df['family_size'] * grouping_df['count']
    grouping_df = grouping_df.loc[grouping_df['family_size'] > 1]
    max_row = grouping_df.loc[grouping_df['total_reads_per_family'].idxmax()]
    max_family_size = max_row['family_size']
    
    ## Get the total reads in the library
    total_read_count = pd.to_numeric(yield_df.loc[yield_df['fraction'] == 1, "read_pairs"]).values[0] * 2
    
    sample_dict = {'sample': sample, 'ptfs': max_family_size, 'total_reads': total_read_count}
    sample_list.append(sample_dict)
    
## Convert sample library information to a dataframe
sample_df = pd.DataFrame(sample_list)

## Merge the sample library information with the merged_df
merged_df = pd.merge(merged_df, sample_df, on="sample", how="inner")    

## Make sure all of samples are in the final_sample_set
# assert set(samples) == set(final_sample_set), "Samples in the sample list do not match the samples in the final sample set."
final_sample_set = merged_df['sample'].tolist()
assert set(samples).issubset(set(final_sample_set)), "Samples in the sample list do not match the samples in the final sample set."

## Write to output_fle
merged_df.to_csv(output_file, sep="\t", index=False, header=True)


