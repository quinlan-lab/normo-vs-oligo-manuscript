## Load libraries
import pandas as pd
import numpy as np

# ---------- Setup ----------

sys.stdout = sys.stderr = open(snakemake.log[0], "wt")

output_file: str = snakemake.output.combined_mutations

metadata_file = snakemake.input.metadata_file[0]
mutation_file_dir: str = snakemake.params.mutation_file_dir
samples: str = snakemake.params.samples
# samples = [f"{sample}.1" for sample in samples]

print(f"Samples: {samples}")
print(f"Mutation file dir: {mutation_file_dir}")
print(f"Output file: {output_file}")

# ---------- Combine mutation files ----------

dfs = []  # list to store DataFrames

## Iterate through each sample
for sample in samples: 
    print(f"Processing mutation file for sample: {sample}")
    sample_file = f"{mutation_file_dir}/{sample}.pysam_processed.txt"
    print(sample_file)
    sample_mutation_df = pd.read_csv(sample_file, sep="\t", header=0)
    sample_mutation_df['sample'] = sample
    dfs.append(sample_mutation_df)  # add DataFrame to list

# Concatenate all DataFrames in the list
combined_df = pd.concat(dfs, axis=0)

## Define conditions for allele frequency groups
conditions = [
    (combined_df['AF_v2'] >= 0.2),
    (combined_df['AF_v2'] < 0.2) & (combined_df['AF_v2'] >= 0.05),
    (combined_df['AF_v2'] < 0.05) & (combined_df['alt_count_v2'] > 1),
    (combined_df['AF_v2'] < 0.05) & (combined_df['alt_count_v2'] == 1)
]
outputs = ['inherited', 'gonosomal', 'clonal', 'nonclonal']
combined_df['af_group'] = np.select(conditions, outputs, 'other')

print(combined_df.head())

# ---------- Add clinical/coverage information to the mutation file ----------

## merge metadata df with mutation df
print("Merging metadata with mutation data")
metadata_df = pd.read_csv(metadata_file, sep="\t", header=0)
combined_df = pd.merge(combined_df, metadata_df[["sample", "donor_id", "tissue", "fertility_status", "timepoint"]], on=["sample"], how="inner")

combined_df.to_csv(output_file, sep="\t", index=False, header=True)
