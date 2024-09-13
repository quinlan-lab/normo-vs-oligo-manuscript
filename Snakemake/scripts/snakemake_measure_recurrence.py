## Libraries
import polars as pl
import pandas as pd
from tqdm import tqdm
import time
tqdm.pandas()


# ---------- Setup ----------

sys.stdout = sys.stderr = open(snakemake.log[0], "wt")

mutation_file: str = snakemake.input.mutation_file[0]
output_file: str = snakemake.output.recurrence_file

"""
I have a mutation in a fertile sperm donor. These are the recurrnences I want to analyze: 

1. Within a donor --> check the donor's blood
2. Within a cohort --> check all other fertile donor's blood
3. Across cohort --> check all the subfertile and fertile donor's blood

"""

## Read in all mutations
all_mutations = pl.read_csv(mutation_file, has_header=True, separator="\t")
# all_mutations = all_mutations.sample(n = 500)

## Get all unique mutation calls
mutations_dedup = all_mutations.unique(["sample", "donor_id", "var_coord", "tissue", "fertility_status"])

num_unique_mutation_calls = mutations_dedup.shape[0]
print(f"Measuring mutation recurrence. There are {num_unique_mutation_calls} mutations to process.")

# ---------- Functions ----------

## Function to get the recurrence of a mutation
def get_var_recurrence(row, df):
    
    ## Get mutation-donor combination information    
    var_coord = row['var_coord'] ## Mutation coordinate
    donor_id = row['donor_id'] ## Donor with the mutation
    tissue = row['tissue'] ## Tissue of the mutation
    fertility_status = row['fertility_status']
    
    # ## Remove after testing
    # var_coord = "chr1:45332618-45332619-C>A"
    # donor_id = "1SC5"
    # tissue = "Sperm_DNA"
    # fertility_status = "Normozoospermic"
    # var_coord = "chr16:2091363-2091364-C>T"
    # donor_id = "3CM1"
    # tissue = "Blood_DNA"
    # fertility_status = "Normozoospermic"
    # df = all_mutations
    # var_coord = "chr11:108330187-108330188-G>A"
    # donor_id = "D130"
    # tissue = "Sperm_DNA"
    # fertility_status = "Normozoospermic"
    
    if tissue == "Sperm_DNA":
        recurrent_tissue = "Blood_DNA"
    else:
        recurrent_tissue = "Sperm_DNA"
        
        
    ## Get fertile donors for recurrent tissue (not tissue mutation is found in)
    fertile_nonclonal_donors = df.filter((pl.col("var_coord") == var_coord) & 
                                         (pl.col("tissue") == recurrent_tissue) & 
                                         (pl.col("fertility_status") == "Normozoospermic") & 
                                         (pl.col("alt_count_v2") <= 1)).get_column("donor_id").unique().to_list()  
    fertile_clonal_donors = df.filter((pl.col("var_coord") == var_coord) &
                                        (pl.col("tissue") == recurrent_tissue) & 
                                        (pl.col("fertility_status") == "Normozoospermic") & 
                                        (pl.col("alt_count_v2") > 1) & 
                                        (pl.col("AF_v2") < 0.05)).get_column("donor_id").unique().to_list()  
    fertile_gonosomal_donors = df.filter((pl.col("var_coord") == var_coord) &
                                         (pl.col("tissue") == recurrent_tissue) & 
                                         (pl.col("fertility_status") == "Normozoospermic") & 
                                         (pl.col("AF_v2") <= 0.2) & 
                                         (pl.col("AF_v2") >= 0.05)).get_column("donor_id").unique().to_list()  
    fertile_inherited_donors = df.filter((pl.col("var_coord") == var_coord) &
                                         (pl.col("tissue") == recurrent_tissue) & 
                                         (pl.col("fertility_status") == "Normozoospermic") & 
                                         (pl.col("AF_v2") >= 0.2)).get_column("donor_id").unique().to_list()
    
    ## Get subfertile blood donors for recurrent tissue (not tissue mutation is found in)
    subfertile_nonclonal_donors = df.filter((pl.col("var_coord") == var_coord) &
                                            (pl.col("tissue") == recurrent_tissue) & 
                                            (pl.col("fertility_status") == "Oligozoospermic") & 
                                            (pl.col("alt_count_v2") <= 1)).get_column("donor_id").unique().to_list()
    subfertile_clonal_donors = df.filter((pl.col("var_coord") == var_coord) &
                                         (pl.col("tissue") == recurrent_tissue) & 
                                         (pl.col("fertility_status") == "Oligozoospermic") & 
                                         (pl.col("alt_count_v2") > 1)).get_column("donor_id").unique().to_list()
    subfertile_gonosomal_donors = df.filter((pl.col("var_coord") == var_coord) &
                                            (pl.col("tissue") == recurrent_tissue) & 
                                            (pl.col("fertility_status") == "Oligozoospermic") & 
                                            (pl.col("AF_v2") <= 0.2)).get_column("donor_id").unique().to_list()
    subfertile_inherited_donors = df.filter((pl.col("var_coord") == var_coord) &
                                             (pl.col("tissue") == recurrent_tissue) & 
                                             (pl.col("fertility_status") == "Oligozoospermic") & 
                                             (pl.col("AF_v2") >= 0.2)).get_column("donor_id").unique().to_list()
        

    ## Get fertile donors for recurrent tissue (not tissue mutation is found in)
    # fertile_nonclonal_donors = df.loc[(df['var_coord'] == var_coord) & (df['tissue'] == recurrent_tissue) & (df['fertility_status'] == "Normozoospermic") & (df['alt_count_v2'] <= 1)]['donor_id'].unique()
    # fertile_clonal_donors = df.loc[(df['var_coord'] == var_coord) & (df['tissue'] == recurrent_tissue) & (df['fertility_status'] == "Normozoospermic") & (df['alt_count_v2'] > 1) & (df['AF_v2'] < 0.05)]['donor_id'].unique()
    # fertile_gonosomal_donors = df.loc[(df['var_coord'] == var_coord) & (df['tissue'] == recurrent_tissue) & (df['fertility_status'] == "Normozoospermic") & (df['AF_v2'] <= 0.2) & (df['AF_v2'] >= 0.05)]['donor_id'].unique()
    # fertile_inherited_donors = df.loc[(df['var_coord'] == var_coord) & (df['tissue'] == recurrent_tissue) & (df['fertility_status'] == "Normozoospermic") & (df['AF_v2'] >= 0.2)]['donor_id'].unique()

    ## Get subfertile blood donors for recurrent tissue (not tissue mutation is found in)
    # subfertile_nonclonal_donors = df.loc[(df['var_coord'] == var_coord) & (df['tissue'] == recurrent_tissue) & (df['fertility_status'] == "Oligozoospermic") & (df['alt_count_v2'] <= 1)]['donor_id'].unique()
    # subfertile_clonal_donors = df.loc[(df['var_coord'] == var_coord) & (df['tissue'] == recurrent_tissue) & (df['fertility_status'] == "Oligozoospermic") & (df['alt_count_v2'] > 1)]['donor_id'].unique()
    # subfertile_gonosomal_donors = df.loc[(df['var_coord'] == var_coord) & (df['tissue'] == recurrent_tissue) & (df['fertility_status'] == "Oligozoospermic") & (df['AF_v2'] <= 0.2)]['donor_id'].unique()
    # subfertile_inherited_donors = df.loc[(df['var_coord'] == var_coord) & (df['tissue'] == recurrent_tissue) & (df['fertility_status'] == "Oligozoospermic") & (df['AF_v2'] >= 0.2)]['donor_id'].unique()

    ## Gather fertile and subfertile donor lists (donors with the variant)
    fertile_donor_lists = [fertile_nonclonal_donors, fertile_clonal_donors, fertile_gonosomal_donors, fertile_inherited_donors]
    subfertile_donor_lists = [subfertile_nonclonal_donors, subfertile_clonal_donors, subfertile_gonosomal_donors, subfertile_inherited_donors]

    ## Get within donor, within cohort, and across cohort recurrence
    if fertility_status == "Normozoospermic":
        
        within_donor_recurrent = any(donor_id in donor_list for donor_list in fertile_donor_lists)
        across_cohort_recurrent = any(len(donor_list) > 0 for donor_list in subfertile_donor_lists)
        
        ## Remove donor_id from fertile donor lists
        new_fertile_donor_lists = [[donor for donor in donor_list if donor != donor_id] if len(donor_list) > 0 else donor_list for donor_list in fertile_donor_lists]
        new_lengths = [len(donor_list) for donor_list in new_fertile_donor_lists] 
        within_cohort_recurrent = any(length > 0 for length in new_lengths)

    if fertility_status == "Oligozoospermic": 
        within_donor_recurrent = any(donor_id in donor_list for donor_list in subfertile_donor_lists)
        across_cohort_recurrent = any(len(donor_list) > 0 for donor_list in fertile_donor_lists)
        
        ## Remove donor_id from subfertile donor lists
        new_subfertile_donor_lists = [[donor for donor in donor_list if donor != donor_id] if len(donor_list) > 0 else donor_list for donor_list in subfertile_donor_lists]
        new_lengths = [len(donor_list) for donor_list in new_subfertile_donor_lists] 
        within_cohort_recurrent = any(length > 0 for length in new_lengths)        

    var_info = {'var_coord': var_coord, 
             'donor_id': donor_id,
             'fertility_status': fertility_status,
             'tissue': tissue,
             'recurrent_tissue': recurrent_tissue,
             "within_donor_recurrent": within_donor_recurrent, 
             "within_cohort_recurrent": within_cohort_recurrent, 
             "across_cohort_recurrent": across_cohort_recurrent}

    return var_info
            
# ---------- Measure recurrence for mutations in each donor ----------
 
## Iterate through unique mutation-donor combination to get recurrence
## The comma after all_mutations_subset is necessary because the args parameter of the apply() function expects a tuple. 
## In Python, a tuple with a single element is defined by including a trailing comma.
res_list = []
counter = 0
start_time = time.time()
for row in mutations_dedup.rows(named=True): 
    # print(row)
    counter += 1
    if counter % 100 == 0:
        stop_time = time.time()
        elapsed_time = stop_time - start_time
        print(f"Processed {counter} mutations. Total elapsed time: {elapsed_time} seconds.")
    res = get_var_recurrence(row, all_mutations)
    # print(res)
    res_list.append(res)
    
recurrence_df = pl.DataFrame(res_list)    
print(recurrence_df)

merged_mutations = all_mutations.join(recurrence_df, on = ["var_coord", "donor_id", "tissue", "fertility_status"], how="left", coalesce=True)
print(merged_mutations)

print("Writing to output file")
merged_mutations.write_csv(output_file, separator="\t", include_header=True)

# var_coord_recurrence = mutations_dedup.progress_apply(get_var_recurrence, axis=1, args=(all_mutations,))    
# var_coord_recurrence_df = pd.DataFrame(var_coord_recurrence.tolist())

## Merge recurrence information with original mutations
# merged_mutations = all_mutations.merge(var_coord_recurrence_df, on=['var_coord', 'donor_id', 'tissue', 'fertility_status'], how='left')

## Write to file
# merged_mutations.to_csv(output_file, sep="\t", index=False, header=True)
 