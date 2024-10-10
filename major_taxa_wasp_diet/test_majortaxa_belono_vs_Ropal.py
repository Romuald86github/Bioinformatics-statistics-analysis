import pandas as pd
from scipy import stats

# Load the wasp gut data with larvae origins and metadata
data_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/wasp_gut_results_with_origins.tsv'
origin_df = pd.read_csv(data_path, sep='\t')

# Step 1: Define the nests for S and T
nests_st = [col for col in origin_df.columns if col.startswith('S') or col.startswith('T')]

# Step 2: Define the other nests
other_nests = [col for col in origin_df.columns if col not in nests_st and col not in origin_df.columns[-8:] and col != 'ASV']

# Step 3: Identify larvae for Nests S and T
larvae_st = [col for col in nests_st if col in origin_df.columns]

# Step 4: Identify larvae for other nests
larvae_others = [col for col in other_nests if col in origin_df.columns]

# Step 5: Count the number of larvae in each group
num_larvae_st = len(larvae_st)
num_larvae_others = len(larvae_others)

# Step 6: Group by genus and calculate the sum of ASV counts for Nests S and T (divided by number of larvae)
genus_column = 'genus'  
asv_counts_st = origin_df.groupby(genus_column)[larvae_st].sum().sum(axis=1) / num_larvae_st

# Group by genus and calculate the sum of ASV counts for all other nests (divided by number of larvae)
asv_counts_others = origin_df.groupby(genus_column)[larvae_others].sum().sum(axis=1) / num_larvae_others

# Step 7: Perform the Mann-Whitney U Test for all genera
u_stat_all, p_value_all = stats.mannwhitneyu(asv_counts_st, asv_counts_others)

# Output the results for all genera
print(f"All Genera Mann-Whitney U Test: U-statistic = {u_stat_all}, p-value = {p_value_all}")
