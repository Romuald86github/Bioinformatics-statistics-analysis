import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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
asv_counts_st = origin_df.groupby(genus_column)[larvae_st].sum().sum(axis=1) / num_larvae_st

# Group by genus and calculate the sum of ASV counts for all other nests (divided by number of larvae)
asv_counts_others = origin_df.groupby(genus_column)[larvae_others].sum().sum(axis=1) / num_larvae_others

# Step 7: Combine the results into a single DataFrame for comparison
genus_comparison = pd.DataFrame({
    'Nests_S_T': asv_counts_st,
    'Other_Nests': asv_counts_others
}).fillna(0)

# Step 8: Find the top 10 genera based on total ASV counts
genus_comparison['Total'] = genus_comparison['Nests_S_T'] + genus_comparison['Other_Nests']
top_10_genera = genus_comparison.nlargest(10, 'Total').index

# Step 9: Filter the DataFrame to include only the top 10 genera
genus_comparison_top10 = genus_comparison.loc[top_10_genera]

# Step 10: Create a bar plot comparing average ASV counts per larvae between Nests S and T and other nests
genus_comparison_top10.reset_index(inplace=True)

plt.figure(figsize=(12, 8))
sns.barplot(x='genus', y='value', hue='variable', 
            data=pd.melt(genus_comparison_top10, id_vars='genus', value_vars=['Nests_S_T', 'Other_Nests']),
            palette="coolwarm")

# Add labels and title
plt.title('Top 10 Genera for Nests S and T Compared to Other Nests (Average per Larvae)', fontsize=16)
plt.xlabel('Genus', fontsize=14)
plt.ylabel('Average ASV Count per Larvae', fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

# Save the plot
plot_output_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/top_5_species_proportions_per_genus_st_vs_other.png'
plt.savefig(plot_output_path)


# Show the plot
plt.show()
