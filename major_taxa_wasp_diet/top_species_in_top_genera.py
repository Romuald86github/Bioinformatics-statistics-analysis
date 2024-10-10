import pandas as pd
import matplotlib.pyplot as plt

species_counts_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/wasp_gut_results_with_origins.tsv'

# Load the origin_df (species counts and taxonomy)
origin_df = pd.read_csv(species_counts_path, sep='\t')

# (1) Count the number of species per genus, excluding 'sp.' and 'nr.'
genus_species_counts = origin_df[~origin_df['species'].str.contains('sp.|nr.', case=False)].groupby(['genus', 'species']).size().reset_index(name='count')

# (2) Select the top 10 genera
top_10_genera = genus_species_counts.groupby('genus')['count'].sum().nlargest(10).index

# (3) Prepare the stacked bar plot
fig, ax = plt.subplots(figsize=(12, 8))
x = np.arange(len(top_10_genera))  # the label locations
width = 0.8 / 5  # the width of the bars

for i, genus in enumerate(top_10_genera):
    genus_df = genus_species_counts[genus_species_counts['genus'] == genus]
    genus_df = genus_df.sort_values('count', ascending=False).head(5)
    
    species = genus_df['species'].tolist()
    counts = genus_df['count'].tolist()
    
    for j, count in enumerate(counts):
        ax.bar(x[i] + j * width, count, width, label=species[j])

# (4) Customize the plot
ax.set_xlabel('Genus')
ax.set_ylabel('Number of Species')
ax.set_title('Top 5 Species per Top 10 Genera (Excluding "sp." and "nr.")')
ax.set_xticks(x)
ax.set_xticklabels(top_10_genera, rotation=90)
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()

# Save the plot
plot_output_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/top_5_species_per_genus_true.png'
plt.savefig(plot_output_path)

plt.show()

#print(f"Stacked bar plot saved to: {plot_output_path}")