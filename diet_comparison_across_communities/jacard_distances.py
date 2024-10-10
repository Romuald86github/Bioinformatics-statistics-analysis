import pandas as pd
import numpy as np
from skbio.diversity import beta_diversity
from skbio import DistanceMatrix
import matplotlib.pyplot as plt
import seaborn as sns

# Load data (species counts and taxonomy)
species_counts_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/wasp_gut_results_with_origins.tsv'
metadata_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/metadata.tsv'

# Load the origin_df
origin_df = pd.read_csv(species_counts_path, sep='\t')

# Drop the ASV column
origin_df = origin_df.drop('ASV', axis=1)

# Extract the taxonomy information
taxonomy_cols = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies']
taxonomy_df = origin_df[taxonomy_cols]

# Extract the species counts (samples)
samples_df = origin_df.drop(taxonomy_cols, axis=1)

# Load the metadata
metadata = pd.read_csv(metadata_path, sep='\t', index_col='Larvae')

# Assign locations based on the nest information
metadata['Location'] = metadata['Nest'].apply(lambda x: 'Apir' if x in list('ABCDEFGHIJ') else
                                             'Ihala' if x in list('KLMNOPQRST') else
                                             'University of Agriculture' if x in list('U') else
                                             'Nsukka' if x in list('WX') else
                                             'Unknown')

# (1) Calculate beta diversity (Jaccard distance) between all larvae
beta_diversity_all = beta_diversity('jaccard', samples_df.T.values, ids=samples_df.columns)
beta_diversity_all_df = pd.DataFrame(beta_diversity_all.data, index=beta_diversity_all.ids, columns=beta_diversity_all.ids)

# Visualize beta diversity between all larvae
plt.figure(figsize=(10, 8))
sns.clustermap(beta_diversity_all_df, cmap='viridis', row_cluster=True, col_cluster=True)
plt.title("Beta Diversity/similarity (Jaccard Distance) Between All Larvae")

# Save the plot
plot_output_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/beta_diversity_clustermap.png'
plt.savefig(plot_output_path)

plt.show()

print(f"Beta diversity clustermap saved to: {plot_output_path}")