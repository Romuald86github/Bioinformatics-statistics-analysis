import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from skbio.diversity import alpha_diversity
from sklearn.decomposition import PCA

# Load the wasp gut data with larvae origins and metadata
data_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/wasp_gut_results_with_origins.tsv'
data_df = pd.read_csv(data_path, sep='\t')

# Load metadata containing nests and locations (larvae IDs are set as index)
metadata_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/metadata.tsv'
metadata_df = pd.read_csv(metadata_path, sep='\t', index_col='Larvae')

# Remove the ASV and taxonomy columns (keeping only the larvae abundance data)
taxonomy_cols = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies']
sample_cols = [col for col in data_df.columns if col not in ['ASV'] + taxonomy_cols]

# Transpose the larvae abundance data so that larvae IDs are rows
abundance_df = data_df[sample_cols].T  # Transpose
abundance_df.index.name = 'Larvae'  # Set index name to match metadata

# Merge the transposed abundance data with the metadata using larvae IDs
merged_df = pd.merge(abundance_df, metadata_df, left_index=True, right_index=True)

# Check the number of metadata columns (to exclude them later)
num_metadata_cols = len(metadata_df.columns)

# Extract the abundance matrix (selecting by column positions rather than names)
abundance_matrix = merged_df.iloc[:, :-num_metadata_cols].values

# Assign nest and location information for grouping
nests = merged_df['Nest'].unique()
locations = merged_df['Location'].unique()

# Verify the result
print("Abundance matrix shape:", abundance_matrix.shape)
print("Nests:", nests)
print("Locations:", locations)




# PCA
pca = PCA(n_components=2)
pca_results = pca.fit_transform(abundance_matrix)

# Plotting PCA to visualize prey preferences
plt.figure(figsize=(8, 6))
plt.scatter(pca_results[:, 0], pca_results[:, 1], c='blue')
plt.title('Dietary Preferences (PCA)')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.savefig('/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/pca_plot.png')
plt.show()

# Calculate and visualize Levin's Niche Breadth
def levins_niche_breadth(abundances):
    p = abundances / abundances.sum()
    return 1 / np.sum(p ** 2)

niche_breadths = [levins_niche_breadth(sample) for sample in abundance_matrix]
plt.figure(figsize=(8, 6))
plt.hist(niche_breadths, bins=10)
plt.title('Levin\'s Niche Breadth Distribution')
plt.xlabel('Niche Breadth')
plt.ylabel('Frequency')
plt.savefig('/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/niche_breadth_histogram.png')
plt.show()