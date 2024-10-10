from skbio.diversity import beta_diversity
from skbio.stats.distance import permanova
from skbio import DistanceMatrix
import pandas as pd
from skbio.stats.ordination import pcoa

# Load the origin data (which contains sample origins, nests, and larvae)
origin_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/wasp_gut_results_with_origins.tsv'
origin_df = pd.read_csv(origin_path, sep='\t')

# (1) Remove the ASV column (first column)
origin_df = origin_df.drop(columns=['ASV'])

# (2) Extract the larvae samples (columns excluding the last 8 taxonomy columns)
taxonomy_cols = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies']
sample_cols = [col for col in origin_df.columns if col not in taxonomy_cols]

# Extract the species matrix (ASVs as rows, sample origins as columns)
species_matrix = origin_df[sample_cols]



species_matrix = species_matrix.join(taxonomy_df['species'])
species_matrix = species_matrix.groupby('species').sum()




# Calculate Bray-Curtis dissimilarity between all larvae
bray_curtis = beta_diversity('braycurtis', species_matrix.T, ids=species_matrix.columns)

# Add the dissimilarity matrix to a dataframe for better understanding
bray_curtis_df = pd.DataFrame(bray_curtis.data, index=bray_curtis.ids, columns=bray_curtis.ids)
bray_curtis_df.head()  # Visualize the Bray-Curtis dissimilarity matrix







# Convert the Bray-Curtis distance matrix to a DataFrame for easier manipulation
bray_curtis_df = pd.DataFrame(bray_curtis.data, index=bray_curtis.ids, columns=bray_curtis.ids)

# Replace inf and -inf with NaNs
bray_curtis_df_clean = bray_curtis_df.replace([np.inf, -np.inf], np.nan)

# Identify the rows/columns that contain NaN values
nan_samples = bray_curtis_df_clean.index[bray_curtis_df_clean.isna().any(axis=1) | bray_curtis_df_clean.isna().any(axis=0)]

# Drop these rows and columns symmetrically to maintain a square matrix
bray_curtis_clean = bray_curtis_df_clean.drop(index=nan_samples, columns=nan_samples)

# Ensure the matrix is now square
assert bray_curtis_clean.shape[0] == bray_curtis_clean.shape[1], "Matrix is not square after cleaning!"

# Convert back to a DistanceMatrix object
bray_curtis_clean_dm = DistanceMatrix(bray_curtis_clean.values, ids=bray_curtis_clean.index)

# Perform PCoA on the cleaned Bray-Curtis matrix
pcoa_results = pcoa(bray_curtis_clean_dm)
pcoa_df = pcoa_results.samples

# Add metadata (nest and location) for coloring in the plot
pcoa_df['Nest'] = metadata['Nest'].reindex(pcoa_df.index)
pcoa_df['Location'] = metadata['Location'].reindex(pcoa_df.index)

# Plot PCoA results, coloring by Location

# Access the first two principal components (PC1 and PC2) using .iloc
plt.figure(figsize=(10, 8))

# Scatterplot for PC1 vs PC2, coloring by location and styling by nest
sns.scatterplot(x=pcoa_df.iloc[:, 0], y=pcoa_df.iloc[:, 1], hue=pcoa_df['Location'], style=pcoa_df['Nest'])

# Add titles and labels for the axes
plt.title("PCoA of Beta Diversity (Bray-Curtis) Between Locations")
plt.xlabel(f"PC1 ({pcoa_results.proportion_explained[0]:.2%} variance)")
plt.ylabel(f"PC2 ({pcoa_results.proportion_explained[1]:.2%} variance)")

# Position the legend outside the plot to the right side
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="Legend")

# Adjust layout to make space for the legend
plt.tight_layout(rect=[0, 0, 0.85, 1])

# Display the plot
plt.show()
