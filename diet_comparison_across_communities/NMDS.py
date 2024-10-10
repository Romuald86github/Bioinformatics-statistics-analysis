import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from skbio.diversity import alpha_diversity
from skbio.diversity import beta_diversity
from sklearn.manifold import MDS
from matplotlib.patches import Ellipse

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













# Compute Bray-Curtis dissimilarity matrix
bray_curtis_dm = beta_diversity('braycurtis', abundance_matrix)

# Convert DistanceMatrix to a NumPy array
bray_curtis_matrix = bray_curtis_dm.data

# Handle NaN values by removing rows and columns with NaNs
nan_mask = np.isnan(bray_curtis_matrix).any(axis=1)
filtered_indices = ~nan_mask
filtered_matrix = bray_curtis_matrix[filtered_indices][:, filtered_indices]

# Perform NMDS
nmds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
nmds_results = nmds.fit_transform(filtered_matrix)

# NMDS coordinates
nmds_coords = pd.DataFrame(nmds_results, columns=['NMDS1', 'NMDS2'], index=merged_df.index[filtered_indices])

# Unique markers for each location
markers = ['o', 's', 'D', '^', 'v', '<', '>', 'p', 'h', 'H', 'X', '*']  # List of different shapes (markers)

# Plot NMDS results with ellipses
plt.figure(figsize=(12, 10))
unique_locations = np.unique(locations)
colors = plt.cm.get_cmap('tab20', len(unique_locations))

for i, location in enumerate(unique_locations):
    loc_indices = np.where(np.isin(merged_df.index[filtered_indices], merged_df.index[locations == location]))[0]
    coords = nmds_coords.iloc[loc_indices]
    
    # Plot points with unique colors and shapes
    plt.scatter(coords['NMDS1'], coords['NMDS2'], label=location, color=colors(i), marker=markers[i % len(markers)])
    
    # Compute mean and covariance
    mean = coords.mean()
    cov = np.cov(coords.T)
    
    # Calculate the 2 standard deviation ellipse to cover more data
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    radii = np.sqrt(eigenvalues) * 2  # 2 standard deviations
    angle = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))
    
    ellipse = Ellipse(mean, width=radii[0]*2, height=radii[1]*2, angle=angle, edgecolor=colors(i), facecolor='none', linestyle='--')
    plt.gca().add_patch(ellipse)

plt.xlabel('NMDS1')
plt.ylabel('NMDS2')
plt.title('NMDS Plot by Location with Different Shapes and Colors')
plt.legend()
plt.savefig('/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/nmds_plot_with_ellipses_and_shapes.png')
plt.show()
